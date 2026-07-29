# packages/HimalayaUI/src/ingest.jl
#
# scan_and_group!: the transactional ingest orchestrator (spec §9.1).
#
# Discipline cloned from _reingest_inner! (cli.jl) + analyze_exposure! (pipeline.jl):
#   - Insert-only for loads/samples/exposures (no clobber on rescan).
#   - The exposure dedup key is (experiment_id, filename) — also enforced by the
#     `exposures_unique_filename` UNIQUE index (db.jl migrate_exposures_unique_filename!).
#   - Geometry is FILL-ONLY: it writes a field only when still unset ('default');
#     any established source (setup/prp/computed, or 'user') is left untouched, so
#     geometry derived once (at the funnel preview) is never re-derived (spec §4).
#   - Every DB write is wrapped in `_DB_WRITE_LOCK` (server.jl) + a SQLite.transaction
#     so a mid-scan failure rolls back the whole insert batch (no partial write).
#   - analyze_exposure! is called OUTSIDE the write transaction (same contract as
#     cli_init_with_db!): an analysis crash must not roll back the ingested rows.
#
# Note: `_reingest_inner!` is not called through; scan_and_group! is net-new and
# clones the discipline because _reingest_inner! requires a manifest (returns
# :no_manifest early) and clobbers sample fields, neither of which fits the
# manifest-free, never-clobber ingest path.

"""
    _find_setup_files(analysis_dir) -> Vector{String}

Locate the experiment's `setup_info_*.txt` (AgBe calibration: beam center +
flight distance). The configured `analysis_dir` is frequently a leaf such as
`.../analysis/automatic_analysis/tot_files`, while `setup_info` lives in the
analysis root, so search `analysis_dir` and walk up its ancestors, returning the
matches from the first directory that has any.

# ponytail: cap the ascent at 3 ancestors — covers `analysis/automatic_analysis/
# tot_files`. Widen the bound if a deeper analysis layout ever appears.
"""
function _find_setup_files(analysis_dir::AbstractString)
    dir = analysis_dir
    for _ in 1:4   # analysis_dir + 3 ancestors
        if isdir(dir)
            hits = String[joinpath(dir, f) for f in readdir(dir)
                          if startswith(f, "setup_info_") && endswith(f, ".txt")]
            isempty(hits) || return hits
        end
        parent = dirname(dir)
        parent == dir && break   # reached the filesystem root
        dir = parent
    end
    return String[]
end

"""
    scan_and_group!(db, experiment_id; analyze=true, on_progress=nothing)

Full ingest of a beamtime directory into `db` under `experiment_id`. The scan
root is the experiment's own `data_dir`/`analysis_dir` (resolved from its row).

Steps:
  1. Resolve `data_dir` and `analysis_dir` from the `experiments` row.
  2. Scan: enumerate TIF+PRP+DAT triplets with `scan_directory`.
  3. Geometry: FILL-ONLY — derive from the sampled PRPs + setup_info, but write
     only fields still unset; never overwrite geometry already established (e.g.
     committed from the funnel preview, which derived it with the confirmed setup
     file). Derived once, not re-derived each scan.
  4. Group: `group_into_samples` → loads/samples/slots.
  5. Persist: insert `loads`, `samples`, `exposures` rows inside ONE write
     transaction (insert-only — dedup loads by (experiment_id, load_index),
     samples by (experiment_id, load_id, slot_index), exposures by
     (experiment_id, filename) — so a clean rescan is a true no-op).
  6. Analyze: run `analyze_exposure!` for every newly inserted exposure, OUTSIDE
     the write transaction (same contract as `cli_init_with_db!`). Skipped when
     `analyze=false` (test-only).

The insert-only discipline makes every scan additive: a clean rescan is a no-op,
new files extend the existing loads/samples/exposures. Returns a NamedTuple summary
`(status, added_loads, added_samples, added_exposures)`.
"""
function scan_and_group!(
    db           ::SQLite.DB,
    experiment_id::Int;
    analyze      ::Bool   = true,
    on_progress  ::Union{Function,Nothing} = nothing,
)
    # Resolve data_dir, analysis_dir, and per-experiment pattern overrides from the row.
    exp_row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT data_dir, analysis_dir, image_pattern, metadata_pattern, integration_pattern FROM experiments WHERE id = ?", [experiment_id])))
    data_dir     = String(exp_row.data_dir)
    analysis_dir = String(exp_row.analysis_dir)

    # Read file patterns from DB columns (set via PATCH /api/experiments/:id); fall back to
    # legacy {name}.<ext> defaults when the column is NULL.
    tif_pattern = coalesce(exp_row.image_pattern,       "{name}.tif")
    prp_pattern = coalesce(exp_row.metadata_pattern,    "{name}.prp")
    dat_pattern = coalesce(exp_row.integration_pattern, "{name}.dat")

    # -----------------------------------------------------------------------
    # 1. Scan: enumerate TIF+PRP+DAT triplets
    # -----------------------------------------------------------------------
    metas = scan_directory(data_dir, analysis_dir;
        tif_pattern = tif_pattern,
        prp_pattern = prp_pattern,
        dat_pattern = dat_pattern)

    isempty(metas) && return (status = :empty, added_loads = 0, added_samples = 0, added_exposures = 0)

    # Publish the DENOMINATOR as soon as discovery knows it. Everything below
    # (geometry derive, grouping, the persist txn) is silent, and the analyze
    # loop used to be the first and only thing that ever reported a total — so
    # the UI sat on "0 / ~0" through the entire slow half of the scan and only
    # then grew a scale. One tick here means the bar has an honest denominator
    # from the moment the file list is known.
    n_found = length(metas)
    on_progress === nothing || on_progress(0, n_found)

    # -----------------------------------------------------------------------
    # 2. Geometry: derive + write to experiments row (never-clobber human fields)
    # -----------------------------------------------------------------------
    prp_paths   = String[m.prp_path for m in metas if m.prp_path !== nothing]
    setup_files = _find_setup_files(analysis_dir)

    geo, _disc = derive_geometry(prp_paths, setup_files)

    # FILL-ONLY: write only geometry fields still unset ('default' source). Fields
    # already established at create (committed from the funnel preview, which
    # derived them with the confirmed setup file) or manually edited ('user') are
    # left untouched — geometry is derived once, not re-derived on every scan.
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            _fill_unset_geometry!(db, experiment_id, geo)
        end
    end

    # -----------------------------------------------------------------------
    # 3. Group: loads / samples / slots
    # -----------------------------------------------------------------------
    grouping = group_into_samples(metas)

    # -----------------------------------------------------------------------
    # 4. Persist: loads → samples → exposures (all inside ONE write transaction
    #    so a mid-scan failure rolls back the whole batch — no partial write).
    # -----------------------------------------------------------------------
    new_exposure_ids = Int[]

    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            # Existing exposure filenames for this experiment (dedup key).
            existing_filenames = Set{String}(String(r.filename)
                for r in Tables.rowtable(DBInterface.execute(db,
                    "SELECT filename FROM exposures WHERE experiment_id = ? AND filename IS NOT NULL",
                    [experiment_id])))

            for gl in grouping.loads
                # Dedup the load row by (experiment_id, load_index). The loads table
                # has only a NON-unique index (loads_experiment_idx), so an
                # unconditional INSERT would mint a fresh load_id on every rescan —
                # that re-keys the sample dedup below onto a NEW load_id, finds
                # nothing, and creates duplicate samples + an orphan load row.
                # Mirror the sample/exposure insert-only discipline: reuse the
                # existing load_id when present; insert only when absent. This makes
                # a clean rescan a true no-op.
                existing_load = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM loads WHERE experiment_id = ? AND load_index = ?",
                    [experiment_id, gl.load_index]))
                if isempty(existing_load)
                    load_id = Int(DBInterface.lastrowid(DBInterface.execute(db, """
                        INSERT INTO loads (experiment_id, load_index, frame_count,
                                           start_time, end_time, session_id)
                        VALUES (?, ?, ?, ?, ?, ?)
                    """, [experiment_id, gl.load_index, gl.frame_count,
                          ismissing(gl.start_time) ? nothing : Dates.format(gl.start_time, "yyyy-mm-ddTHH:MM:SS"),
                          ismissing(gl.end_time)   ? nothing : Dates.format(gl.end_time,   "yyyy-mm-ddTHH:MM:SS"),
                          gl.session_index])))
                else
                    load_id = Int(first(existing_load).id)
                end

                for gs in gl.samples
                    # Sample dedup by (experiment_id, load_id, slot_index). Reuse the
                    # existing row on rescan rather than re-inserting (also avoids
                    # tripping the samples_unique_name UNIQUE index).
                    existing_sample = Tables.rowtable(DBInterface.execute(db,
                        "SELECT id FROM samples WHERE experiment_id = ? AND load_id = ? AND slot_index = ?",
                        [experiment_id, load_id, gs.slot_index]))
                    if isempty(existing_sample)
                        sample_id = create_sample!(db;
                            experiment_id   = experiment_id,
                            name            = gs.name,
                            load_id         = load_id,
                            slot_index      = gs.slot_index,
                            grouping_source = gs.grouping_source,
                            name_source     = gs.name_source)
                    else
                        sample_id = Int(first(existing_sample).id)
                    end

                    for ge in gs.exposures
                        ge.stem in existing_filenames && continue  # insert-only
                        eid = create_exposure!(db;
                            experiment_id       = experiment_id,
                            sample_id           = sample_id,
                            filename            = ge.stem,
                            # image_path is required by the image route + prewarm_thumbnails!
                            # (both filter WHERE image_path IS NOT NULL).
                            image_path          = ge.tif_path,
                            prp_path            = ge.prp_path,
                            timestamp           = ismissing(ge.timestamp) ? nothing :
                                                  Dates.format(ge.timestamp, "yyyy-mm-ddTHH:MM:SS"),
                            exposure_time       = ismissing(ge.exposure_time_s) ? nothing : ge.exposure_time_s,
                            horizontal_position = ismissing(ge.horizontal_position_mm) ? nothing : ge.horizontal_position_mm,
                            scan_id             = ismissing(ge.scan_id)  ? nothing : ge.scan_id,
                            frame_no            = ismissing(ge.frame_no) ? nothing : ge.frame_no,
                            load_id             = load_id)
                        push!(new_exposure_ids, eid)
                        push!(existing_filenames, ge.stem)
                    end
                end
            end
        end
    end

    # -----------------------------------------------------------------------
    # 5. Analyze new exposures OUTSIDE the write transaction (same contract
    #    as cli_init_with_db!: a crash mid-analyze must not roll back ingest).
    # -----------------------------------------------------------------------
    if analyze
        n_new = length(new_exposure_ids)
        # Report against n_found (the discovery denominator), NOT n_new. Files
        # that already had rows are done the instant the persist txn commits, so
        # they count as processed. This keeps ONE scale for the whole scan and
        # fixes the clean-rescan case: with nothing new, n_new == 0, the loop
        # below never runs, and reporting `i / n_new` meant the bar showed
        # "0 / ~0" from start to finish and then simply vanished.
        base = max(0, n_found - n_new)
        # Coalesce to ~1% steps. One frame per exposure meant a 600-exposure
        # scan fired 600 frames into a 64-slot channel, and each frame costs the
        # client two invalidateQueries → a refetch of the growing loads payload.
        # The client falls behind, the channel saturates, and the subscriber is
        # evicted mid-scan (see _fanout_frame! in events.jl).
        step = max(1, n_new ÷ 100)
        for (i, eid) in enumerate(new_exposure_ids)
            try
                analyze_exposure!(db, eid)
            catch e
                @warn "scan_and_group!: analyze_exposure! failed" exposure_id=eid exception=e
            end
            # Progress fires OUTSIDE the structural txn (already committed above).
            if on_progress !== nothing && (i % step == 0 || i == n_new)
                on_progress(base + i, n_found)
            end
        end
        # Always land on a full bar, including the zero-new rescan that skips the
        # loop entirely.
        on_progress === nothing || on_progress(n_found, n_found)

        # Prewarm the thumbnail disk cache for the freshly-ingested exposures so the
        # first contact-sheet visit is fast (cold lazy generation staggers badly over
        # SMB — the exact case prewarm exists for, issue #261). The scan→ingest rewrite
        # dropped the old init/reingest call sites; this restores it on the live path.
        # overwrite=true defeats whole-second mtime granularity on a re-scan. Prewarm
        # is a non-essential cache warm — an unreadable TIFF must never fail the ingest
        # (the thumb then just generates lazily on first view).
        try
            # Scoped to the exposures this scan just inserted — the unscoped call
            # re-rendered every thumbnail in every experiment on every scan.
            prewarm_thumbnails!(db; overwrite = true, exposure_ids = new_exposure_ids)
        catch e
            @warn "scan_and_group!: thumbnail prewarm failed (non-fatal)" exception=e
        end
    end

    return (
        status          = :ok,
        added_loads     = length(grouping.loads),
        added_samples   = sum(length(l.samples) for l in grouping.loads; init = 0),
        added_exposures = length(new_exposure_ids),
    )
end

"""
    _DryRunRollback

Private sentinel thrown at the end of `regroup_experiment!`'s write transaction
when `dry_run=true`. SQLite.transaction rolls back on any throw, so all DB
writes in that transaction are undone. The sentinel is caught by the outer
`try/catch` in `regroup_experiment!` itself and swallowed (only real errors
propagate). The counters declared OUTSIDE the transaction closure survive the
rollback and are returned unchanged, so the caller sees what WOULD have
happened without any persistent effect.
"""
struct _DryRunRollback <: Exception end

"""
    regroup_experiment!(db, experiment_id; dry_run=false) -> NamedTuple

Retrofit-persist: make a pre-rework experiment's `loads`/`samples`/`exposures`
rows match the NEW auto-grouping partition (`scan_directory` → `group_into_samples`)
*while preserving curation*. The goal is a result indistinguishable from a fresh
`scan_and_group!` of the same directory. It ADOPTS existing rows where they map to
the partition: relinks existing exposures (keeping their ids → all exposure-keyed
curation survives; rewriting `filename` to the full scan stem so future rescans
dedup correctly), greedily reuses existing sample rows (keeping their human
names/notes), and for derived cells with NO existing row it INSERTS exposures +
creates samples exactly as `scan_and_group!` would (un-manifested files become
real, analyzed samples). Displaced-sample metadata is carried onto the absorbing
owner, then the now-empty displaced rows are deleted.

**No existing exposure is ever deleted** — adopted rows are `UPDATE`/relinked,
un-manifested files are `INSERT`ed (then analyzed outside the write transaction,
unless `analyze=false`). The only sample deletes are genuinely-displaced empty
rows, whose metadata is carried first (dedup-then-repoint, mirroring the merge
route §9.3). Idempotent: a second run is a no-op (loads dedup on
`(experiment_id, load_index)`, owners already claimed by lowest id, inserted
files now have rows so they relink, every UPDATE re-writes the same value).

When `dry_run=true` the transaction is rolled back via `_DryRunRollback` after
all writes — the DB is left untouched but the returned summary reflects what a
real run would have done.

Returns a NamedTuple summary:
`(status, loads_created, cells, samples_retrofitted, samples_created,
  samples_displaced, exposures_relinked, exposures_inserted, exposures_no_file,
  reshoots, geometry, discrepancies)` where `reshoots` = old samples split across
≥2 cells and `exposures_inserted` = un-manifested files newly ingested.
On an empty scan returns `(status=:empty, …)` and
writes nothing.
"""
function regroup_experiment!(db::SQLite.DB, experiment_id::Int; dry_run::Bool = false, analyze::Bool = true)
    # -----------------------------------------------------------------------
    # 1. Derive (read-only) — mirror scan_and_group! lines 68-98.
    # -----------------------------------------------------------------------
    exp_row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT data_dir, analysis_dir, image_pattern, metadata_pattern, integration_pattern FROM experiments WHERE id = ?", [experiment_id])))
    data_dir = String(exp_row.data_dir)

    # Resolve analysis_dir + file patterns + setup from DISK with the SAME resolver
    # the funnel runs at create (resolve_experiment_layout), so a migrated row is
    # indistinguishable from a freshly-ingested one. A pre-rework row arrives with
    # NULL patterns + a manifest-era analysis_dir; without this the migration falls
    # back to {name}.* defaults — which name the wrong integration files ({name}.dat
    # vs the real {name}_tot.dat) and produce full-stem filenames where a fresh scan
    # produces short ones. FILL-ONLY: a value already set on the row (a prior run /
    # user edit) wins over detection. Persisted at 3e′ so a rescan reuses them.
    # ponytail: this resolves one exposure per acquisition (single-frame, the
    # detected {name}_0_001.* pattern). Multi-frame is out of scope by design.
    _pick(rowval, detected, default) =
        (rowval === missing || rowval === nothing) ? something(detected, default) : String(rowval)
    dd     = rstrip(data_dir, '/')
    root   = basename(dd) == "data" ? dirname(dd) : dd
    layout = resolve_experiment_layout(root)
    analysis_dir = something(layout.analysis_dir, String(exp_row.analysis_dir))
    tif_pattern  = _pick(exp_row.image_pattern,       layout.image_pattern,       "{name}.tif")
    prp_pattern  = _pick(exp_row.metadata_pattern,    layout.metadata_pattern,    "{name}.prp")
    dat_pattern  = _pick(exp_row.integration_pattern, layout.integration_pattern, "{name}.dat")

    metas = scan_directory(data_dir, analysis_dir;
        tif_pattern = tif_pattern, prp_pattern = prp_pattern, dat_pattern = dat_pattern)

    isempty(metas) && return (
        status = :empty, loads_created = 0, cells = 0, samples_retrofitted = 0,
        samples_created = 0, samples_displaced = 0, exposures_relinked = 0,
        exposures_inserted = 0, exposures_no_file = 0, reshoots = 0,
        geometry = nothing, discrepancies = String[])

    prp_paths   = String[m.prp_path for m in metas if m.prp_path !== nothing]
    setup_files = layout.setup_file === nothing ? String[] : String[layout.setup_file]
    geo, disc = derive_geometry(prp_paths, setup_files)

    result = group_into_samples(metas)

    # -----------------------------------------------------------------------
    # 2. Build the cell list + byfile map (keyed on ge.stem == exposures.filename).
    # -----------------------------------------------------------------------
    cells = NamedTuple[]  # (load_index, slot_index, name, name_source, grouping_source, frame_count, stems)
    byfile = Dict{String, NamedTuple}()
    load_meta = Dict{Int, NamedTuple}()  # load_index => (frame_count, start_time, end_time)
    for gl in result.loads
        load_meta[gl.load_index] = (frame_count = gl.frame_count,
                                    start_time = gl.start_time, end_time = gl.end_time,
                                    session_index = gl.session_index)
        for gs in gl.samples
            push!(cells, (load_index = gl.load_index, slot_index = gs.slot_index,
                          name = gs.name, name_source = gs.name_source,
                          grouping_source = gs.grouping_source,
                          stems = String[ge.stem for ge in gs.exposures]))
            for ge in gs.exposures
                byfile[ge.stem] = (load_index = gl.load_index, slot_index = gs.slot_index,
                    tif_path = ge.tif_path, prp_path = ge.prp_path, timestamp = ge.timestamp,
                    exposure_time = ge.exposure_time_s,
                    horizontal_position = ge.horizontal_position_mm,
                    scan_id = ge.scan_id, frame_no = ge.frame_no)
            end
        end
    end
    # Deterministic, stable claim order: (load_index, slot_index).
    sort!(cells; by = c -> (c.load_index, c.slot_index))

    # Match each scan stem to the DB filename that represents it. The manifest-era
    # pipeline stored `exposures.filename` as the on-disk stem with the trailing
    # `_<rep>_<frame>` suffix STRIPPED (real dev-db: scan stem minus `_0_001`);
    # scan_directory keeps the suffix. A row may therefore be stored under the full
    # stem (a fresh scan or a prior post-fix re-run) OR its truncated form (manifest
    # era) — match either. On relink (3c) we rewrite `filename` to the full stem so
    # future rescans dedup on the same `(experiment_id, filename)` key.
    # ponytail: assumes one DB row per scan stem (curated dev-db is uniform
    # single-frame, verified 1:1). Multi-frame collapse (several full stems →
    # one truncated key) would need per-frame disambiguation; add when real data shows it.
    manifest_key(s) = replace(s, r"_\d+_\d+$" => "")
    existing_fn = Set(String(r.filename) for r in Tables.rowtable(DBInterface.execute(db,
        "SELECT DISTINCT filename FROM exposures WHERE experiment_id = ? AND filename IS NOT NULL",
        [experiment_id])))
    dbkey_of_stem = Dict{String, String}()  # scan stem => DB filename, only for stems with a row
    claimed_dk    = Dict{String, String}()  # DB filename => the scan stem that claimed it
    for s in keys(byfile)
        dk = (s in existing_fn) ? s :
             (manifest_key(s) in existing_fn ? manifest_key(s) : nothing)
        dk === nothing && continue
        # Fail loudly on a multi-frame collapse (two distinct full stems truncating
        # to one manifest key) rather than silently relinking only one frame and
        # dropping the other. Single-frame is the supported case (verified 1:1);
        # this is the fail-fast tripwire for when multi-frame data first appears.
        prior = get(claimed_dk, dk, nothing)
        prior === nothing || error("regroup_experiment!: scan stems \"$prior\" and " *
            "\"$s\" both map to DB filename \"$dk\" (multi-frame collapse); per-frame " *
            "relink is unsupported — single-frame only.")
        claimed_dk[dk] = s
        dbkey_of_stem[s] = dk
    end

    # Reshoot count (read-only): old samples whose stems span ≥2 derived cells — a
    # specimen the fresh-scan partition splits across loads. This is the
    # operator-eyeball metric (≈28 on the dev-db), NOT the total-load count.
    cell_of_dbkey = Dict{String, Int}()   # DB filename => cell index (via dbkey_of_stem)
    for (ci, cell) in enumerate(cells), s in cell.stems
        dk = get(dbkey_of_stem, s, nothing)
        dk === nothing || (cell_of_dbkey[dk] = ci)
    end
    old_cells = Dict{Int, Set{Int}}()   # old sample_id => set of cell indices its stems touch
    if !isempty(cell_of_dbkey)
        ph = join(fill("?", length(cell_of_dbkey)), ", ")
        for r in Tables.rowtable(DBInterface.execute(db,
                "SELECT filename, sample_id FROM exposures " *
                "WHERE experiment_id = ? AND sample_id IS NOT NULL AND filename IN ($ph)",
                vcat(Any[experiment_id], collect(keys(cell_of_dbkey)))))
            ci = get(cell_of_dbkey, String(r.filename), 0)
            ci == 0 && continue
            push!(get!(old_cells, Int(r.sample_id), Set{Int}()), ci)
        end
    end
    reshoots = count(cs -> length(cs) > 1, values(old_cells))

    # -----------------------------------------------------------------------
    # 3. Retrofit persist — all writes in ONE lock + transaction (atomic rollback).
    # -----------------------------------------------------------------------
    samples_retrofitted = 0
    samples_created     = 0
    samples_displaced   = 0
    exposures_relinked  = 0
    exposures_inserted  = 0
    exposures_no_file   = 0
    loads_created       = 0
    new_exposure_ids    = Int[]   # un-manifested files inserted this run (analyzed post-txn)

    try
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            # 3a. Loads — dedup on (experiment_id, load_index); reuse on re-run.
            # Persist frame_count + start/end times so a migrated load is identical
            # to a freshly-scanned one (scan_and_group! sets the same fields).
            fmt(t) = ismissing(t) ? nothing : Dates.format(t, "yyyy-mm-ddTHH:MM:SS")
            load_id_of = Dict{Int, Int}()
            for (li, m) in load_meta
                st = fmt(m.start_time); et = fmt(m.end_time)
                existing = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM loads WHERE experiment_id = ? AND load_index = ?",
                    [experiment_id, li]))
                if isempty(existing)
                    lid = create_load!(db; experiment_id = experiment_id,
                                       load_index = li, frame_count = m.frame_count,
                                       start_time = st, end_time = et,
                                       session_id = m.session_index)
                    loads_created += 1
                else
                    lid = Int(first(existing).id)
                    # Keep derived values fresh (idempotent re-write).
                    DBInterface.execute(db,
                        "UPDATE loads SET frame_count = ?, start_time = ?, end_time = ?, session_id = ? WHERE id = ?",
                        [m.frame_count, st, et, m.session_index, lid])
                end
                load_id_of[li] = lid
            end

            # 3b. Assign each cell a sample row — greedy reuse, lowest-id owner.
            claimed     = Set{Int}()
            cell_owner  = Vector{Int}(undef, length(cells))
            for (ci, cell) in enumerate(cells)
                lid = load_id_of[cell.load_index]
                # Distinct existing sample_ids of this cell's stems, minus claimed.
                # Match on the DB filename each stem maps to (full or truncated).
                candidates = Int[]
                cell_dbkeys = String[dbkey_of_stem[s] for s in cell.stems if haskey(dbkey_of_stem, s)]
                if !isempty(cell_dbkeys)
                    placeholders = join(fill("?", length(cell_dbkeys)), ", ")
                    rows = Tables.rowtable(DBInterface.execute(db,
                        "SELECT DISTINCT sample_id FROM exposures " *
                        "WHERE experiment_id = ? AND sample_id IS NOT NULL AND filename IN ($placeholders)",
                        vcat([experiment_id], cell_dbkeys)))
                    candidates = sort(Int[Int(r.sample_id) for r in rows if !(Int(r.sample_id) in claimed)])
                end
                if !isempty(candidates)
                    owner = first(candidates)  # lowest id
                    push!(claimed, owner)
                    # Stamp name_source='user' ONLY on a genuine first adoption — a
                    # pre-rework manifest sample carries a human label and has load_id
                    # NULL (no loads existed before regroup). A sample THIS migration
                    # auto-created on a prior run already has load_id set; re-claiming it
                    # must NOT flip its 'auto' label to 'user' (idempotency). SQLite reads
                    # the pre-update load_id in the CASE, so this is the right discriminator.
                    # Re-derive grouping_source too (symmetric with the create branch
                    # below) so a retrofitted sample is indistinguishable from a fresh one.
                    DBInterface.execute(db,
                        "UPDATE samples SET load_id = ?, slot_index = ?, grouping_source = ?, " *
                        "name_source = CASE WHEN load_id IS NULL THEN 'user' ELSE name_source END " *
                        "WHERE id = ?",
                        [lid, cell.slot_index, cell.grouping_source, owner])
                    samples_retrofitted += 1
                else
                    # Mirror scan_and_group! exactly (incl. name_source) so a created
                    # cell is indistinguishable from a freshly-scanned one.
                    owner = create_sample!(db; experiment_id = experiment_id,
                        name = cell.name, grouping_source = cell.grouping_source,
                        name_source = cell.name_source,
                        load_id = lid, slot_index = cell.slot_index)
                    push!(claimed, owner)
                    samples_created += 1
                end
                cell_owner[ci] = owner
            end

            # Stem → owner lookup for the relink.
            owner_of_stem = Dict{String, Int}()
            for (ci, cell) in enumerate(cells)
                for s in cell.stems
                    owner_of_stem[s] = cell_owner[ci]
                end
            end

            # 3c. Relink exposures. Capture old_sample_id → absorbing-owner tally
            # BEFORE sample_id is overwritten (it can't be reconstructed after).
            absorb_tally = Dict{Int, Dict{Int, Int}}()  # old_sid => (owner => count)
            for stem in keys(byfile)
                bf  = byfile[stem]
                own = owner_of_stem[stem]
                dk = get(dbkey_of_stem, stem, nothing)
                if dk === nothing
                    # Scanned file with no pre-rework row — INSERT it so it becomes a
                    # real sample (ingest-everything), exactly as scan_and_group! would.
                    eid = create_exposure!(db; experiment_id = experiment_id,
                        sample_id = own, filename = stem, image_path = bf.tif_path,
                        prp_path = ismissing(bf.prp_path) ? nothing : bf.prp_path,
                        timestamp = ismissing(bf.timestamp) ? nothing : Dates.format(bf.timestamp, "yyyy-mm-ddTHH:MM:SS"),
                        exposure_time = ismissing(bf.exposure_time) ? nothing : bf.exposure_time,
                        horizontal_position = ismissing(bf.horizontal_position) ? nothing : bf.horizontal_position,
                        scan_id = ismissing(bf.scan_id) ? nothing : bf.scan_id,
                        frame_no = ismissing(bf.frame_no) ? nothing : bf.frame_no,
                        load_id = load_id_of[bf.load_index])
                    push!(new_exposure_ids, eid)
                    exposures_inserted += 1
                    continue
                end
                # Old sample_id(s) currently on this stem's exposures (pre-overwrite).
                olds = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id, sample_id FROM exposures WHERE experiment_id = ? AND filename = ?",
                    [experiment_id, dk]))
                for r in olds
                    if r.sample_id !== missing && r.sample_id !== nothing
                        osid = Int(r.sample_id)
                        if osid != own
                            d = get!(absorb_tally, osid, Dict{Int, Int}())
                            d[own] = get(d, own, 0) + 1
                        end
                    end
                end
                # Rewrite filename to the full scan stem so future rescans dedup on
                # the same (experiment_id, filename) key (manifest era stored it truncated).
                DBInterface.execute(db, """
                    UPDATE exposures SET sample_id = ?, load_id = ?, filename = ?, prp_path = ?,
                        timestamp = ?, exposure_time = ?, horizontal_position = ?,
                        scan_id = ?, frame_no = ?
                    WHERE experiment_id = ? AND filename = ?
                """, [own, load_id_of[bf.load_index], stem,
                      ismissing(bf.prp_path) ? nothing : bf.prp_path,
                      ismissing(bf.timestamp) ? nothing : Dates.format(bf.timestamp, "yyyy-mm-ddTHH:MM:SS"),
                      ismissing(bf.exposure_time) ? nothing : bf.exposure_time,
                      ismissing(bf.horizontal_position) ? nothing : bf.horizontal_position,
                      ismissing(bf.scan_id) ? nothing : bf.scan_id,
                      ismissing(bf.frame_no) ? nothing : bf.frame_no,
                      experiment_id, dk])
                exposures_relinked += length(olds)
            end

            # Exposures with a row but no on-disk file: leave, count. Relinked rows
            # now carry the full scan stem (== a byfile key); only truly-gone files
            # keep a non-matching filename.
            allrows = Tables.rowtable(DBInterface.execute(db,
                "SELECT filename FROM exposures WHERE experiment_id = ? AND filename IS NOT NULL",
                [experiment_id]))
            exposures_no_file = count(r -> !haskey(byfile, String(r.filename)), allrows)

            # 3d. Displaced old samples — those whose every exposure was absorbed by a
            # lower-id sibling (so they now hold none). Restrict to ids we actually
            # absorbed FROM (a row emptied by some other cause is not ours to reap).
            now_empty = Set(Int(r.id) for r in Tables.rowtable(DBInterface.execute(db,
                """SELECT id FROM samples
                   WHERE experiment_id = ? AND merged_into_id IS NULL
                     AND id NOT IN (SELECT DISTINCT sample_id FROM exposures
                                    WHERE sample_id IS NOT NULL)""", [experiment_id])))
            for (old_sid, owners) in absorb_tally
                (old_sid in now_empty) || continue
                # Absorbing owner = the one that took the most of this old sample's exposures.
                owner = first(sort(collect(owners); by = p -> (-p.second, p.first))).first
                # Carry metadata (dedup-then-repoint, mirroring routes_grouping.jl:150-203).
                # series_samples: drop displaced membership where owner already in that series.
                DBInterface.execute(db,
                    """DELETE FROM series_samples
                       WHERE sample_id = ?
                         AND series_id IN (SELECT series_id FROM series_samples WHERE sample_id = ?)""",
                    [old_sid, owner])
                # sample_tags: drop displaced tag where owner already holds the key.
                DBInterface.execute(db,
                    """DELETE FROM sample_tags
                       WHERE sample_id = ?
                         AND key IN (SELECT key FROM sample_tags WHERE sample_id = ?)""",
                    [old_sid, owner])
                # Re-point the survivors of the dedup.
                DBInterface.execute(db,
                    "UPDATE series_samples SET sample_id = ? WHERE sample_id = ?", [owner, old_sid])
                DBInterface.execute(db,
                    "UPDATE sample_tags SET sample_id = ? WHERE sample_id = ?", [owner, old_sid])
                DBInterface.execute(db,
                    "UPDATE sample_messages SET sample_id = ? WHERE sample_id = ?", [owner, old_sid])
                # notes (a column on the samples row — NOT touched by the block above).
                dnote = first(Tables.rowtable(DBInterface.execute(db,
                    "SELECT notes FROM samples WHERE id = ?", [old_sid]))).notes
                if dnote !== missing && dnote !== nothing
                    DBInterface.execute(db,
                        "UPDATE samples SET notes = COALESCE(notes, ?) WHERE id = ?", [dnote, owner])
                end
                # Delete the now-empty displaced row.
                DBInterface.execute(db, "DELETE FROM samples WHERE id = ?", [old_sid])
                samples_displaced += 1
            end

            # 3d′. Reap dedup-orphaned ghosts. A same-image_path collision deleted by
            # migrate_exposures_experiment_id! (schema migration) can empty a sample
            # whose sole exposure was the non-survivor frame. That row is empty but
            # never in absorb_tally (its exposure was deleted, not absorbed) and never
            # relinked (no stem left to match), so it would linger as a ghost (load_id
            # NULL, 0 exposures) — invisible in load rollups but counted by
            # _experiment_stats and listed by /api/samples (both key on experiment_id),
            # over-reporting the migrated experiment. Reap any still-empty, never-grouped
            # row, dropping its sample-keyed children first (FK enforcement is ON).
            for r in Tables.rowtable(DBInterface.execute(db,
                    """SELECT id FROM samples
                       WHERE experiment_id = ? AND merged_into_id IS NULL AND load_id IS NULL
                         AND id NOT IN (SELECT DISTINCT sample_id FROM exposures
                                        WHERE sample_id IS NOT NULL)""", [experiment_id]))
                gid = Int(r.id)
                DBInterface.execute(db, "DELETE FROM series_samples WHERE sample_id = ?", [gid])
                DBInterface.execute(db, "DELETE FROM sample_tags    WHERE sample_id = ?", [gid])
                DBInterface.execute(db, "DELETE FROM sample_messages WHERE sample_id = ?", [gid])
                DBInterface.execute(db, "DELETE FROM samples WHERE id = ?", [gid])
                # Not counted as `displaced` — a reaped dedup-orphan was never part of
                # the new partition (no exposures to absorb), a distinct concept from a
                # sample whose frames were absorbed by a sibling.
            end

            # 3e. Geometry (never-clobber 'user').
            _fill_unset_geometry!(db, experiment_id, geo)

            # 3e′. Persist the file-resolved analysis_dir + patterns (FILL-ONLY via
            # COALESCE — a value already set is never clobbered) so the migrated row
            # is indistinguishable from a fresh ingest and later rescans reuse them.
            DBInterface.execute(db,
                """UPDATE experiments SET
                       analysis_dir        = COALESCE(analysis_dir, ?),
                       image_pattern       = COALESCE(image_pattern, ?),
                       metadata_pattern    = COALESCE(metadata_pattern, ?),
                       integration_pattern = COALESCE(integration_pattern, ?)
                   WHERE id = ?""",
                [layout.analysis_dir, layout.image_pattern, layout.metadata_pattern,
                 layout.integration_pattern, experiment_id])

            # 3f. Mark the experiment complete.
            DBInterface.execute(db,
                "UPDATE experiments SET ingest_status = 'complete', last_scanned_at = ? WHERE id = ?",
                [Dates.format(Dates.now(), "yyyy-mm-ddTHH:MM:SS"), experiment_id])

            # Dry-run rollback: throw the sentinel AFTER all writes so the
            # transaction rolls back everything, but the counters (declared
            # outside this closure) are already populated for the return value.
            dry_run && throw(_DryRunRollback())
        end
    end
    catch e
        e isa _DryRunRollback || rethrow()
    end

    # Analyze inserted exposures OUTSIDE the write transaction (same contract as
    # scan_and_group!: a crash mid-analyze must not roll back the retrofit). Skipped
    # on dry-run — those inserts were rolled back, so the ids no longer exist.
    if analyze && !dry_run
        for eid in new_exposure_ids
            try
                analyze_exposure!(db, eid)
            catch e
                @warn "regroup_experiment!: analyze_exposure! failed" exposure_id=eid exception=e
            end
        end
    end

    return (
        status              = :ok,
        loads_created       = loads_created,
        cells               = length(cells),
        samples_retrofitted = samples_retrofitted,
        samples_created     = samples_created,
        samples_displaced   = samples_displaced,
        exposures_relinked  = exposures_relinked,
        exposures_inserted  = exposures_inserted,
        exposures_no_file   = exposures_no_file,
        reshoots            = reshoots,
        geometry            = geo,
        discrepancies       = disc,
    )
end

"""
    cheap_change_check(db, experiment_id) -> Bool

Cheap "has the directory changed since the last ingest?" probe for the Phase-C
auto-rescan scheduler and the `POST /api/experiments/{id}/scan` route.

Returns `true` when the data directory appears to hold files not yet persisted
(so a `scan_and_group!` is warranted), `false` when it looks unchanged (the
scheduler tick can stay quiet / back off — spec §9.4).

**Cheapness contract:** this does NOT parse PRP files or run grouping. It counts
matching image files in the experiment's `data_dir` (a single `readdir` + suffix
filter) and compares against `COUNT(*)` of already-persisted exposures. Additive
ingest dedups on `(experiment_id, filename)`, so "more files on disk than rows in
the DB" is exactly "there is new data to ingest". The image suffix is read from
the experiment row's `image_pattern` column (same source `scan_and_group!` uses).

**Bias:** any ambiguity returns `true` (an extra scan is a harmless no-op via the
insert-only dedup; a missed scan would silently drop data). The only `false`
shortcuts are (a) a vanished/unreadable data dir — nothing to ingest, and a
scheduler tick must never crash on a missing volume — and (b) on-disk count ≤
persisted count.
"""
function cheap_change_check(
    db::SQLite.DB,
    experiment_id::Int,
)::Bool
    # Resolve the experiment's data_dir and image_pattern (the authoritative scan root + suffix).
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT data_dir, image_pattern FROM experiments WHERE id = ?", [experiment_id]))
    isempty(rows) && return true   # unknown experiment: let the caller's scan surface the error
    exp_row = first(rows)
    data_dir = String(exp_row.data_dir)

    isdir(data_dir) || return false  # vanished/unreadable dir: nothing to ingest, never crash

    # Cheap on-disk count: number of files whose name matches the image suffix.
    # image_pattern is "{name}<suffix>"; everything after "{name}" is the literal suffix.
    # Read from the experiment row — same source scan_and_group! uses.
    image_pattern = coalesce(exp_row.image_pattern, "{name}.tif")
    suffix = replace(image_pattern, "{name}" => "")
    on_disk = try
        count(f -> endswith(f, suffix), readdir(data_dir))
    catch
        return true   # readdir failed for any reason: assume changed (safe direction)
    end

    # Persisted count.
    persisted = Int(first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS n FROM exposures WHERE experiment_id = ?", [experiment_id]))).n)

    return on_disk > persisted
end

"""
    _fill_unset_geometry!(db, experiment_id, geo)

FILL-ONLY geometry write: write a derived field to the `experiments` row ONLY when
that field is still UNSET (its current `*_source` is `'default'`). Any field that
already has an established source — `'setup'`/`'prp'`/`'computed'` (derived once,
e.g. committed from the funnel preview) or `'user'` (a manual override) — is left
untouched. Also skips a field whose derived value is `missing`. Called inside a
write transaction.

Geometry is derived ONCE, at the funnel preview (with the confirmed setup file),
and committed at create; this scan step is the fallback that fills geometry for
callers that did not supply it, and it never silently re-derives over an
established value (so a rescan presents what was committed, not a fresh — and
possibly worse — re-derivation). The user can always re-edit in Configuration.
"""
function _fill_unset_geometry!(db::SQLite.DB, experiment_id::Int, geo::NamedTuple)
    current = first(Tables.rowtable(DBInterface.execute(db,
        """SELECT energy_kev_source, flight_path_m_source,
                  beam_center_x_source, beam_center_y_source, pixel_size_um_source
             FROM experiments WHERE id = ?""", [experiment_id])))

    updates = Pair{String,Any}[]

    function maybe!(field, source_field, val, source)
        curr_src = getproperty(current, Symbol(source_field))
        # Only fill a field that is still unset ('default'); never overwrite an
        # established source (setup/prp/computed/user).
        if !ismissing(curr_src) && String(curr_src) != "default"
            return
        end
        val === missing && return
        push!(updates, field => val)
        push!(updates, source_field => source)
    end

    maybe!("energy_kev",    "energy_kev_source",    geo.energy_kev,    geo.energy_kev_source)
    maybe!("flight_path_m", "flight_path_m_source", geo.flight_path_m, geo.flight_path_m_source)
    maybe!("beam_center_x", "beam_center_x_source", geo.beam_center_x, geo.beam_center_x_source)
    maybe!("beam_center_y", "beam_center_y_source", geo.beam_center_y, geo.beam_center_y_source)
    maybe!("pixel_size_um", "pixel_size_um_source", geo.pixel_size_um, geo.pixel_size_um_source)

    isempty(updates) && return

    cols = join(first.(updates) .* " = ?", ", ")
    vals = last.(updates)
    DBInterface.execute(db,
        "UPDATE experiments SET $cols WHERE id = ?",
        vcat(vals, [experiment_id]))
end
