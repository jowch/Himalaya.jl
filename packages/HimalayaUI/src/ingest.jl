# packages/HimalayaUI/src/ingest.jl
#
# scan_and_group!: the transactional ingest orchestrator (spec §9.1).
#
# Discipline cloned from _reingest_inner! (cli.jl) + analyze_exposure! (pipeline.jl):
#   - Insert-only for loads/samples/exposures (no clobber on rescan).
#   - The exposure dedup key is (experiment_id, filename) — also enforced by the
#     `exposures_unique_filename` UNIQUE index (db.jl migrate_exposures_unique_filename!).
#   - Geometry never clobbers a `*_source='human'` field (spec §4).
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
    scan_and_group!(db, experiment_id, root_dir; additive=true, analyze=true,
                    tif_pattern="{name}.tif", prp_pattern="{name}.prp", dat_pattern="{name}.dat")

Full ingest of a beamtime directory into `db` under `experiment_id`.

Steps:
  1. Resolve `data_dir` and `analysis_dir` from the `experiments` row.
  2. Scan: enumerate TIF+PRP+DAT triplets with `scan_directory`.
  3. Geometry: derive from the sampled PRPs + the latest setup_info file; write
     to the `experiments` row (respects `*_source='human'` never-clobber).
  4. Group: `group_into_samples` → loads/samples/slots.
  5. Persist: insert `loads`, `samples`, `exposures` rows inside ONE write
     transaction (insert-only — dedup loads by (experiment_id, load_index),
     samples by (experiment_id, load_id, slot_index), exposures by
     (experiment_id, filename) — so a clean rescan is a true no-op).
  6. Analyze: run `analyze_exposure!` for every newly inserted exposure, OUTSIDE
     the write transaction (same contract as `cli_init_with_db!`). Skipped when
     `analyze=false` (tests + the HTTP "scan-only" path).

`additive` is accepted for API symmetry with the rescan caller (Task 10/11); the
insert-only discipline means the orchestrator is always additive, so the flag is
currently advisory. Returns a NamedTuple summary
`(status, added_loads, added_samples, added_exposures)`.
"""
function scan_and_group!(
    db           ::SQLite.DB,
    experiment_id::Int,
    root_dir     ::AbstractString;
    additive     ::Bool   = true,
    analyze      ::Bool   = true,
    tif_pattern  ::String = "{name}.tif",
    prp_pattern  ::String = "{name}.prp",
    dat_pattern  ::String = "{name}.dat",
)
    root_dir = abspath(root_dir)

    # Resolve data_dir and analysis_dir from the experiment row.
    exp_row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT data_dir, analysis_dir FROM experiments WHERE id = ?", [experiment_id])))
    data_dir     = String(exp_row.data_dir)
    analysis_dir = String(exp_row.analysis_dir)

    # -----------------------------------------------------------------------
    # 1. Scan: enumerate TIF+PRP+DAT triplets
    # -----------------------------------------------------------------------
    metas = scan_directory(data_dir, analysis_dir;
        tif_pattern = tif_pattern,
        prp_pattern = prp_pattern,
        dat_pattern = dat_pattern)

    isempty(metas) && return (status = :empty, added_loads = 0, added_samples = 0, added_exposures = 0)

    # -----------------------------------------------------------------------
    # 2. Geometry: derive + write to experiments row (never-clobber human fields)
    # -----------------------------------------------------------------------
    prp_paths   = String[m.prp_path for m in metas if m.prp_path !== nothing]
    setup_files = isdir(analysis_dir) ?
        String[joinpath(analysis_dir, f)
                for f in readdir(analysis_dir)
                if startswith(f, "setup_info_") && endswith(f, ".txt")] :
        String[]

    geo, _disc = derive_geometry(prp_paths, setup_files)

    # Persist geometry with never-clobber: only write fields whose current source
    # is not 'human' (spec §4).
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            _update_geometry_if_not_human!(db, experiment_id, geo)
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
                                           start_time, end_time)
                        VALUES (?, ?, ?, ?, ?)
                    """, [experiment_id, gl.load_index, gl.frame_count,
                          ismissing(gl.start_time) ? nothing : Dates.format(gl.start_time, "yyyy-mm-ddTHH:MM:SS"),
                          ismissing(gl.end_time)   ? nothing : Dates.format(gl.end_time,   "yyyy-mm-ddTHH:MM:SS")])))
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
        for eid in new_exposure_ids
            try
                analyze_exposure!(db, eid)
            catch e
                @warn "scan_and_group!: analyze_exposure! failed" exposure_id=eid exception=e
            end
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
    cheap_change_check(db, experiment_id, root_dir) -> Bool

Cheap "has the directory changed since the last ingest?" probe for the Phase-C
auto-rescan scheduler and the `POST /api/experiments/{id}/scan` route (both
resolve this by name via `isdefined(HimalayaUI, :cheap_change_check)`).

Returns `true` when the data directory appears to hold files not yet persisted
(so a `scan_and_group!` is warranted), `false` when it looks unchanged (the
scheduler tick can stay quiet / back off — spec §9.4).

**Cheapness contract:** this does NOT parse PRP files or run grouping. It counts
matching image files in the experiment's `data_dir` (a single `readdir` + suffix
filter) and compares against `COUNT(*)` of already-persisted exposures. Additive
ingest dedups on `(experiment_id, filename)`, so "more files on disk than rows in
the DB" is exactly "there is new data to ingest".

**Bias:** any ambiguity returns `true` (an extra scan is a harmless no-op via the
insert-only dedup; a missed scan would silently drop data). The only `false`
shortcuts are (a) a vanished/unreadable data dir — nothing to ingest, and a
scheduler tick must never crash on a missing volume — and (b) on-disk count ≤
persisted count.

`root_dir` is accepted to match the Phase-C call contract
(`cheap_change_check(db, experiment_id, root_dir)`); the authoritative scan root
is the experiment's stored `data_dir`, so `root_dir` is currently advisory only.
"""
function cheap_change_check(
    db::SQLite.DB,
    experiment_id::Int,
    root_dir::AbstractString;
    image_pattern::String = "{name}.tif",
)::Bool
    # Resolve the experiment's data_dir (the authoritative scan root).
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT data_dir FROM experiments WHERE id = ?", [experiment_id]))
    isempty(rows) && return true   # unknown experiment: let the caller's scan surface the error
    data_dir = String(first(rows).data_dir)

    isdir(data_dir) || return false  # vanished/unreadable dir: nothing to ingest, never crash

    # Cheap on-disk count: number of files whose name matches the image suffix.
    # image_pattern is "{name}<suffix>"; everything after "{name}" is the literal suffix.
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
    _update_geometry_if_not_human!(db, experiment_id, geo)

Write derived geometry fields to the `experiments` row, skipping any field whose
current `*_source` is `'human'` (never-clobber contract, spec §4) and any field
whose derived value is `missing`. Called inside a write transaction.
"""
function _update_geometry_if_not_human!(db::SQLite.DB, experiment_id::Int, geo::NamedTuple)
    current = first(Tables.rowtable(DBInterface.execute(db,
        """SELECT energy_kev_source, flight_path_m_source,
                  beam_center_x_source, beam_center_y_source, pixel_size_um_source
             FROM experiments WHERE id = ?""", [experiment_id])))

    updates = Pair{String,Any}[]

    function maybe!(field, source_field, val, source)
        curr_src = getproperty(current, Symbol(source_field))
        if !ismissing(curr_src) && String(curr_src) == "human"
            return  # never overwrite a human-set field
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
