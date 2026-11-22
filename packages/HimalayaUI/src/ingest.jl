# packages/HimalayaUI/src/ingest.jl
#
# scan_and_group!: the transactional ingest orchestrator (spec §9.1).
#
# Discipline cloned from _reingest_inner! (cli.jl) + analyze_exposure! (pipeline.jl):
#   - Insert-only for loads/samples/exposures (no clobber on rescan).
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

_nothing_if_missing(x) = ismissing(x) ? nothing : x

function backfill_monitors!(db::SQLite.DB, experiment_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT id, prp_path FROM exposures
          WHERE experiment_id = ? AND prp_path IS NOT NULL AND monitors_read = 0",
        [experiment_id]))
    isempty(rows) && return 0
    parsed = Tuple{Int,NamedTuple}[]
    for r in rows
        path = String(r.prp_path)
        isfile(path) || continue
        prp = try
            parse_prp(path)
        catch err
            @warn "backfill_monitors!: could not read" path exception = err
            continue
        end
        push!(parsed, (Int(r.id), _monitors(prp)))
    end
    isempty(parsed) && return 0
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            for (id, m) in parsed
                DBInterface.execute(db, """
                    UPDATE exposures
                       SET transmitted_counts = COALESCE(transmitted_counts, ?),
                           incident_monitor   = COALESCE(incident_monitor, ?),
                           counting_time      = COALESCE(counting_time, ?),
                           monitors_read      = 1
                     WHERE id = ?
                """, [_nothing_if_missing(m.transmitted_counts),
                      _nothing_if_missing(m.incident_monitor),
                      _nothing_if_missing(m.counting_time_s), id])
            end
        end
    end
    length(parsed)
end

"""
    scan_and_group!(db, experiment_id; analyze=true, on_progress=nothing)

Full ingest of every included setup of `experiment_id`, in setup order. Each
setup is scanned, grouped, persisted and analyzed on its own (`_scan_setup!`).
Returns a NamedTuple summary `(status, added_loads, added_samples, added_exposures)`;
`status` is `:empty` only when no setup found any image.
"""
function scan_and_group!(
    db           ::SQLite.DB,
    experiment_id::Int;
    analyze      ::Bool   = true,
    retry_failed ::Bool   = false,
    on_progress  ::Union{Function,Nothing} = nothing,
)
    results = [_scan_setup!(db, experiment_id, setup; analyze, retry_failed, on_progress)
               for setup in setups_of(db, experiment_id; included_only = true)]
    ok = filter(r -> r.status == :ok, results)
    isempty(ok) && return (status = :empty, added_loads = 0, added_samples = 0, added_exposures = 0)
    (status          = :ok,
     added_loads     = sum(r.added_loads for r in ok),
     added_samples   = sum(r.added_samples for r in ok),
     added_exposures = sum(r.added_exposures for r in ok))
end

function _scan_setup!(db::SQLite.DB, experiment_id::Int, setup;
                      analyze::Bool, retry_failed::Bool, on_progress)
    setup_id     = Int(setup.id)
    data_dir     = String(setup.data_dir)
    analysis_dir = String(setup.analysis_dir)
    tif_pattern = coalesce(setup.image_pattern,       "{name}.tif")
    prp_pattern = coalesce(setup.metadata_pattern,    "{name}.prp")
    dat_pattern = coalesce(setup.integration_pattern, "{name}.dat")

    # -----------------------------------------------------------------------
    # 1. Scan: enumerate TIF+PRP+DAT triplets
    # -----------------------------------------------------------------------
    # `on_progress(processed, total, stage)` — stage names the segment the UI
    # should be filling ("discovery" / "analyzing"). Discovery is
    # the expensive silent half (two stats + a PRP parse per exposure over SMB),
    # so it reports from inside scan_directory rather than only announcing its
    # total afterwards.
    metas = scan_directory(data_dir, analysis_dir;
        tif_pattern = tif_pattern,
        prp_pattern = prp_pattern,
        dat_pattern = dat_pattern,
        on_progress = on_progress)

    any(m -> m.prp !== nothing, metas) && filter!(m -> m.prp !== nothing, metas)
    isempty(metas) && return (status = :empty, added_loads = 0, added_samples = 0, added_exposures = 0)

    setup_files = _find_setup_files(analysis_dir)

    # Reuse the PRPs scan_directory already parsed. The path-based
    # derive_geometry re-reads and re-parses every one of them — a second full
    # pass over every PRP in the experiment, per ingest.
    # NamedTuple[...] required: ExposureMeta.prp is Union{NamedTuple,Nothing}, so
    # when NO exposure has a PRP (tif-without-prp layouts, or a misconfigured
    # prp_pattern) the filtered comprehension collects to Vector{Any} and
    # derive_geometry MethodErrors -- turning "geometry falls back to defaults"
    # into "the ingest crashes". Same trap as geometry.jl's delegation.
    parsed_prps = NamedTuple[m.prp for m in metas if m.prp !== nothing]

    geo, _disc = derive_geometry(parsed_prps, setup_files)

    # FILL-ONLY: write only geometry fields still unset ('default' source). Fields
    # already established at create (committed from the funnel preview, which
    # derived them with the confirmed setup file) or manually edited ('user') are
    # left untouched — geometry is derived once, not re-derived on every scan.
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            _fill_unset_geometry!(db, setup_id, geo)
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
            existing_filenames = Set{String}(String(r.filename)
                for r in Tables.rowtable(DBInterface.execute(db,
                    "SELECT filename FROM exposures WHERE setup_id = ? AND filename IS NOT NULL",
                    [setup_id])))

            for gl in grouping.loads
                any(ge -> !(ge.stem in existing_filenames), (ge for gs in gl.samples for ge in gs.exposures)) || continue
                existing_load = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM loads WHERE setup_id = ? AND load_index = ?",
                    [setup_id, gl.load_index]))
                if isempty(existing_load)
                    load_id = Int(DBInterface.lastrowid(DBInterface.execute(db, """
                        INSERT INTO loads (experiment_id, setup_id, load_index, frame_count,
                                           start_time, end_time, session_id)
                        VALUES (?, ?, ?, ?, ?, ?, ?)
                    """, [experiment_id, setup_id, gl.load_index, gl.frame_count,
                          ismissing(gl.start_time) ? nothing : Dates.format(gl.start_time, "yyyy-mm-ddTHH:MM:SS"),
                          ismissing(gl.end_time)   ? nothing : Dates.format(gl.end_time,   "yyyy-mm-ddTHH:MM:SS"),
                          gl.session_index])))
                else
                    load_id = Int(first(existing_load).id)
                end

                for gs in gl.samples
                    any(ge -> !(ge.stem in existing_filenames), gs.exposures) || continue
                    existing_sample = Tables.rowtable(DBInterface.execute(db,
                        "SELECT id FROM samples WHERE setup_id = ? AND load_id = ? AND slot_index = ?",
                        [setup_id, load_id, gs.slot_index]))
                    if isempty(existing_sample)
                        sample_id = create_sample!(db;
                            experiment_id   = experiment_id,
                            setup_id        = setup_id,
                            name            = gs.name,
                            load_id         = load_id,
                            slot_index      = gs.slot_index,
                            grouping_source = gs.grouping_source,
                            name_source     = gs.name_source)
                    else
                        sample_id = Int(first(existing_sample).id)
                    end

                    for ge in gs.exposures
                        ge.stem in existing_filenames && continue
                        eid = create_exposure!(db;
                            experiment_id       = experiment_id,
                            setup_id            = setup_id,
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
                            load_id             = load_id,
                            transmitted_counts  = _nothing_if_missing(ge.transmitted_counts),
                            incident_monitor    = _nothing_if_missing(ge.incident_monitor),
                            counting_time       = _nothing_if_missing(ge.counting_time_s))
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
        to_analyze = vcat(new_exposure_ids,
            setdiff(_pending_setup_ids(db, setup; retry_failed), new_exposure_ids))
        n_new = length(to_analyze)
        # The "analyzing" segment reports against n_new, its OWN denominator — the
        # segmented bar gives each stage its own track, so stages no longer have to
        # share one scale (which is what forced the old cross-stage bookkeeping here).
        #
        # Cap the frame count at ~100 per stage. One frame per exposure meant a
        # 600-exposure scan fired 600 frames into a 64-slot channel, and each
        # frame costs the client a cache invalidation. NOTE this only thins large
        # scans: below n_new = 100 the step is 1, i.e. still one frame per
        # exposure, and 65..100 exposures still exceeds the 64-slot channel.
        # Saturation is reduced, not eliminated — correctness under saturation
        # rests on _fanout_frame! (events.jl) closing an evicted subscriber.
        step = max(1, n_new ÷ 100)
        for (i, eid) in enumerate(to_analyze)
            analyze_or_record!(db, eid, analysis_dir)
            # Progress fires OUTSIDE the structural txn (already committed above).
            if on_progress !== nothing && (i % step == 0 || i == n_new)
                on_progress(i, n_new, "analyzing")
            end
        end
        # Always land the segment full, including the clean rescan where n_new == 0
        # and the loop above never runs (a 0/0 segment reads as complete, not as a
        # stalled bar).
        on_progress === nothing || on_progress(n_new, n_new, "analyzing")
    end

    return (
        status          = :ok,
        added_loads     = length(grouping.loads),
        added_samples   = sum(length(l.samples) for l in grouping.loads; init = 0),
        added_exposures = length(new_exposure_ids),
    )
end

const NO_TRACE = "no trace"

function _file_finder(base::AbstractString, pattern::AbstractString; strip_acquisition::Bool = false)
    pattern_dir = dirname(pattern)
    if occursin("{name}", pattern_dir)
        return function (stem)
            p = joinpath(base, replace(pattern, "{name}" => stem))
            isfile(p) ? p : nothing
        end
    end
    dir = joinpath(base, pattern_dir)
    names = try
        Set(readdir(dir))
    catch
        Set{String}()
    end
    leaf = basename(pattern)
    function (stem)
        candidate = replace(leaf, "{name}" => stem)
        candidate in names && return joinpath(dir, candidate)
        if strip_acquisition
            stripped = acquisition_stem(stem)
            if stripped != stem
                alt = replace(leaf, "{name}" => stripped)
                alt in names && return joinpath(dir, alt)
            end
        end
        nothing
    end
end

_trace_finder(setup) =
    _file_finder(String(setup.analysis_dir), integration_pattern_for(setup); strip_acquisition = true)

function _stat_signature(path::AbstractString)
    st = stat(path)
    string(st.size, ":", st.mtime)
end

function analyze_or_record!(db::SQLite.DB, exposure_id::Int, analysis_dir::AbstractString)
    row = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT filename FROM exposures WHERE id = ?", [exposure_id])))
    trace = resolve_trace_path(analysis_dir, integration_pattern_for_exposure(db, exposure_id), String(row.filename))
    error_text = if trace === nothing
        NO_TRACE
    else
        try
            analyze_exposure!(db, exposure_id, String(analysis_dir))
            nothing
        catch e
            @warn "analyze_exposure! failed" exposure_id exception = e
            sprint(showerror, e)
        end
    end
    error_text === nothing && return true
    signature = error_text == NO_TRACE ? nothing : (try _stat_signature(trace) catch; nothing end)
    lock(_DB_WRITE_LOCK) do
        DBInterface.execute(db,
            "UPDATE exposures SET analysis_error = ?, analysis_version = ?, analysis_trace_stat = ? WHERE id = ?",
            [first(error_text, 500), error_text == NO_TRACE ? nothing : ANALYSIS_VERSION, signature, exposure_id])
    end
    false
end

function ingest_problems(db::SQLite.DB, experiment_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT e.id, e.filename, e.analysis_error, e.sample_id, s.name AS sample_name
          FROM exposures e LEFT JOIN samples s ON s.id = e.sample_id
         WHERE e.experiment_id = ? AND e.analysis_error IS NOT NULL
         ORDER BY e.sample_id, e.filename, e.id""", [experiment_id]))
    [Dict(:exposure_id => Int(r.id),
          :filename => coalesce(r.filename, nothing),
          :sample_id => ismissing(r.sample_id) ? nothing : Int(r.sample_id),
          :sample_name => coalesce(r.sample_name, nothing),
          :problem => r.analysis_error == NO_TRACE ? "no_trace" : "failed",
          :reason => String(r.analysis_error)) for r in rows]
end

pending_analysis_ids(db::SQLite.DB, experiment_id::Int; retry_failed::Bool = false) =
    reduce(vcat, (_pending_setup_ids(db, st; retry_failed)
                  for st in setups_of(db, experiment_id; included_only = true)); init = Int[])

function _pending_setup_ids(db::SQLite.DB, setup; retry_failed::Bool = false)
    rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT id, filename, analysis_error, analysis_version, analysis_trace_stat FROM exposures
         WHERE setup_id = ?2 AND filename IS NOT NULL
           AND (trace_hash IS NULL OR analysis_version IS NULL OR analysis_version != ?1)
         ORDER BY id""", [ANALYSIS_VERSION, Int(setup.id)]))
    isempty(rows) && return Int[]
    find_trace = _trace_finder(setup)
    function due(r)
        ismissing(r.analysis_error) && return true
        trace = find_trace(String(r.filename))
        r.analysis_error == NO_TRACE && return trace !== nothing
        coalesce(r.analysis_version, "") != ANALYSIS_VERSION && return true
        trace === nothing && return false
        retry_failed && return true
        signature = try _stat_signature(trace) catch; nothing end
        signature !== nothing && signature != coalesce(r.analysis_trace_stat, "")
    end
    Int[Int(r.id) for r in rows if due(r)]
end

_insertable_stems(db::SQLite.DB, experiment_id::Int) =
    reduce(vcat, (_insertable_setup_stems(db, st)
                  for st in setups_of(db, experiment_id; included_only = true)); init = String[])

function _insertable_setup_stems(db::SQLite.DB, r)
    data_dir = String(r.data_dir)
    listing = try
        tif_listing(data_dir, coalesce(r.image_pattern, "{name}.tif"))
    catch
        nothing
    end
    listing === nothing && return String[]
    _, tif_files, tif_suffix = listing
    stems = String[f[1:end-length(tif_suffix)] for f in tif_files]
    known = Set(String(x.filename) for x in Tables.rowtable(DBInterface.execute(db,
        "SELECT filename FROM exposures WHERE setup_id = ? AND filename IS NOT NULL", [Int(r.id)])))
    find_prp = _file_finder(data_dir, coalesce(r.metadata_pattern, "{name}.prp"))
    has_prp(stem) = find_prp(stem) !== nothing
    any_prp = any(has_prp, stems)
    String[stem for stem in stems if !(stem in known) && (!any_prp || has_prp(stem))]
end

"""
    ingest_has_work(db, experiment_id) -> Bool

True when a pass would change something: a tif (with its prp, when the folder
has prps) that has no exposure row yet, or an exposure `pending_analysis_ids`
returns. A few directory reads plus a stat per failed exposure, so it is cheap
enough to run inside a request.
"""
function ingest_has_work(db::SQLite.DB, experiment_id::Int)::Bool
    !isempty(_insertable_stems(db, experiment_id)) || !isempty(pending_analysis_ids(db, experiment_id)) ||
        !isempty(new_setup_roots(db, experiment_id))
end

"""
    _fill_unset_geometry!(db, setup_id, geo)

FILL-ONLY geometry write: write a derived field to the `setups` row ONLY when
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
function _fill_unset_geometry!(db::SQLite.DB, setup_id::Int, geo::NamedTuple)
    current = first(Tables.rowtable(DBInterface.execute(db,
        """SELECT energy_kev_source, flight_path_m_source,
                  beam_center_x_source, beam_center_y_source, pixel_size_um_source
             FROM setups WHERE id = ?""", [setup_id])))

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
        "UPDATE setups SET $cols WHERE id = ?",
        vcat(vals, [setup_id]))
end
