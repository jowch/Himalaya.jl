function _find_setup_files(analysis_dir::AbstractString)
    dir = analysis_dir
    for _ in 1:4
        if isdir(dir)
            hits = String[joinpath(dir, f) for f in list_names(dir)
                          if startswith(f, "setup_info_") && endswith(f, ".txt")]
            isempty(hits) || return hits
        end
        parent = dirname(dir)
        parent == dir && break
        dir = parent
    end
    return String[]
end

_nothing_if_missing(x) = ismissing(x) ? nothing : x

function backfill_monitors!(db::SQLite.DB, experiment_id::Int)
    read_trace_normalisations!(db; experiment_id)
    rows = sqlrows(DBInterface.execute(db,
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
                           pipe_length_m      = COALESCE(pipe_length_m, ?),
                           monitors_read      = 1
                     WHERE id = ?
                """, [_nothing_if_missing(m.transmitted_counts),
                      _nothing_if_missing(m.incident_monitor),
                      _nothing_if_missing(m.counting_time_s),
                      _nothing_if_missing(m.pipe_length_m), id])
            end
        end
    end
    length(parsed)
end

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
    setup = settle_setup_mode!(db, setup)
    analyze = analyze && !is_autosampler(setup)
    setup_id     = Int(setup.id)
    data_dir     = String(setup.data_dir)
    analysis_dir = String(setup.analysis_dir)
    tif_pattern = coalesce(setup.image_pattern,       "{name}.tif")
    prp_pattern = coalesce(setup.metadata_pattern,    "{name}.prp")
    dat_pattern = coalesce(setup.integration_pattern, "{name}.dat")

    metas = scan_directory(data_dir, analysis_dir;
        tif_pattern = tif_pattern,
        prp_pattern = prp_pattern,
        dat_pattern = dat_pattern,
        on_progress = on_progress)

    any(m -> m.prp !== nothing, metas) && filter!(m -> m.prp !== nothing, metas)
    isempty(metas) && return (status = :empty, added_loads = 0, added_samples = 0, added_exposures = 0)

    setup_files = _find_setup_files(analysis_dir)

    parsed_prps = NamedTuple[m.prp for m in metas if m.prp !== nothing]

    geo, _disc = derive_geometry(parsed_prps, setup_files)

    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            _fill_unset_geometry!(db, setup_id, geo)
        end
    end

    grouping = is_autosampler(setup) ? group_autosampler(metas, data_dir) : group_into_samples(metas)

    new_exposure_ids = Int[]
    reuse_samples = !is_autosampler(setup)

    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            existing_filenames = Set{String}(String(r.filename)
                for r in sqlrows(DBInterface.execute(db,
                    "SELECT filename FROM exposures WHERE setup_id = ? AND filename IS NOT NULL",
                    [setup_id])))

            for gl in grouping.loads
                any(ge -> !(ge.stem in existing_filenames), (ge for gs in gl.samples for ge in gs.exposures)) || continue
                existing_load = sqlrows(DBInterface.execute(db,
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
                    existing_sample = reuse_samples ? sqlrows(DBInterface.execute(db,
                        "SELECT id FROM samples WHERE setup_id = ? AND load_id = ? AND slot_index = ?",
                        [setup_id, load_id, gs.slot_index])) : NamedTuple[]
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
                            counting_time       = _nothing_if_missing(ge.counting_time_s),
                            pipe_length_m       = _nothing_if_missing(ge.pipe_length_m))
                        push!(new_exposure_ids, eid)
                        push!(existing_filenames, ge.stem)
                    end
                end
            end
        end
    end

    _keep_custom_selections!(db, new_exposure_ids)
    if analyze
        to_analyze = vcat(new_exposure_ids,
            setdiff(_pending_setup_ids(db, setup; retry_failed), new_exposure_ids))
        n_new = length(to_analyze)
        step = max(1, n_new ÷ 100)
        for (i, eid) in enumerate(to_analyze)
            analyze_or_record!(db, eid, analysis_dir)
            if on_progress !== nothing && (i % step == 0 || i == n_new)
                on_progress(i, n_new, "analyzing")
            end
        end
        on_progress === nothing || on_progress(n_new, n_new, "analyzing")
        deselect_beam_off!(db, to_analyze)
        for r in sqlrows(DBInterface.execute(db,
                "SELECT id FROM samples WHERE setup_id = ? ORDER BY id", [setup_id]))
            refresh_selection!(db, Int(r.id))
        end
    end

    return (
        status          = :ok,
        added_loads     = length(grouping.loads),
        added_samples   = sum(length(l.samples) for l in grouping.loads; init = 0),
        added_exposures = length(new_exposure_ids),
    )
end

const NO_TRACE = "no trace"

_file_finder(base::AbstractString, pattern::AbstractString; strip_acquisition::Bool = false) =
    _file_finder(String(base), String(pattern); strip_acquisition)

function _file_finder(base::String, pattern::String; strip_acquisition::Bool = false)
    pattern_dir = dirname(pattern)
    if occursin("{name}", pattern_dir)
        return function (stem)
            p = joinpath(base, replace(pattern, "{name}" => String(stem)))
            isfile(p) ? p : nothing
        end
    end
    dir = joinpath(base, pattern_dir)
    names = Set(list_names(dir))
    leaf = basename(pattern)
    function (stem)
        stem = String(stem)
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
    row = first(sqlrows(DBInterface.execute(db,
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
    rows = sqlrows(DBInterface.execute(db, """
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
    is_autosampler(setup) && return Int[]
    rows = sqlrows(DBInterface.execute(db, """
        SELECT id, filename, analysis_error, analysis_version, analysis_trace_stat FROM exposures
         WHERE setup_id = ?2 AND filename IS NOT NULL
           AND (reviewed_at IS NULL
                AND (trace_hash IS NULL OR analysis_version IS NULL OR analysis_version != ?1)
             OR reviewed_at IS NOT NULL AND trace_hash IS NULL
                AND NOT EXISTS (SELECT 1 FROM assignment_members m WHERE m.exposure_id = exposures.id)
                AND NOT EXISTS (SELECT 1 FROM auto_peaks a WHERE a.exposure_id = exposures.id))
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
    layout = is_autosampler(r) && !_setup_has_exposures(db, Int(r.id)) ? autosampler_layout(String(r.root_dir)) : nothing
    image_pattern = layout === nothing ? coalesce(r.image_pattern, "{name}.tif") : layout.image_pattern
    metadata_pattern = layout === nothing ? coalesce(r.metadata_pattern, "{name}.prp") : layout.metadata_pattern
    listing = try
        tif_listing(data_dir, image_pattern)
    catch
        nothing
    end
    listing === nothing && return String[]
    _, tif_files, tif_suffix = listing
    stems = String[f[1:end-length(tif_suffix)] for f in tif_files]
    known = Set(String(x.filename) for x in sqlrows(DBInterface.execute(db,
        "SELECT filename FROM exposures WHERE setup_id = ? AND filename IS NOT NULL", [Int(r.id)])))
    find_prp = _file_finder(data_dir, metadata_pattern)
    has_prp(stem) = find_prp(stem) !== nothing
    any_prp = any(has_prp, stems)
    String[stem for stem in stems if !(stem in known) && (!any_prp || has_prp(stem))]
end

function ingest_has_work(db::SQLite.DB, experiment_id::Int)::Bool
    !isempty(_insertable_stems(db, experiment_id)) || !isempty(pending_analysis_ids(db, experiment_id)) ||
        !isempty(new_setup_roots(db, experiment_id))
end

function _fill_unset_geometry!(db::SQLite.DB, setup_id::Int, geo::NamedTuple)
    current = first(sqlrows(DBInterface.execute(db,
        """SELECT energy_kev_source, flight_path_m_source,
                  beam_center_x_source, beam_center_y_source, pixel_size_um_source
             FROM setups WHERE id = ?""", [setup_id])))

    updates = Pair{String,Any}[]

    function maybe!(field, source_field, val, source)
        curr_src = getproperty(current, Symbol(source_field))
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

function _keep_custom_selections!(db::SQLite.DB, new_ids)
    isempty(new_ids) && return nothing
    ph = join(fill("?", length(new_ids)), ",")
    custom = Set{Int}()
    for r in sqlrows(DBInterface.execute(db, """
            SELECT o.sample_id AS sid, o.trace_median, o.analysis_error, o.incident_monitor FROM exposures o
             WHERE o.sample_id IN (SELECT n.sample_id FROM exposures n WHERE n.id IN ($ph) AND n.sample_id IS NOT NULL)
               AND o.kind = 'file' AND o.deselected = 1
               AND (o.status IS NULL OR o.status != 'rejected')
               AND o.id NOT IN ($ph)""", vcat(new_ids, new_ids)))
        is_beam_off(r.trace_median, r.analysis_error, r.incident_monitor) || push!(custom, Int(r.sid))
    end
    for sid in sort!(collect(custom))
        DBInterface.execute(db,
            "UPDATE exposures SET deselected = 1 WHERE sample_id = ? AND id IN ($ph)", vcat(sid, new_ids))
    end
    nothing
end
