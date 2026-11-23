struct SkipExposure <: Exception
    reason::Symbol
end

const SKIP_TEXT = Dict(
    :no_group => "no experiment or sample id",
    :no_setup => "no setup",
    :no_trace_file => "trace file not found",
    :unreadable_trace => "trace file not parseable")

function _load_label_trace(db::SQLite.DB, r)
    setup = exposure_setup(db, r.id)
    setup === nothing && throw(SkipExposure(:no_setup))
    path = resolve_trace_path(String(setup.analysis_dir), integration_pattern_for(setup), String(r.filename))
    path === nothing && throw(SkipExposure(:no_trace_file))
    try
        load_dat(path)
    catch err
        err isa Union{ErrorException, ArgumentError, DimensionMismatch, BoundsError} || rethrow()
        throw(SkipExposure(:unreadable_trace))
    end
end

function _latest_snapshot(db::SQLite.DB, exposure_id)
    rows = sqlrows(DBInterface.execute(db,
        "SELECT peaks_json FROM peak_reviews WHERE exposure_id = ? AND direct = 1 ORDER BY id DESC LIMIT 1", [exposure_id]))
    isempty(rows) ? nothing : JSON3.read(String(rows[1].peaks_json))
end

function labeled_exposures(db::SQLite.DB; group_by::Symbol = :experiment,
                           include_form_factor::Bool = true, include_unlabelled::Bool = false, io::IO = stderr)
    group_by in (:experiment, :sample) || throw(ArgumentError("group_by must be :experiment or :sample"))
    _table_exists(db, "peak_reviews") ||
        error("the database has no review snapshots; open it once with the current HimalayaUI to migrate it")
    out = HimalayaPeakModel.LabeledExposure[]
    skips = Dict{Symbol, Int}()
    rows = sqlrows(DBInterface.execute(db, """
        SELECT e.id, e.filename, e.sample_id, e.experiment_id,
               EXISTS (SELECT 1 FROM peak_reviews r WHERE r.exposure_id = e.id AND r.direct = 1) AS labelled,
               COALESCE((SELECT sc.state = 'form_factor' FROM sample_calls sc
                          WHERE sc.sample_id = e.sample_id), 0) AS ff
          FROM exposures e
         WHERE e.status IS NOT 'rejected' AND e.kind = 'file'
           AND ($(include_unlabelled ? 1 : 0) OR labelled OR ($(include_form_factor ? 1 : 0) AND ff))
         ORDER BY e.id"""))
    for r in rows
        try
            key = group_by === :experiment ? r.experiment_id : r.sample_id
            ismissing(key) && throw(SkipExposure(:no_group))
            q, I, σ = _load_label_trace(db, r)
            snap = _latest_snapshot(db, r.id)
            seen = snap === nothing ? [(q = Float64(x.q), found_by = "hybrid") for x in sqlrows(
                DBInterface.execute(db, "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY id", [r.id]))] :
                   [(q = Float64(p.q), found_by = String(p.found_by)) for p in snap]
            aq = Float64[p.q for p in seen]
            fl = curation_flags(db, r.id, aq)
            left = .!fl.excluded .& .!fl.superseded
            curated = r.labelled == 1 || r.ff == 1
            kept = [p.q for (p, ok) in zip(seen, left) if curated && ok && p.found_by != "model"]
            kept_model = [p.q for (p, ok) in zip(seen, left) if curated && ok && p.found_by == "model"]
            excl = curated ? Float64[x.q for x in sqlrows(DBInterface.execute(db,
                "SELECT q FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude' ORDER BY id", [r.id]))] : Float64[]
            push!(out, HimalayaPeakModel.LabeledExposure(r.id, Int(key), q, I, σ,
                curated ? Float64[a.q for a in fl.adds] : Float64[], kept, excl, r.labelled == 1, r.ff == 1, kept_model))
        catch err
            err isa SkipExposure || rethrow()
            skips[err.reason] = get(skips, err.reason, 0) + 1
        end
    end
    if !isempty(skips)
        println(io, "skipped $(sum(values(skips))) exposures: ",
                join(["$(n) ($(SKIP_TEXT[k]))" for (k, n) in sort!(collect(skips))], ", "))
    end
    out
end
