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

function labeled_exposures(db::SQLite.DB; group_by::Symbol = :experiment,
                           include_form_factor::Bool = true, io::IO = stderr)
    group_by in (:experiment, :sample) || throw(ArgumentError("group_by must be :experiment or :sample"))
    out = HimalayaPeakModel.LabeledExposure[]
    skips = Dict{Symbol, Int}()
    rows = sqlrows(DBInterface.execute(db, """
        SELECT e.id, e.filename, e.sample_id, e.experiment_id,
               EXISTS(SELECT 1 FROM peak_curations c WHERE c.exposure_id = e.id) AS curated,
               COALESCE((SELECT sc.state = 'form_factor' FROM sample_calls sc
                          WHERE sc.sample_id = e.sample_id), 0) AS ff
          FROM exposures e
         WHERE e.status IS NOT 'rejected' AND e.kind = 'file'
           AND (curated OR ($(include_form_factor ? 1 : 0) AND ff))
         ORDER BY e.id"""))
    for r in rows
        try
            key = group_by === :experiment ? r.experiment_id : r.sample_id
            ismissing(key) && throw(SkipExposure(:no_group))
            q, I, σ = _load_label_trace(db, r)
            auto = sqlrows(DBInterface.execute(db,
                "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY id", [r.id]))
            aq = Float64[x.q for x in auto]
            fl = curation_flags(db, r.id, aq)
            kept = aq[.!fl.excluded .& .!fl.superseded]
            excl = Float64[x.q for x in sqlrows(DBInterface.execute(db,
                "SELECT q FROM peak_curations WHERE exposure_id = ? AND kind = 'exclude' ORDER BY id", [r.id]))]
            push!(out, HimalayaPeakModel.LabeledExposure(r.id, Int(key), q, I, σ,
                Float64[a.q for a in fl.adds], kept, excl, r.curated == 1, r.ff == 1))
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
