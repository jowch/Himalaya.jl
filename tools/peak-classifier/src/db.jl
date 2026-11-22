using SQLite, DBInterface, Tables

"""
    load_labeled_exposures(db_path; group_by = :experiment, include_form_factor = true)
        -> Vector{LabeledExposure}

Reads curated exposures (and form-factor exposures) from a HimalayaUI
database. The file is opened read-only without migration. Curations are
matched to auto peaks with `HimalayaUI.curation_flags`, the same code the app
uses. Exposures whose trace file is missing or unreadable are skipped and
counted on stderr.
"""
function load_labeled_exposures(db_path::AbstractString; group_by::Symbol = :experiment,
                                include_form_factor::Bool = true)
    group_by in (:experiment, :sample) || throw(ArgumentError("group_by must be :experiment or :sample"))
    db = SQLite.DB(db_path)
    DBInterface.execute(db, "PRAGMA query_only = ON")
    rows = Tables.rowtable(DBInterface.execute(db, """
        SELECT e.id, e.filename, e.sample_id, e.experiment_id,
               EXISTS(SELECT 1 FROM peak_curations c WHERE c.exposure_id = e.id) AS curated,
               COALESCE((SELECT sc.state = 'form_factor' FROM sample_calls sc
                          WHERE sc.sample_id = e.sample_id), 0) AS ff
          FROM exposures e
         WHERE e.status IS NOT 'rejected' AND e.kind = 'file'
           AND (curated OR ($(include_form_factor ? 1 : 0) AND ff))
         ORDER BY e.id"""))
    out = LabeledExposure[]
    skipped = 0
    for r in rows
        try
            setup = HimalayaUI.exposure_setup(db, r.id)
            setup === nothing && throw(ErrorException("no setup"))
            path = HimalayaUI.resolve_trace_path(String(setup.analysis_dir),
                HimalayaUI.integration_pattern_for(setup), String(r.filename))
            path === nothing && throw(ErrorException("no trace file"))
            q, I, σ = HimalayaUI.load_dat(path)
            auto = Tables.rowtable(DBInterface.execute(db,
                "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY id", [r.id]))
            aq = Float64[x.q for x in auto]
            fl = HimalayaUI.curation_flags(db, r.id, aq)
            kept = aq[.!fl.excluded .& .!fl.superseded]
            push!(out, LabeledExposure(r.id, group_by === :experiment ? r.experiment_id : r.sample_id,
                q, I, σ, Float64[a.q for a in fl.adds], kept, aq[fl.excluded],
                r.curated == 1, r.ff == 1))
        catch err
            skipped += 1
        end
    end
    close(db)
    skipped > 0 && println(stderr, "skipped $skipped exposures (trace file missing or unreadable)")
    out
end
