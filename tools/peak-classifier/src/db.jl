using SQLite, DBInterface, Tables

# An exposure that cannot be loaded for a known reason. Any other exception is
# a bug or an incompatible database and propagates.
struct SkipExposure <: Exception
    reason::Symbol
end

const SKIP_TEXT = Dict(
    :no_group => "no experiment or sample id",
    :no_setup => "no setup",
    :no_trace_file => "trace file not found",
    :unreadable_trace => "trace file not parseable")

# Percent-encode the characters that are special in a SQLite URI path.
_uri_path(p) = replace(abspath(p), "%" => "%25", "?" => "%3f", "#" => "%23")

"""
    open_readonly(db_path) -> SQLite.DB

Opens with the URI `mode=ro`, so a missing file is an error rather than a new
empty database, and closing the connection never checkpoints a `-wal` file.
`PRAGMA query_only` is a second guard. A WAL database still needs its `-shm`
file or a writable directory to be opened read-only.
"""
function open_readonly(db_path::AbstractString)
    db = SQLite.DB("file:" * _uri_path(db_path) * "?mode=ro")
    try
        DBInterface.execute(db, "PRAGMA query_only = ON")
    catch
        close(db)
        rethrow()
    end
    db
end

function _load_trace(db, r)
    setup = HimalayaUI.exposure_setup(db, r.id)
    setup === nothing && throw(SkipExposure(:no_setup))
    path = HimalayaUI.resolve_trace_path(String(setup.analysis_dir),
        HimalayaUI.integration_pattern_for(setup), String(r.filename))
    path === nothing && throw(SkipExposure(:no_trace_file))
    try
        HimalayaUI.load_dat(path)
    catch err
        # readdlm and load_dat signal a malformed file with these types only
        err isa Union{ErrorException, ArgumentError, DimensionMismatch, BoundsError} || rethrow()
        throw(SkipExposure(:unreadable_trace))
    end
end

"""
    load_labeled_exposures(db_path; group_by = :experiment, include_form_factor = true,
                           io = stderr) -> Vector{LabeledExposure}

Reads curated exposures (and form-factor exposures) from a HimalayaUI
database. The file is opened read-only (`open_readonly`) without migration.
Curations are matched to auto peaks with `HimalayaUI.curation_flags`, the same
code the app uses. Exposures with no group id, no setup, or a missing or
unparseable trace file are skipped and counted by reason on `io`; any other
error propagates.
"""
function load_labeled_exposures(db_path::AbstractString; group_by::Symbol = :experiment,
                                include_form_factor::Bool = true, io::IO = stderr)
    group_by in (:experiment, :sample) || throw(ArgumentError("group_by must be :experiment or :sample"))
    db = open_readonly(db_path)
    out = LabeledExposure[]
    skips = Dict{Symbol, Int}()
    try
        rows = Tables.rowtable(DBInterface.execute(db, """
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
                q, I, σ = _load_trace(db, r)
                auto = Tables.rowtable(DBInterface.execute(db,
                    "SELECT q FROM auto_peaks WHERE exposure_id = ? ORDER BY id", [r.id]))
                aq = Float64[x.q for x in auto]
                fl = HimalayaUI.curation_flags(db, r.id, aq)
                kept = aq[.!fl.excluded .& .!fl.superseded]
                push!(out, LabeledExposure(r.id, Int(key), q, I, σ, Float64[a.q for a in fl.adds],
                    kept, aq[fl.excluded], r.curated == 1, r.ff == 1))
            catch err
                err isa SkipExposure || rethrow()
                skips[err.reason] = get(skips, err.reason, 0) + 1
            end
        end
    finally
        close(db)
    end
    if !isempty(skips)
        println(io, "skipped $(sum(values(skips))) exposures: ",
                join(["$(n) ($(SKIP_TEXT[k]))" for (k, n) in sort!(collect(skips))], ", "))
    end
    out
end
