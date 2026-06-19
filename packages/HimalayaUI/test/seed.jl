# seed.jl — shared test fixture helper.
#
# `cli_init_with_db!` (the manifest-driven CLI ingest) was deleted by the
# ingestion redesign; the HTTP scan path is now the sole ingestion entry point.
# Several non-CLI tests used `cli_init_with_db!` purely as SETUP — to stand up a
# DB with one experiment + sample + exposures pointing at on-disk .dat fixtures —
# and then exercised something else (analyze, fast-skip, _resolve_experiment).
#
# `seed_experiment!` reproduces that DB shape directly via the Phase-A writers
# (create_experiment! / create_sample! / create_exposure!), with no manifest and
# no scan. Callers that need the trace bytes on disk (to run analyze_exposure!)
# copy the fixture .dat into analysis_dir themselves, exactly as before — this
# helper only writes DB rows.

using HimalayaUI: create_experiment!, create_sample!, create_exposure!
using SQLite, DBInterface, Tables

"""
    seed_experiment!(db, exp_dir; name, analysis_dir, data_dir, sample_name, stems,
                     config=nothing, experiment_type=nothing) -> (exp_id, sample_id, exposure_ids)

Register one experiment (rooted at `exp_dir`) with a single sample owning one
exposure per stem in `stems`. Each exposure's `filename` is the bare stem (no
extension) — matching what the old manifest-driven ingest stored, so downstream
filename-based lookups (`WHERE filename = stem`) and pattern resolution
(`{name}.dat`) keep working.

Pure DB writes — does NOT touch the filesystem. Callers that subsequently run
`analyze_exposure!` must place the fixture `.dat` files under `analysis_dir`.
"""
function seed_experiment!(db::SQLite.DB, exp_dir::AbstractString;
        name::AbstractString,
        analysis_dir::AbstractString,
        data_dir::AbstractString = joinpath(exp_dir, "data"),
        sample_name::AbstractString = "D1",
        stems::AbstractVector = String["ST001"],
        config = nothing,
        experiment_type = nothing)
    exp_id = create_experiment!(db;
        name            = String(name),
        path            = String(exp_dir),
        data_dir        = String(data_dir),
        analysis_dir    = String(analysis_dir),
        config          = config,
        experiment_type = experiment_type)
    s_id = create_sample!(db; experiment_id = exp_id, name = String(sample_name))
    exposure_ids = Int[]
    for stem in stems
        e_id = create_exposure!(db;
            experiment_id = exp_id,
            sample_id     = s_id,
            filename      = String(stem))
        push!(exposure_ids, e_id)
    end
    return (exp_id = exp_id, sample_id = s_id, exposure_ids = exposure_ids)
end
