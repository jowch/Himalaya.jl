# Shared cross-file test helpers for the HimalayaUI test suite.
#
# Kept in one file so that parallel GROUP buckets can each include this file once without hitting
# UndefVarError from severed serial-include ordering.
#
# Canonical rule: ONE definition per helper here. Callers that previously held
# a local copy must use this one instead.
#
# Do NOT add @testset blocks here — this file is included once per bucket,
# not run as a test file itself.

using DBInterface, Tables
using HimalayaUI

# Build an analyzed exposure. Returns a NamedTuple:
#   (db, experiment_id, sample_id, exposure_id, analysis_dir)
#
# This is the 5-field canonical form (superset of the 4-field version that
# used to live in test_route_response_shapes.jl). Callers that only need a
# subset of these fields can safely ignore the extras — a NamedTuple with
# an extra field is backward-compatible.
function _setup_analyzed_exposure(tmp::String;
                                   datfile::String = "example_tot.dat",
                                   filename::String = "example_tot")
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", datfile),
       joinpath(analysis_dir, datfile))
    db = open_prepared_clone(tmp)
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id = HimalayaUI.create_exposure!(db; experiment_id=exp_id, sample_id=s_id, filename=filename)
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    (db = db, experiment_id = exp_id, sample_id = s_id,
     exposure_id = e_id, analysis_dir = analysis_dir)
end

"""
Like `_setup_analyzed_exposure` but UPDATEs the rows to known slug-resolvable
names ("test-exp" / "S1" / "JC001-007") and captures `experiment_id` for tests
that need it. Used by `test_routes_resolve.jl` and the resolve-shape rows in
`test_route_response_shapes.jl`.
"""
function _setup_for_resolve(tmp::String)
    ctx = _setup_analyzed_exposure(tmp)
    DBInterface.execute(ctx.db, "UPDATE experiments SET name = 'test-exp'")
    DBInterface.execute(ctx.db, "UPDATE samples SET name = 'S1' WHERE id = ?",
                        [ctx.sample_id])
    DBInterface.execute(ctx.db, "UPDATE exposures SET filename = 'JC001-007' WHERE id = ?",
                        [ctx.exposure_id])
    exp_row = Tables.rowtable(DBInterface.execute(ctx.db,
        "SELECT id FROM experiments LIMIT 1"))[1]
    return (db = ctx.db,
            experiment_id = Int(exp_row.id),
            sample_id = ctx.sample_id,
            exposure_id = ctx.exposure_id)
end

function align_next_curation_id!(db, id::Integer)
    taken = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT COALESCE(MAX(id), 0) AS m FROM peak_curations"))).m
    taken < id || error("peak_curations already holds id $(taken) ≥ $(id)")
    DBInterface.execute(db, "DELETE FROM sqlite_sequence WHERE name = 'peak_curations'")
    DBInterface.execute(db,
        "INSERT INTO sqlite_sequence (name, seq) VALUES ('peak_curations', ?)", [id - 1])
    nothing
end
