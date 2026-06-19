# Shared cross-file test helpers for the HimalayaUI test suite.
#
# Extracted from test_comparisons.jl and test_route_response_shapes.jl so that
# future parallel GROUP buckets can each include this file once without hitting
# UndefVarError from severed serial-include ordering.
#
# Canonical rule: ONE definition per helper here. Callers that previously held
# a local copy must use this one instead.
#
# Do NOT add @testset blocks here — this file is included once per bucket,
# not run as a test file itself.

using DBInterface, Tables
using HimalayaUI

# Pre-mint a comparisons row at a known id so the dispatcher can run UPDATEs
# (matches the route handler's two-step "INSERT placeholder, then dispatcher
# fills in" pattern). Mirrors `_premint_comparison!` in test_events.jl.
# Uses NULL placeholders to mirror the post-#67 route — the dispatcher's
# `COALESCE(col, ?)` then stamps real values on first fold.
function _premint_cmp!(db, id::Int)
    DBInterface.execute(db,
        """INSERT INTO comparisons (id, title, content_hash, created_at, updated_at)
           VALUES (?, NULL, NULL, NULL, NULL)""", [id])
    nothing
end

function _member_payload(; id=nothing, exposure_id=nothing, display_order::Int=0,
                          band_height::Float64=1.0, y_offset::Float64=0.0,
                          normalization::String="none",
                          color_override=nothing, label_override=nothing,
                          q_window_min=nothing, q_window_max=nothing,
                          peak_display=nothing,
                          snapshot=Dict(:effective_peaks => [],
                                        :confirmed_index => nothing,
                                        :analysis_inputs_hash => "sha256:zero"))
    Dict{Symbol,Any}(
        :id             => id,
        :exposure_id    => exposure_id,
        :display_order  => display_order,
        :band_height    => band_height,
        :y_offset       => y_offset,
        :normalization  => normalization,
        :color_override => color_override,
        :label_override => label_override,
        :q_window_min   => q_window_min,
        :q_window_max   => q_window_max,
        :peak_display   => peak_display,
        :snapshot       => snapshot,
    )
end

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
    e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename=filename)
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
