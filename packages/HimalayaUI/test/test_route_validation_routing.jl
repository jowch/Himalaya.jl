using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# Failure-class routing parametric test (gap-1 from review-discipline).
#
# The frontend's failure-class router (lib/queue/errors.ts +
# InfrastructureBanner.tsx) routes 4xx → validation toast, 5xx →
# InfrastructureBanner. Per spec, MALFORMED-BODY errors must surface
# as toasts (the user typed something wrong), not as the banner (the
# server is broken). Smoke Bug B was exactly this — `Float64(body.q)`
# threw an unguarded MethodError → 500 → mis-routed to the banner.
#
# This testset hits every mutating route that accepts a JSON body
# with a malformed/missing-field body and asserts the response status
# is in [400, 500). A future regression that re-introduces an
# un-guarded conversion throw will fail here.
#
# Routes that don't accept a body (DELETE, PATCH /select) are listed
# at the bottom for completeness — they can't have malformed-body
# regressions but their handling-of-bad-id is a separate concern.

function _setup_full_fixture(tmp::String)
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

    # Ensure there's at least one auto-peak so PATCH /peaks/:id has a
    # real id to target.
    auto_id = first(Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM auto_peaks WHERE exposure_id = ? LIMIT 1", [e_id]))).id

    (db = db, experiment_id = exp_id, sample_id = s_id, exposure_id = e_id,
     auto_peak_id = Int(auto_id), analysis_dir = analysis_dir)
end

# Bodies designed to be malformed in the route's required-field sense.
# Each route's body field expectations live near the route handler; the
# malformed body should miss/scramble the required keys WITHOUT being
# JSON-syntactically broken (we want to test the route's validation,
# not the JSON parser).
@testset "Validation routing: every mutating route returns 4xx on malformed body" begin
    mktempdir() do tmp
        ctx = _setup_full_fixture(tmp)
        with_inproc_routes(ctx.db) do call
            base_headers = ["Content-Type" => "application/json",
                            "X-Username"   => "alice"]

            cases = [
                # (method, path, body_dict, label)
                ("POST",  "/api/users",
                 Dict(:NOT_username => "x"), "POST /users — missing username"),
                ("PATCH", "/api/experiments/$(ctx.experiment_id)",
                 Dict(:NOT_a_field => "x"), "PATCH /experiments — no patchable fields"),
                ("PATCH", "/api/samples/$(ctx.sample_id)",
                 Dict(:NOT_a_field => "x"), "PATCH /samples — no patchable fields"),
                ("POST",  "/api/samples/$(ctx.sample_id)/tags",
                 Dict(:value => "v"), "POST /samples/:id/tags — missing key"),
                ("POST",  "/api/samples/tags/batch",
                 Dict(:tags => [Dict(:sample_id => ctx.sample_id, :value => "v")]),
                 "POST /samples/tags/batch — missing key"),
                ("POST",  "/api/samples/$(ctx.sample_id)/messages",
                 Dict(:NOT_body => "hi"), "POST /samples/:id/messages — missing body"),
                ("POST",  "/api/samples/$(ctx.sample_id)/messages",
                 Dict(:body => "   "), "POST /samples/:id/messages — empty body"),
                ("POST",  "/api/exposures/$(ctx.exposure_id)/tags",
                 Dict(:value => "v"), "POST /exposures/:id/tags — missing key"),
                ("PATCH", "/api/exposures/$(ctx.exposure_id)/status",
                 Dict(:status => "not-a-valid-status"),
                 "PATCH /exposures/:id/status — invalid enum"),
                ("POST",  "/api/exposures/$(ctx.exposure_id)/peaks",
                 Dict(:NOT_q => 0.5), "POST /exposures/:id/peaks — missing q"),
                ("POST",  "/api/exposures/$(ctx.exposure_id)/peaks",
                 Dict(:q => "not-a-number"),
                 "POST /exposures/:id/peaks — non-numeric q"),
                ("PATCH", "/api/peaks/$(ctx.auto_peak_id)",
                 Dict(:NOT_excluded => true),
                 "PATCH /peaks/:id — missing excluded"),
                ("POST",  "/api/exposures/$(ctx.exposure_id)/speculative",
                 Dict(:phase => "Pn3m"),
                 "POST /speculative — missing anchor_peak_id"),
                ("POST",  "/api/exposures/$(ctx.exposure_id)/speculative",
                 Dict(:phase => "Bogus", :anchor_peak_id => 1, :anchor_ratio => 1),
                 "POST /speculative — unknown phase"),
            ]

            for (method, path, body, label) in cases
                @testset "$label" begin
                    r = call(method, path;
                        headers = base_headers,
                        body = Vector{UInt8}(JSON3.write(body)))
                    @test 400 <= r.status < 500
                end
            end

            # D-10: POST /groups/:id/members (and its missing-index_id validation)
            # was retired. The assignment-native POST /assignment/members carries
            # the equivalent validation, exercised in test_assignments.jl.
        end
    end
end
