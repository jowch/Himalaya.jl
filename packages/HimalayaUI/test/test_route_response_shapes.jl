using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# Route response-shape contract test (deep-scan follow-up).
#
# Catches the bug class found in the post-fourth-pass deep scan: a route's
# JSON response shape silently drifting out of alignment with its
# TypeScript counterpart. The Vitest unit tests' fetch mocks are
# hand-maintained — they encode the type, not the route — so a route
# that changes its keys won't fail any frontend test until production.
#
# Each subtest pins the EXACT set of top-level keys (and optionally key
# types) on the queue-affecting routes' responses. To change a response
# shape: update the test here AND the corresponding TypeScript interface
# in `frontend/src/api.ts` together. Anyone bumping one without the
# other will fail this test.

"""
Assert the `actual` JSON object's top-level keys exactly match `expected`.
A mismatch in either direction (extra or missing) is a contract change.
"""
function assert_keys(actual, expected::Vector{Symbol})
    actual_keys = Set(Symbol(string(k)) for k in keys(actual))
    expected_set = Set(expected)
    @test actual_keys == expected_set
end

# Fixture: a fully-analyzed exposure with peaks + indices. Reuses the
# example_tot.dat trace (same as test_routes_analysis.jl). Returns a tuple
# (db, exposure_id, sample_id, analysis_dir).
function _setup_analyzed_exposure(tmp::String)
    analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
    mkpath(analysis_dir)
    cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
       joinpath(analysis_dir, "example_tot.dat"))
    db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
    exp_id = HimalayaUI.init_experiment!(db; path=tmp,
        data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
    s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
    e_id   = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
    HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)
    (db = db, exposure_id = e_id, sample_id = s_id, analysis_dir = analysis_dir)
end

@testset "Route response shapes (queue contract)" begin

    @testset "POST /api/exposures/:id/peaks → PeakAddResponse (flat Peak & metadata)" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.post("$base/api/exposures/$(ctx.exposure_id)/peaks";
                    body = JSON3.write(Dict(:q => 0.5)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 201
                body = JSON3.read(String(r.body))
                # Full Peak shape (intensity/prominence/sharpness null on
                # manual peaks — but explicitly present in the response, NOT
                # omitted) plus metadata. PeakAddResponse extends Peak in
                # api.ts so missing fields would silently leave the cache
                # entry with `undefined` for those columns. Issue #17 made
                # this contract failure explicit; this test pins it.
                assert_keys(body, [
                    :id, :exposure_id, :q,
                    :intensity, :prominence, :sharpness,
                    :source, :excluded,
                    :event_id, :view_row_id, :analysis_inputs_hash,
                ])
                @test body.id isa Integer
                @test body.q == 0.5
                @test body.source == "manual"
                @test body.excluded === false
                @test body.intensity === nothing
                @test body.prominence === nothing
                @test body.sharpness === nothing
                @test body.event_id isa Integer
                @test body.view_row_id isa Integer
                @test body.analysis_inputs_hash isa AbstractString
            end
        end
    end

    @testset "PATCH /api/peaks/:id → PeakUpdatedResponse (extends Peak)" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            auto = first(Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT id FROM auto_peaks WHERE exposure_id = ? LIMIT 1",
                [ctx.exposure_id])))
            with_test_server(ctx.db) do port, base
                r = HTTP.patch("$base/api/peaks/$(auto.id)";
                    body = JSON3.write(Dict(:excluded => true)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 200
                body = JSON3.read(String(r.body))
                # Same shape as PeakAddResponse: Peak fields + metadata.
                # All Peak interface fields should be present (intensity,
                # prominence, sharpness exist on auto peaks but null on
                # manual peaks; here we patched an auto peak so non-null).
                expected = [
                    :id, :exposure_id, :q, :intensity, :prominence, :sharpness,
                    :source, :excluded,
                    :event_id, :view_row_id, :analysis_inputs_hash,
                ]
                assert_keys(body, expected)
                @test body.excluded === true
            end
        end
    end

    @testset "DELETE /api/peaks/:id (manual) → PeakRemoveResponse" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                # Add a manual peak first so we have something deletable.
                r = HTTP.post("$base/api/exposures/$(ctx.exposure_id)/peaks";
                    body = JSON3.write(Dict(:q => 0.5)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                manual_id = JSON3.read(String(r.body)).id

                r2 = HTTP.delete("$base/api/peaks/$manual_id";
                    headers = ["X-Username" => "alice"])
                @test r2.status == 200
                body = JSON3.read(String(r2.body))
                assert_keys(body, [:event_id, :view_row_id, :analysis_inputs_hash])
                @test body.event_id isa Integer
                @test body.analysis_inputs_hash isa AbstractString
            end
        end
    end

    @testset "POST /api/exposures/:id/analyze → ReanalyzeResponse" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Force a non-no-op so the slow path runs.
            DBInterface.execute(ctx.db,
                "UPDATE exposures SET analysis_inputs_hash = NULL WHERE id = ?",
                [ctx.exposure_id])
            with_test_server(ctx.db) do port, base
                r = HTTP.post("$base/api/exposures/$(ctx.exposure_id)/analyze";
                    headers = ["X-Username"     => "alice",
                               "X-Client-Op-Id" => "uuid-shape-reanalyze"])
                @test r.status == 200
                body = JSON3.read(String(r.body))
                assert_keys(body, [:id, :analyzed, :analysis_inputs_hash])
                @test body.analyzed === true
                @test body.analysis_inputs_hash isa AbstractString
            end
        end
    end

    @testset "POST /api/groups/:id/members → GroupMutationResponse" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.get("$base/api/exposures/$(ctx.exposure_id)/groups")
                groups = JSON3.read(String(r1.body))
                gid = groups[1].id

                r2 = HTTP.get("$base/api/exposures/$(ctx.exposure_id)/indices")
                indices = JSON3.read(String(r2.body))
                idx = first(filter(i -> !(i.id in groups[1].members), indices))

                r3 = HTTP.post("$base/api/groups/$gid/members";
                    body = JSON3.write(Dict(:index_id => idx.id)),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r3.status == 200
                body = JSON3.read(String(r3.body))
                # GroupMutationResponse = GroupEntry & {event_id, view_row_id}
                # GroupEntry: {id, exposure_id, kind, active, members}
                assert_keys(body, [
                    :id, :exposure_id, :kind, :active, :members,
                    :event_id, :view_row_id,
                ])
                @test body.event_id isa Integer
            end
        end
    end

    @testset "DELETE /api/groups/:id/members/:idx → GroupMutationResponse" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r1 = HTTP.get("$base/api/exposures/$(ctx.exposure_id)/groups")
                groups = JSON3.read(String(r1.body))
                gid = groups[1].id
                isempty(groups[1].members) && return
                idx_id = first(groups[1].members)

                r2 = HTTP.delete("$base/api/groups/$gid/members/$idx_id";
                    headers = ["X-Username" => "alice"])
                @test r2.status == 200
                body = JSON3.read(String(r2.body))
                assert_keys(body, [
                    :id, :exposure_id, :kind, :active, :members,
                    :event_id, :view_row_id,
                ])
            end
        end
    end

    @testset "GET /api/exposures/:id/indices → IndexEntry[] full shape" begin
        # Suggestion #11: pin the GET /indices shape so a future SELECT
        # change can't silently leak a server-internal column into every
        # cache write or SSE post_state frame.
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/exposures/$(ctx.exposure_id)/indices")
                @test r.status == 200
                indices = JSON3.read(String(r.body))
                @test !isempty(indices)
                expected = [
                    :id, :exposure_id, :phase, :basis, :score, :r_squared,
                    :lattice_d, :ngc, :status, :kind, :inputs_hash,
                    :peaks, :predicted_q,
                ]
                # Every entry must have exactly these top-level keys —
                # IndexEntry is closed.
                for ix in indices
                    assert_keys(ix, expected)
                end
            end
        end
    end

    @testset "GET /api/indices/:id → single IndexEntry full shape" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                rs = HTTP.get("$base/api/exposures/$(ctx.exposure_id)/indices")
                first_id = JSON3.read(String(rs.body))[1].id
                r = HTTP.get("$base/api/indices/$first_id")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                expected = [
                    :id, :exposure_id, :phase, :basis, :score, :r_squared,
                    :lattice_d, :ngc, :status, :kind, :inputs_hash,
                    :peaks, :predicted_q,
                ]
                assert_keys(body, expected)
            end
        end
    end

    @testset "POST /api/exposures/:id/speculative → IndexEntry full shape (incl. ngc)" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            peaks = Tables.rowtable(DBInterface.execute(ctx.db,
                "SELECT id FROM auto_peaks WHERE exposure_id = ? ORDER BY q LIMIT 2",
                [ctx.exposure_id]))
            p1 = Int(peaks[1].id); p2 = Int(peaks[2].id)
            with_test_server(ctx.db) do port, base
                body_in = Dict(:phase => "Pn3m",
                               :anchor_peak_id => p1, :anchor_ratio => 1,
                               :additional => [Dict(:ratio_position => 2, :peak_id => p2)])
                r = HTTP.post("$base/api/exposures/$(ctx.exposure_id)/speculative";
                    body = JSON3.write(body_in),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 200
                body = JSON3.read(String(r.body))
                # IndexEntry interface in api.ts: id, exposure_id, phase,
                # basis, score, r_squared, lattice_d, ngc, status, kind,
                # inputs_hash, peaks, predicted_q.
                # `ngc` was missing pre-deep-scan fix; pin it explicitly.
                expected = [
                    :id, :exposure_id, :phase, :basis, :score, :r_squared,
                    :lattice_d, :ngc, :status, :kind, :inputs_hash,
                    :peaks, :predicted_q,
                ]
                assert_keys(body, expected)
                @test body.kind == "speculative"
            end
        end
    end

    @testset "PATCH /api/samples/:id → bare samples row (NO tags field — known)" begin
        # Pinning the fact that this route does NOT include tags. The
        # frontend `updateSampleMutator` knows this and merges only the
        # patched fields into the existing cache entry to preserve tags.
        # If a future change starts including tags, update the mutator
        # AND this test together.
        mktempdir() do tmp
            db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
            with_test_server(db) do port, base
                r = HTTP.patch("$base/api/samples/$s_id";
                    body = JSON3.write(Dict(:notes => "n")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test :tags ∉ keys(body)
                @test :id in keys(body)
                @test :notes in keys(body)
            end
        end
    end

    @testset "POST /api/samples/:id/tags → SampleTag-shaped (with sample_id allowed extra)" begin
        mktempdir() do tmp
            db     = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
            s_id   = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
            with_test_server(db) do port, base
                r = HTTP.post("$base/api/samples/$s_id/tags";
                    body = JSON3.write(Dict(:key => "k", :value => "v")),
                    headers = ["Content-Type" => "application/json",
                               "X-Username"   => "alice"])
                @test r.status == 201
                body = JSON3.read(String(r.body))
                # Frontend SampleTag interface = {id, key, value, source}.
                # The route includes an extra `sample_id` (and `exposure_id`
                # for the exposure variant). They're harmless cache
                # pollution today; documented here so a future tightening
                # of the shape is a deliberate change.
                for required in (:id, :key, :value, :source)
                    @test required in keys(body)
                end
            end
        end
    end
end
