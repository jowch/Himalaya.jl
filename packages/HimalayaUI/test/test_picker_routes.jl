using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# ─────────────────────────────────────────────────────────────────────────────
# Phase 5 Task 5.2 — picker support routes.
#
# Three GET endpoints feed the picker modal with a "Recently used" section
# and a tag filter dropdown:
#   GET /api/users/:id/recently-picked-exposures?limit=20
#   GET /api/experiments/:eid/sample-tags
#   GET /api/sample-tags              (corpus-wide sibling of the route above)
#
# All are read-only — no with_idempotency, no SSE, no event log row.
# ─────────────────────────────────────────────────────────────────────────────

@testset "picker routes" begin
    @testset "GET /api/users/:id/recently-picked-exposures: ordered list" begin
        mktempdir() do tmp
            ctx   = _setup_analyzed_exposure(tmp)
            alice = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            e2    = HimalayaUI.create_exposure!(ctx.db; sample_id=ctx.sample_id, filename="e2")
            e3    = HimalayaUI.create_exposure!(ctx.db; sample_id=ctx.sample_id, filename="e3")
            _premint_cmp!(ctx.db, 200)
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(ctx.db, req;
                kind="comparison_created", entity_type="comparison", entity_id=200,
                payload=Dict(:title => "T", :description => nothing,
                             :forked_from_id => nothing, :forked_at_hash => nothing,
                             :members => [
                                _member_payload(exposure_id=ctx.exposure_id, display_order=0),
                                _member_payload(exposure_id=e2, display_order=1),
                                _member_payload(exposure_id=e3, display_order=2),
                             ]))
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/users/$alice/recently-picked-exposures")
                @test r.status == 200
                ids = JSON3.read(String(r.body))
                # All three exposures should be present.
                @test sort(collect(ids)) == sort([ctx.exposure_id, e2, e3])

                # Limit honored.
                r = HTTP.get("$base/api/users/$alice/recently-picked-exposures?limit=2")
                @test r.status == 200
                ids2 = JSON3.read(String(r.body))
                @test length(ids2) == 2
            end
        end
    end

    @testset "GET /api/users/:id/recently-picked-exposures: empty list when no membership" begin
        mktempdir() do tmp
            ctx   = _setup_analyzed_exposure(tmp)
            alice = HimalayaUI.get_or_create_user!(ctx.db, "alice")
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/users/$alice/recently-picked-exposures")
                @test r.status == 200
                @test JSON3.read(String(r.body)) == []
            end
        end
    end

    @testset "GET /api/users/:id/recently-picked-exposures: 404 for unknown user id" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/users/9999/recently-picked-exposures";
                             status_exception = false)
                @test r.status == 404
            end
        end
    end

    @testset "GET /api/experiments/:eid/sample-tags: distinct (key, value) pairs" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Seed two more samples so we cover multi-sample aggregation, and
            # populate sample_tags rows. Critical: two tags share the same VALUE
            # but different KEYS — both must surface as separate filter options.
            s2 = HimalayaUI.create_sample!(ctx.db; experiment_id=ctx.experiment_id, name="D2")
            s3 = HimalayaUI.create_sample!(ctx.db; experiment_id=ctx.experiment_id, name="D3")
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [ctx.sample_id, "lipid", "DOPC"])
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [s2, "lipid", "DOPE"])
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [s3, "lipid", "DOPC"]) # duplicate (key, value) → distinct collapses
            # Same VALUE "DOPC" but a different KEY "control" — must NOT collapse.
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [s3, "control", "DOPC"])

            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/experiments/$(ctx.experiment_id)/sample-tags")
                @test r.status == 200
                tags = JSON3.read(String(r.body))
                pairs = Set([(String(t.key), String(t.value)) for t in tags])
                @test ("lipid", "DOPC") in pairs
                @test ("lipid", "DOPE") in pairs
                # Critical: same value, different key → both present as
                # separate options. Distinct collapses on the (key, value)
                # pair, not on value alone.
                @test ("control", "DOPC") in pairs
                # Distinct: the duplicate (lipid, DOPC) on s3 collapses.
                @test length(tags) == 3
            end
        end
    end

    @testset "GET /api/experiments/:eid/sample-tags: only tags in scope" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Build a second experiment with its own sample + tag — must not
            # leak into the first experiment's sample-tags response.
            e2_id = HimalayaUI.init_experiment!(ctx.db; path=tmp * "/e2",
                data_dir=tmp * "/e2/data", analysis_dir=tmp * "/e2/analysis")
            s_other = HimalayaUI.create_sample!(ctx.db; experiment_id=e2_id, name="OTHER")
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [s_other, "leaky", "value"])

            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [ctx.sample_id, "lipid", "DOPC"])

            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/experiments/$(ctx.experiment_id)/sample-tags")
                @test r.status == 200
                tags = JSON3.read(String(r.body))
                keys_seen = Set(String(t.key) for t in tags)
                @test "lipid" in keys_seen
                @test !("leaky" in keys_seen)
            end
        end
    end

    @testset "GET /api/experiments/:eid/sample-tags: empty list when no tags" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/experiments/$(ctx.experiment_id)/sample-tags")
                @test r.status == 200
                @test JSON3.read(String(r.body)) == []
            end
        end
    end

    @testset "GET /api/sample-tags: corpus-wide tags across experiments" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            # Second experiment with its own sample + tag. The corpus route
            # MUST include it — this is the inverse of the experiment-scoped
            # route's "only tags in scope" test.
            e2_id = HimalayaUI.init_experiment!(ctx.db; path=tmp * "/e2",
                data_dir=tmp * "/e2/data", analysis_dir=tmp * "/e2/analysis")
            s_other = HimalayaUI.create_sample!(ctx.db; experiment_id=e2_id, name="OTHER")
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [ctx.sample_id, "lipid", "DOPC"])
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [s_other, "buffer", "PBS"])
            # Same VALUE "DOPC" as the lipid tag above but a different KEY.
            # DISTINCT collapses on the (key, value) PAIR, not on value alone,
            # so this must surface as its own corpus entry.
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [s_other, "buffer", "DOPC"])
            # A duplicate (key, value) on a third sample in experiment 1 —
            # DISTINCT must collapse it to a single corpus entry.
            s3 = HimalayaUI.create_sample!(ctx.db; experiment_id=ctx.experiment_id, name="D3")
            DBInterface.execute(ctx.db,
                "INSERT INTO sample_tags (sample_id, key, value, source) VALUES (?, ?, ?, 'manual')",
                [s3, "lipid", "DOPC"])

            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/sample-tags")
                @test r.status == 200
                tags = JSON3.read(String(r.body))
                pairs = Set([(String(t.key), String(t.value)) for t in tags])
                # Tags from BOTH experiments are present.
                @test ("lipid", "DOPC") in pairs
                @test ("buffer", "PBS") in pairs
                # Same value, different key → surfaces separately; DISTINCT is
                # on the (key, value) pair, matching the docstring contract.
                @test ("buffer", "DOPC") in pairs
                # DISTINCT: the duplicate (lipid, DOPC) collapses to one entry.
                @test length(tags) == 3
            end
        end
    end

    @testset "GET /api/sample-tags: empty list when no tags" begin
        mktempdir() do tmp
            ctx = _setup_analyzed_exposure(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/sample-tags")
                @test r.status == 200
                @test JSON3.read(String(r.body)) == []
            end
        end
    end
end
