using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# Spec §3.1 — read-only `/api/resolve` route. Translates slugs → IDs in
# one round trip. Read-only — no with_idempotency, no SSE, no events.
#
# Helpers `_setup_analyzed_exposure` (defined in test_route_response_shapes.jl)
# and `with_test_server` (defined in test_http.jl) are available because
# runtests.jl includes those files before this one. See Step 3 for the
# correct ordering when adding the new include lines.

# `_setup_for_resolve` is defined in test_route_response_shapes.jl
# (alongside `_setup_analyzed_exposure`) so both Task 2 and Task 4 can
# call it. See Task 2 Step 3 for the helper body.

@testset "GET /api/resolve" begin
    @testset "200: experiment + sample + exposure happy path" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=S1&exposure=JC001-007")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.experiment_id == ctx.experiment_id
                @test body.experiment_name == "test-exp"
                @test body.sample_id == ctx.sample_id
                @test body.sample_name == "S1"
                @test body.exposure_id == ctx.exposure_id
                @test body.exposure_filename == "JC001-007"
            end
        end
    end

    @testset "200: experiment-only" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.experiment_id == ctx.experiment_id
                @test body.experiment_name == "test-exp"
                @test !haskey(body, :sample_id) || body.sample_id === nothing
            end
        end
    end

    @testset "200: id-form lookup returns names" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment_id=$(ctx.experiment_id)&sample_id=$(ctx.sample_id)")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.experiment_name == "test-exp"
                @test body.sample_name == "S1"
            end
        end
    end

    @testset "404: missing experiment" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=nope"; status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.error == "not_found"
                @test body.missing == "experiment"
                @test body.missing_value == "nope"
            end
        end
    end

    @testset "404: missing sample (experiment_resolved present)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=nope"; status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.error == "not_found"
                @test body.missing == "sample"
                @test body.missing_value == "nope"
                @test body.experiment_resolved.id == ctx.experiment_id
                @test body.experiment_resolved.name == "test-exp"
            end
        end
    end

    @testset "404: missing exposure (experiment_resolved + sample_resolved present)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=S1&exposure=nope"; status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.error == "not_found"
                @test body.missing == "exposure"
                @test body.missing_value == "nope"
                @test body.experiment_resolved.id == ctx.experiment_id
                @test body.sample_resolved.id == ctx.sample_id
            end
        end
    end

    @testset "400: malformed numeric param returns 400, not 500" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment_id=abc"; status_exception=false)
                @test r.status == 400
                body = JSON3.read(String(r.body))
                @test body.error == "invalid_id"
            end
        end
    end

    @testset "404: stale-name regression — rename experiment, old name 404s" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            DBInterface.execute(ctx.db, "UPDATE experiments SET name = 'test-exp-renamed' WHERE id = ?",
                                [ctx.experiment_id])
            with_test_server(ctx.db) do port, base
                r1 = HTTP.get("$base/api/resolve?experiment=test-exp"; status_exception=false)
                @test r1.status == 404
                r2 = HTTP.get("$base/api/resolve?experiment=test-exp-renamed")
                @test r2.status == 200
            end
        end
    end

    @testset "404: id-form for deleted sample (cold-mount Zustand staleness)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            # FK enforcement is on (open_db sets PRAGMA foreign_keys = ON);
            # exposures.sample_id REFERENCES samples(id) without ON DELETE.
            # Disable FKs around the synthetic delete to simulate the
            # "sample id no longer exists" cold-mount state without
            # tearing down dependent rows.
            DBInterface.execute(ctx.db, "PRAGMA foreign_keys = OFF")
            DBInterface.execute(ctx.db, "DELETE FROM samples WHERE id = ?", [ctx.sample_id])
            DBInterface.execute(ctx.db, "PRAGMA foreign_keys = ON")
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment_id=$(ctx.experiment_id)&sample_id=$(ctx.sample_id)";
                             status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.missing == "sample"
            end
        end
    end

    @testset "400: ambiguous params (experiment + experiment_id)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&experiment_id=$(ctx.experiment_id)";
                             status_exception=false)
                @test r.status == 400
                body = JSON3.read(String(r.body))
                @test body.error == "ambiguous_params"
            end
        end
    end

    @testset "400: ambiguous params (sample + sample_id)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample=S1&sample_id=$(ctx.sample_id)";
                             status_exception=false)
                @test r.status == 400
            end
        end
    end

    @testset "200: mixed name+id across entities is allowed" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            with_test_server(ctx.db) do port, base
                # name-form experiment + id-form sample is fine.
                r = HTTP.get("$base/api/resolve?experiment=test-exp&sample_id=$(ctx.sample_id)")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                @test body.experiment_id == ctx.experiment_id
                @test body.sample_id == ctx.sample_id
            end
        end
    end

    @testset "tiebreaker: duplicate experiment names → lowest id wins" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            # Insert a second experiment with the same name.
            res = DBInterface.execute(ctx.db,
                "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('test-exp','/p','/d','/a')")
            second_id = Int(DBInterface.lastrowid(res))
            @test second_id > ctx.experiment_id

            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment=test-exp")
                @test r.status == 200
                body = JSON3.read(String(r.body))
                # Lowest id wins deterministically.
                @test body.experiment_id == ctx.experiment_id
                @test body.experiment_id < second_id
            end
        end
    end

    @testset "404: id-form lookup of NULL-name sample (no canonical slug)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            DBInterface.execute(ctx.db, "UPDATE samples SET name = NULL WHERE id = ?",
                                [ctx.sample_id])
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment_id=$(ctx.experiment_id)&sample_id=$(ctx.sample_id)";
                             status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.missing == "sample"
                @test body.reason == "sample has no canonical name"
                # experiment_resolved present so the frontend StaleUrlPage
                # can offer "go to the experiment" recovery.
                @test body.experiment_resolved.id == ctx.experiment_id
            end
        end
    end

    @testset "404: id-form lookup of NULL-filename exposure (no canonical slug)" begin
        mktempdir() do tmp
            ctx = _setup_for_resolve(tmp)
            DBInterface.execute(ctx.db, "UPDATE exposures SET filename = NULL WHERE id = ?",
                                [ctx.exposure_id])
            with_test_server(ctx.db) do port, base
                r = HTTP.get("$base/api/resolve?experiment_id=$(ctx.experiment_id)&sample_id=$(ctx.sample_id)&exposure_id=$(ctx.exposure_id)";
                             status_exception=false)
                @test r.status == 404
                body = JSON3.read(String(r.body))
                @test body.missing == "exposure"
                @test body.reason == "exposure has no canonical filename"
                @test body.experiment_resolved.id == ctx.experiment_id
                @test body.sample_resolved.id == ctx.sample_id
            end
        end
    end
end
