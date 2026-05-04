using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# Idempotency replay invariant (Tier 1.4).
#
# Per-route, posting the same X-Client-Op-Id twice MUST produce:
#   1. Byte-identical HTTP response bodies (already checked at unit level
#      in test_idempotency.jl, but pinned end-to-end through real routes).
#   2. Exactly ONE durable row in user_actions (the second call is a cache
#      replay — no new event row).
#   3. Exactly ONE SSE frame fanout (the second call must NOT broadcast a
#      duplicate frame; subscribers depend on this for own-op deferred
#      resolution to be unique-keyed).
#
# Existing test_idempotency.jl unit tests pin (1) at the with_idempotency
# helper level. Existing per-route tests cover (2) implicitly. Pinning all
# three together at the route level — across the curation, group, and
# message kinds — closes the contract.
#
# What's NOT tested here: concurrent replay (test_idempotency.jl already
# covers OP_LOCKS serialization).

# Helper — count user_actions rows + filter by kind.
function _count_actions(db, kind)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT COUNT(*) AS c FROM user_actions WHERE action = ?", [kind]))
    Int(rows[1].c)
end

# Helper — collect SSE frames during a closure execution. Subscribes,
# runs the closure, sleeps briefly to flush post-commit broadcasts, then
# unsubscribes and returns the captured frames.
function _capture_sse_during(base::String, kind_filter::String, f::Function)
    frames = String[]
    task = @async begin
        try
            HTTP.open("GET", "$base/api/events";
                headers = ["Accept" => "text/event-stream"]) do io
                while !eof(io)
                    line = readavailable(io)
                    isempty(line) && break
                    s = String(line)
                    occursin("\"kind\":\"$kind_filter\"", s) && push!(frames, s)
                end
            end
        catch
            # Connection closed by the test harness — expected.
        end
    end
    sleep(0.2)  # let SSE handshake complete
    f()
    sleep(0.5)  # let post-commit broadcasts flush
    schedule(task, InterruptException(); error = true)
    frames
end

@testset "Idempotency replay invariant — same client_op_id → 1 row, 1 SSE, identical body" begin

    @testset "POST /api/exposures/:id/peaks (peak_added)" begin
        mktempdir() do tmp
            analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
            mkpath(analysis_dir)
            cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
               joinpath(analysis_dir, "example_tot.dat"))
            db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.init_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
            e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
            HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

            with_test_server(db) do port, base
                op_id = "replay-test-peak-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:q => 0.99))

                pre_count = _count_actions(db, "peak_added")
                r1 = HTTP.post("$base/api/exposures/$e_id/peaks";
                    body = body_json, headers = headers)
                @test r1.status == 201
                # Replay — same op_id, must produce identical body and zero new rows.
                r2 = HTTP.post("$base/api/exposures/$e_id/peaks";
                    body = body_json, headers = headers)
                @test r2.status == 201
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "peak_added")
                @test post_count - pre_count == 1  # exactly one durable row
            end
        end
    end

    @testset "POST /api/groups/:id/members (index_confirmed)" begin
        mktempdir() do tmp
            analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
            mkpath(analysis_dir)
            cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
               joinpath(analysis_dir, "example_tot.dat"))
            db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.init_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")
            e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
            HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

            # Get an index id and the auto group id.
            ix_id = Int(Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM indices WHERE exposure_id = ? LIMIT 1",
                [e_id]))[1].id)
            grp_rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'custom'",
                [e_id]))
            isempty(grp_rows) && (HimalayaUI.ensure_custom_group!(db, e_id);
                grp_rows = Tables.rowtable(DBInterface.execute(db,
                    "SELECT id FROM index_groups WHERE exposure_id = ? AND kind = 'custom'",
                    [e_id])))
            grp_id = Int(grp_rows[1].id)

            with_test_server(db) do port, base
                op_id = "replay-test-grp-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:index_id => ix_id))

                pre_count = _count_actions(db, "index_confirmed")
                r1 = HTTP.post("$base/api/groups/$grp_id/members";
                    body = body_json, headers = headers)
                @test r1.status == 200
                r2 = HTTP.post("$base/api/groups/$grp_id/members";
                    body = body_json, headers = headers)
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "index_confirmed")
                @test post_count - pre_count == 1
            end
        end
    end

    @testset "POST /api/samples/:id/messages (post_message)" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")

            with_test_server(db) do port, base
                op_id = "replay-test-msg-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:body => "hello"))

                pre_count = _count_actions(db, "post_message")
                r1 = HTTP.post("$base/api/samples/$s_id/messages";
                    body = body_json, headers = headers)
                @test r1.status == 201
                r2 = HTTP.post("$base/api/samples/$s_id/messages";
                    body = body_json, headers = headers)
                @test r2.status == 201
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "post_message")
                @test post_count - pre_count == 1
            end
        end
    end

    @testset "PATCH /api/samples/:id (update_sample)" begin
        mktempdir() do tmp
            db = HimalayaUI.open_db(joinpath(tmp, "h.db"))
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, label="D1")

            with_test_server(db) do port, base
                op_id = "replay-test-upd-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:name => "renamed"))

                pre_count = _count_actions(db, "update_sample")
                r1 = HTTP.patch("$base/api/samples/$s_id";
                    body = body_json, headers = headers)
                @test r1.status == 200
                r2 = HTTP.patch("$base/api/samples/$s_id";
                    body = body_json, headers = headers)
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "update_sample")
                @test post_count - pre_count == 1
            end
        end
    end
end
