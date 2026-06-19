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

# Helper — register a direct in-process SSE subscriber (no HTTP), run the
# closure, then drain captured frames. Subscribers in HimalayaUI.SSE_SUBSCRIBERS
# receive a JSON-encoded `event: curation\ndata: {...}\n\n` string per
# `broadcast_event!` call. Filtering by `kind_filter` (the event kind from
# the JSON) lets us count only the frames the test cares about.
#
# Why this matters: the third invariant in this file's docstring is "exactly
# one SSE fanout per replayed op." If `with_idempotency`'s cache-replay path
# ever stops short-circuiting the broadcast (e.g. a refactor removes the
# `result.cached` skip), this helper is the only thing that catches it —
# byte-equal HTTP body and durable-row count would both still pass.
function _capture_sse_during(f::Function, kind_filter::String)
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        f()
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(pending)
    end
    frames = String[]
    for frame in pending
        # Skip heartbeats (": heartbeat" lines) and frames of other kinds.
        startswith(frame, ":") && continue
        occursin("\"kind\":\"$kind_filter\"", frame) && push!(frames, frame)
    end
    frames
end

@testset "Idempotency replay invariant — same client_op_id → 1 row, 1 SSE, identical body" begin

    @testset "POST /api/exposures/:id/peaks (peak_added)" begin
        mktempdir() do tmp
            analysis_dir = joinpath(tmp, "analysis", "automatic_analysis")
            mkpath(analysis_dir)
            cp(joinpath(@__DIR__, "..", "..", "..", "test", "data", "example_tot.dat"),
               joinpath(analysis_dir, "example_tot.dat"))
            db = open_prepared_clone(tmp)
            exp_id = HimalayaUI.init_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=analysis_dir)
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
            e_id = HimalayaUI.create_exposure!(db; sample_id=s_id, filename="example_tot")
            HimalayaUI.analyze_exposure!(db, e_id, analysis_dir)

            with_inproc_routes(db) do call
                op_id = "replay-test-peak-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:q => 0.99))

                pre_count = _count_actions(db, "peak_added")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("peak_added") do
                    r1 = call("POST", "/api/exposures/$e_id/peaks";
                        headers = headers, body = Vector{UInt8}(body_json))
                    # Replay — same op_id, must produce identical body and zero new rows.
                    r2 = call("POST", "/api/exposures/$e_id/peaks";
                        headers = headers, body = Vector{UInt8}(body_json))
                end
                @test r1.status == 201
                @test r2.status == 201
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "peak_added")
                @test post_count - pre_count == 1  # exactly one durable row
                @test length(frames) == 1          # exactly one SSE fanout
            end
        end
    end

    # D-10: the POST /api/groups/:id/members (index_confirmed) idempotency-replay
    # invariant was retired with the route. Its assignment-native successor — the
    # assignment_add member route — carries the same invariant, directly below.

    @testset "POST /api/exposures/:id/assignment/members (assignment_add)" begin
        # Plan D Task D-3 layer-1/2: the native assignment member route must be
        # idempotent under retry — one durable row, one SSE frame, identical body.
        mktempdir() do tmp
            db = open_prepared_clone(tmp)
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")
            e_id = HimalayaUI.create_exposure!(db; sample_id=s_id)
            DBInterface.execute(db,
                "INSERT INTO indices (id, exposure_id, phase, basis) VALUES (77, ?, 'Pn3m', 0.1)", [e_id])

            with_inproc_routes(db) do call
                op_id = "replay-test-asg-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:index_id => 77))

                pre_count = _count_actions(db, "assignment_add")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("assignment_add") do
                    r1 = call("POST", "/api/exposures/$e_id/assignment/members";
                        headers = headers, body = Vector{UInt8}(body_json))
                    r2 = call("POST", "/api/exposures/$e_id/assignment/members";
                        headers = headers, body = Vector{UInt8}(body_json))
                end
                @test r1.status == 200
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                @test _count_actions(db, "assignment_add") - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    @testset "POST /api/samples/:id/messages (post_message)" begin
        mktempdir() do tmp
            db = open_prepared_clone(tmp)
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")

            with_inproc_routes(db) do call
                op_id = "replay-test-msg-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:body => "hello"))

                pre_count = _count_actions(db, "post_message")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("post_message") do
                    r1 = call("POST", "/api/samples/$s_id/messages";
                        headers = headers, body = Vector{UInt8}(body_json))
                    r2 = call("POST", "/api/samples/$s_id/messages";
                        headers = headers, body = Vector{UInt8}(body_json))
                end
                @test r1.status == 201
                @test r2.status == 201
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "post_message")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    @testset "PATCH /api/samples/:id (update_sample)" begin
        mktempdir() do tmp
            db = open_prepared_clone(tmp)
            exp_id = HimalayaUI.create_experiment!(db; path=tmp,
                data_dir=joinpath(tmp,"data"), analysis_dir=joinpath(tmp,"analysis"))
            s_id = HimalayaUI.create_sample!(db; experiment_id=exp_id, name="D1")

            with_inproc_routes(db) do call
                op_id = "replay-test-upd-$(rand(UInt32))"
                headers = ["Content-Type" => "application/json",
                           "X-Username"   => "alice",
                           "X-Client-Id"  => "tab-1",
                           "X-Client-Op-Id" => op_id]
                body_json = JSON3.write(Dict(:display_name => "Renamed Display"))

                pre_count = _count_actions(db, "update_sample")
                r1 = nothing; r2 = nothing
                frames = _capture_sse_during("update_sample") do
                    r1 = call("PATCH", "/api/samples/$s_id";
                        headers = headers, body = Vector{UInt8}(body_json))
                    r2 = call("PATCH", "/api/samples/$s_id";
                        headers = headers, body = Vector{UInt8}(body_json))
                end
                @test r1.status == 200
                @test r2.status == 200
                @test String(r2.body) == String(r1.body)
                post_count = _count_actions(db, "update_sample")
                @test post_count - pre_count == 1
                @test length(frames) == 1
            end
        end
    end

    # I3.6 (#177): the comparison replay-invariant testsets (comparison_created
    # / _submitted / _deleted / comparison_message / pinned / unpinned + the
    # 409-retry) drove the invariant through the deleted /api/comparisons*
    # routes. Those routes are retired with the Compare page; the invariant
    # itself stays covered by the surviving kinds above (peak_added,
    # assignment_add, post_message, update_sample). The kept comparison_*
    # dispatcher branches' fold is covered by test_events.jl.
    @testset "comparison replay invariants (retired with the Compare routes, #177)" begin
        @test true  # placeholder; routes deleted in I3.6
    end
end

# I3.6 (#177): the "SSE frame includes view_* on comparison_submitted (A-6)"
# testset drove the deleted POST /api/comparisons + /submit routes to assert
# the SSE post_state envelope. The routes are retired with the Compare page;
# the kept `comparison_submitted` dispatcher branch's fold is covered by
# test_events.jl, and the general SSE-fanout-per-op invariant by the surviving
# kinds in the replay-invariant testset above.
