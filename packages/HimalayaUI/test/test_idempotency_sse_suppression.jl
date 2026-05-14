using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI: with_idempotency, open_db, InTransaction, apply_event!

# SSE suppression on rollback / 4xx (issue #125 items 2 & 4).
#
# `with_idempotency` and bare `SQLite.transaction` both rely on
# `_clear_post_commit_broadcasts!` to drop queued SSE frames before exit when
# the tx body either throws or returns 4xx. test_idempotency_replay_invariant.jl
# covers the success path; these testsets cover the suppression paths:
#
#   1. `with_idempotency` body throw — three sites in idempotency.jl
#      (no-op-id catch at line 98, op-id catch at line 139, replayed_cache /
#      ≥400 branch at line 161). All three must produce zero SSE frames
#      even though `apply_event!(InTransaction, ...)` (or a bare
#      `_enqueue_post_commit_broadcast!`) enqueued one.
#   2. Bare `SQLite.transaction` body throw — pins the ordering invariant
#      that post-commit broadcasts only fire AFTER commit. A throw before
#      commit must produce zero delivered frames even though a deferred
#      broadcast was queued, and the rolled-back tx must leave no
#      user_actions row.

# Inline copy of test_idempotency_replay_invariant.jl's helper so this file
# is self-contained (avoids include-order coupling).
function _capture_sse(f::Function, kind_filter::String)
    pending = Channel{String}(64)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        f()
        # The post-commit queue fires synchronously on the same task, so no
        # sleep is needed for the rollback paths (queue is cleared, never
        # drained). Keep a small slack to absorb any future async refactor.
        sleep(0.1)
    finally
        lock(HimalayaUI.SSE_LOCK) do
            filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
        end
        close(pending)
    end
    frames = String[]
    for frame in pending
        startswith(frame, ":") && continue
        occursin("\"kind\":\"$kind_filter\"", frame) && push!(frames, frame)
    end
    frames
end

# Insert the minimal experiment → sample → exposure rows the `peak_added`
# dispatcher's view-update needs to satisfy FK constraints.
function _seed_exposure(db::SQLite.DB)
    DBInterface.execute(db,
        "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
    DBInterface.execute(db,
        "INSERT INTO samples (experiment_id, name) VALUES (1, 'A1')")
    res = DBInterface.execute(db,
        "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")
    Int(DBInterface.lastrowid(res))
end

@testset "with_idempotency: SSE queue cleared on body throw / 4xx (issue #125 item 2)" begin

    @testset "no-op-id path: throw clears queue, zero frames delivered" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "test.db"))
            exp_id = _seed_exposure(db)
            # Header set deliberately omits X-Client-Op-Id → with_idempotency
            # takes the no-op-id branch (line 91 in idempotency.jl).
            req = HTTP.Request("POST", "/",
                ["X-Client-Id" => "tab-A", "X-Username" => "alice"], UInt8[])

            frames = _capture_sse("peak_added") do
                try
                    with_idempotency(db, req) do
                        result = apply_event!(InTransaction(), db, req;
                            kind = "peak_added", entity_type = "exposure",
                            entity_id = exp_id, payload = Dict(:q => 0.5))
                        HimalayaUI._enqueue_broadcast_from_result!(result,
                            "peak_added", "exposure", exp_id)
                        error("body explosion: no-op-id branch")
                    end
                catch
                end
            end
            @test isempty(frames)
            # And the rolled-back tx leaves no durable event row.
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT 1 FROM user_actions WHERE action='peak_added' AND entity_id=?",
                [exp_id]))
            @test isempty(rows)
            SQLite.close(db)
        end
    end

    @testset "op-id path: throw clears queue, zero frames delivered" begin
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "test.db"))
            exp_id = _seed_exposure(db)
            req = HTTP.Request("POST", "/",
                ["X-Client-Id" => "tab-A", "X-Username" => "alice",
                 "X-Client-Op-Id" => "op-sse-throw-1"], UInt8[])

            frames = _capture_sse("peak_added") do
                try
                    with_idempotency(db, req) do
                        result = apply_event!(InTransaction(), db, req;
                            kind = "peak_added", entity_type = "exposure",
                            entity_id = exp_id, payload = Dict(:q => 0.6))
                        HimalayaUI._enqueue_broadcast_from_result!(result,
                            "peak_added", "exposure", exp_id)
                        error("body explosion: op-id branch")
                    end
                catch
                end
            end
            @test isempty(frames)
            # OP_LOCKS entry must be cleaned up too (issue #37 Bug 4 already
            # covered in test_idempotency.jl; re-check here so this test
            # doesn't leak state into subsequent tests).
            @test !haskey(HimalayaUI.OP_LOCKS, "op-sse-throw-1")
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT 1 FROM user_actions WHERE client_op_id='op-sse-throw-1'"))
            @test isempty(rows)
            SQLite.close(db)
        end
    end

    @testset "post-tx 4xx branch: queue cleared even though tx committed" begin
        # Body returns 400 without throwing → tx commits → line 161 in
        # idempotency.jl clears the queue rather than flushing it. The
        # docstring on _clear_post_commit_broadcasts! captures the intent:
        # "Failed-but-not-thrown bodies (status≥400) … any speculative
        # enqueues are dropped." This is an artificial route shape (real
        # 4xx routes validate before any apply_event!), but the suppression
        # mechanism must hold even if a route mis-orders its writes.
        mktempdir() do tmp
            db = open_db(joinpath(tmp, "test.db"))
            exp_id = _seed_exposure(db)
            req = HTTP.Request("POST", "/",
                ["X-Client-Id" => "tab-A",
                 "X-Client-Op-Id" => "op-sse-4xx-1"], UInt8[])

            frames = _capture_sse("peak_added") do
                with_idempotency(db, req) do
                    # Speculative enqueue without a real apply_event! — the
                    # branch we're testing fires before the queue contents
                    # are inspected.
                    HimalayaUI._enqueue_post_commit_broadcast!(
                        999, "peak_added", "exposure", exp_id,
                        nothing, "tab-A", "op-sse-4xx-1", "{\"q\":0.7}")
                    HTTP.Response(400; body = Vector{UInt8}("{\"error\":\"bad\"}"))
                end
            end
            @test isempty(frames)
            # Cache row NOT written (4xx is not cached).
            crows = Tables.rowtable(DBInterface.execute(db,
                "SELECT 1 FROM idempotent_responses WHERE client_op_id='op-sse-4xx-1'"))
            @test isempty(crows)
            # Clean up the OP_LOCKS entry the route helper just registered.
            delete!(HimalayaUI.OP_LOCKS, "op-sse-4xx-1")
            SQLite.close(db)
        end
    end
end

@testset "Deferred broadcast: SQLite.transaction throw before commit → zero SSE (issue #125 item 4)" begin
    # Pins the ordering invariant: a post-commit broadcast enqueued inside a
    # `SQLite.transaction` must NOT fire if the tx body raises before commit,
    # even though the InTransaction variant's INSERT was reached and queued
    # a broadcast via _enqueue_broadcast_from_result!. The catch site is
    # responsible for calling _clear_post_commit_broadcasts! (in production,
    # `with_idempotency` does this; tests must mirror it explicitly).
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        exp_id = _seed_exposure(db)
        req = HTTP.Request("POST", "/",
            ["X-Client-Id" => "tab-A", "X-Username" => "alice"], UInt8[])

        frames = _capture_sse("peak_added") do
            try
                SQLite.transaction(db) do
                    result = apply_event!(InTransaction(), db, req;
                        kind = "peak_added", entity_type = "exposure",
                        entity_id = exp_id, payload = Dict(:q => 0.8))
                    HimalayaUI._enqueue_broadcast_from_result!(result,
                        "peak_added", "exposure", exp_id)
                    error("tx body explosion before commit")
                end
            catch
                # In production, with_idempotency or another wrapper would
                # call this; mirror it so the queue doesn't leak into the
                # surrounding test task and accidentally flush on a later
                # success path.
                HimalayaUI._clear_post_commit_broadcasts!()
            end
        end
        @test isempty(frames)
        # The rolled-back tx leaves no event row — sanity check that the
        # transaction boundary actually rolled the InTransaction INSERT back.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM user_actions WHERE action='peak_added' AND entity_id=?",
            [exp_id]))
        @test isempty(rows)
        SQLite.close(db)
    end
end
