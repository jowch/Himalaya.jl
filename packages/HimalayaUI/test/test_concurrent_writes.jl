# Issue #122 — singleton _DB_REF write-path races.
#
# Pre-fix, two concurrent writers on the singleton could:
#   - Race on `SQLite.transaction(db)`'s TOCTOU and silently nest a SAVEPOINT
#     inside the other's BEGIN (Failure mode B: silent corruption).
#   - Race on the unlocked `db.stmt_wrappers` Dict mutation in `Stmt(db, sql)`.
#
# The fix adds a global `_DB_WRITE_LOCK::ReentrantLock` that serializes every
# `SQLite.transaction(db)` call against the singleton at the Julia level. This
# closes Race 1 entirely and Race 2 for writer-vs-writer (reader-vs-writer
# remains a tracked follow-up).
#
# These tests assert the *contract* of the lock (deterministic, channel-coordinated)
# rather than trying to trigger the underlying race (timing-dependent).

using Test, HTTP, SQLite, DBInterface, Tables
using HimalayaUI: with_idempotency, open_db
import HimalayaUI

@testset "issue #122: distinct-op-id writers serialize via _DB_WRITE_LOCK" begin
    # Two writers with DIFFERENT client_op_ids → different OP_LOCKS entries → the
    # only thing that can serialize them at the Julia level is _DB_WRITE_LOCK.
    # On unfixed code, t2 would silently nest a SAVEPOINT inside t1's BEGIN
    # (Failure mode B from #122) and complete before t1 releases.
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        ready     = Channel{Nothing}(1)
        release   = Channel{Nothing}(1)
        t2_armed  = Channel{Nothing}(1)

        t1 = @async with_idempotency(db,
            HTTP.Request("POST", "/", ["X-Client-Op-Id" => "lock-test-1"], UInt8[])) do
            put!(ready, nothing)
            take!(release)
            HTTP.Response(201; body = "{\"i\":1}")
        end

        take!(ready)  # t1 is inside its transaction, holding the write lock.

        # Different op-id, so per-op-id OP_LOCKS does NOT serialize this against t1.
        # `t2_armed` is published *just before* the with_idempotency call, so a
        # successful `take!(t2_armed)` proves t2 was scheduled and about to enter
        # the wrapper — without it the `!istaskdone(t2)` check would pass trivially
        # on a CI runner where t2 hadn't even started yet.
        t2 = @async begin
            put!(t2_armed, nothing)
            with_idempotency(db,
                HTTP.Request("POST", "/", ["X-Client-Op-Id" => "lock-test-2"], UInt8[])) do
                HTTP.Response(201; body = "{\"i\":2}")
            end
        end

        take!(t2_armed)
        # Small grace window for t2 to traverse the fast-path cache check + op-lock
        # acquire and reach the `_DB_WRITE_LOCK` attempt. The 50 ms is well above
        # the local-cache-check round-trip; if it's not enough on a slow CI host
        # the stress test still catches the BEGIN-vs-BEGIN regression.
        sleep(0.05)
        @test !istaskdone(t2)  # t2 must be waiting on _DB_WRITE_LOCK.

        put!(release, nothing)
        wait(t1); wait(t2)
        @test istaskdone(t1) && istaskdone(t2)

        # Both writes durably committed, no rollback, no orphans.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM idempotent_responses ORDER BY client_op_id"))
        @test [String(r.client_op_id) for r in rows] == ["lock-test-1", "lock-test-2"]
    end
end

@testset "issue #122: stress — N concurrent distinct-op writers all commit cleanly" begin
    # Spray N concurrent writers at the singleton with distinct op_ids. With the
    # write lock, they serialize; without it, two-or-more would race in
    # SQLite.transaction (nested savepoints, lost updates, or
    # "cannot start a transaction within a transaction" 500s).
    #
    # Each task records any exception it caught; we assert the error channel
    # is empty AND that the cache count matches N.
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        N = 25
        errors = Channel{Any}(N + 1)
        @sync for i in 1:N
            Threads.@spawn begin
                try
                    with_idempotency(db,
                        HTTP.Request("POST", "/",
                                     ["X-Client-Op-Id" => "stress-$i"],
                                     UInt8[])) do
                        # Force a yield-prone DB call inside the tx to widen the
                        # race window against any peer task on this connection.
                        DBInterface.execute(db, "SELECT 1")
                        HTTP.Response(201; body = "{\"i\":$i}")
                    end
                catch e
                    put!(errors, e)
                end
            end
        end
        close(errors)
        observed = collect(errors)
        @test isempty(observed)
        if !isempty(observed)
            @info "concurrent-writer stress saw $(length(observed)) errors" first=observed[1]
        end
        n_rows = first(Tables.rowtable(DBInterface.execute(db,
            "SELECT COUNT(*) AS n FROM idempotent_responses"))).n
        @test n_rows == N
    end
end

@testset "issue #122: write lock is reentrant within a single task" begin
    # Sanity: nested with_idempotency (or anything that re-enters the lock from
    # the SAME task) must not deadlock. The lock is ReentrantLock, not Mutex.
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        outer_called = Ref(0)
        inner_called = Ref(0)
        with_idempotency(db,
            HTTP.Request("POST", "/", Pair{String,String}[], UInt8[])) do
            outer_called[] += 1
            # Re-enter from the same task — must succeed (ReentrantLock).
            with_idempotency(db,
                HTTP.Request("POST", "/", Pair{String,String}[], UInt8[])) do
                inner_called[] += 1
                HTTP.Response(200; body = "{}")
            end
        end
        @test outer_called[] == 1
        @test inner_called[] == 1
    end
end
