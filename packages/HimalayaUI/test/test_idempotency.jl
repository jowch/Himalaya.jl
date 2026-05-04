using Test, HTTP, SQLite, DBInterface, Tables
using HimalayaUI: with_idempotency, open_db

@testset "with_idempotency: passthrough when X-Client-Op-Id absent" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", Pair{String,String}[], UInt8[])
        called = Ref(0)
        result = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(200; body = "{\"x\":1}")
        end
        @test called[] == 1
        @test result.status == 200
        # Cache table should be empty.
        rows = Tables.rowtable(DBInterface.execute(db, "SELECT * FROM idempotent_responses"))
        @test isempty(rows)
    end
end

@testset "with_idempotency: cache hit returns prior response" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-1"], UInt8[])
        called = Ref(0)
        # First call.
        r1 = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(201; body = "{\"id\":42}")
        end
        # Second call with same op_id — body should NOT execute.
        r2 = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(500; body = "{\"x\":\"should not appear\"}")
        end
        @test called[] == 1
        @test r1.status == 201
        @test r2.status == 201
        @test String(r2.body) == "{\"id\":42}"
    end
end

@testset "with_idempotency: failures NOT cached" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-2"], UInt8[])
        called = Ref(0)
        # First call returns 4xx.
        r1 = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(400; body = "{\"error\":\"bad\"}")
        end
        # Second call should re-execute body (failures aren't cached).
        r2 = with_idempotency(db, req) do
            called[] += 1
            HTTP.Response(200; body = "{\"id\":99}")
        end
        @test called[] == 2
        @test r2.status == 200
        @test String(r2.body) == "{\"id\":99}"
        # Cache now contains the successful retry.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT * FROM idempotent_responses WHERE client_op_id = 'uuid-2'"))
        @test length(rows) == 1
    end
end

@testset "with_idempotency: concurrent retries serialize via OP_LOCKS" begin
    # Two parallel tasks with same op_id, both miss the cache initially.
    # Per-op-id ReentrantLock serializes; one writes, other reads cache.
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-3"], UInt8[])
        called = Ref(0)
        body_lock = ReentrantLock()
        function delay_body()
            lock(body_lock) do
                called[] += 1
            end
            sleep(0.05)
            HTTP.Response(201; body = "{\"called\":$(called[])}")
        end
        # Spawn two tasks racing.
        t1 = @async with_idempotency(db, req) do; delay_body(); end
        t2 = @async with_idempotency(db, req) do; delay_body(); end
        r1, r2 = fetch(t1), fetch(t2)
        # Body executed exactly once.
        @test called[] == 1
        # Both responses have the same body.
        @test String(r1.body) == String(r2.body)
    end
end

@testset "with_idempotency: response body is preserved for caller" begin
    # Regression: in Julia 1.x, `String(::Vector{UInt8})` takes ownership of
    # the buffer (empties the vector). Real route bodies (e.g. from
    # `Oxygen.json(...)`) are `Vector{UInt8}`, so caching must `copy` first
    # or the caller's response object is left with an empty body.
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-preserve"], UInt8[])
        r = with_idempotency(db, req) do
            HTTP.Response(200; body = Vector{UInt8}("{\"id\":1}"))
        end
        # The caller must still be able to read the body after the wrapper returns.
        @test String(r.body) == "{\"id\":1}"
        # The response from the cache lookup also has a readable body.
        r2 = with_idempotency(db, req) do
            HTTP.Response(500; body = Vector{UInt8}("{\"x\":\"should not appear\"}"))
        end
        @test String(r2.body) == "{\"id\":1}"
    end
end

@testset "gc_idempotent_responses! sweeps rows older than TTL" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        DBInterface.execute(db, """
            INSERT INTO idempotent_responses (client_op_id, status_code, body, created_at)
            VALUES ('old-1', 200, '{}', datetime('now', '-2 hours'))
        """)
        DBInterface.execute(db, """
            INSERT INTO idempotent_responses (client_op_id, status_code, body, created_at)
            VALUES ('new-1', 200, '{}', datetime('now', '-5 minutes'))
        """)
        HimalayaUI.gc_idempotent_responses!(db; ttl_seconds = 3600)
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM idempotent_responses ORDER BY client_op_id"))
        @test [String(r.client_op_id) for r in rows] == ["new-1"]
        SQLite.close(db)
    end
end

@testset "gc_idempotent_responses! also sweeps OP_LOCKS entries with no live response" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        # Trigger creation of two OP_LOCKS entries.
        req1 = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "live-1"], UInt8[])
        req2 = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "stale-1"], UInt8[])

        with_idempotency(db, req1) do
            HTTP.Response(200; body = Vector{UInt8}("{}"))
        end
        with_idempotency(db, req2) do
            HTTP.Response(200; body = Vector{UInt8}("{}"))
        end

        @test haskey(HimalayaUI.OP_LOCKS, "live-1")
        @test haskey(HimalayaUI.OP_LOCKS, "stale-1")

        # Manually expire the "stale-1" row so the GC sweep removes it.
        DBInterface.execute(db, """
            UPDATE idempotent_responses SET created_at = datetime('now', '-2 hours')
            WHERE client_op_id = 'stale-1'
        """)
        HimalayaUI.gc_idempotent_responses!(db; ttl_seconds = 3600)

        @test haskey(HimalayaUI.OP_LOCKS, "live-1")
        @test !haskey(HimalayaUI.OP_LOCKS, "stale-1")
        # Cleanup: don't leak OP_LOCKS state across tests.
        delete!(HimalayaUI.OP_LOCKS, "live-1")
        SQLite.close(db)
    end
end

@testset "with_idempotency: multi-event route cache hit returns identical body" begin
    # Simulate a route that emits N events (speculative POST shape).
    # First call commits 3 events with same client_op_id; second call returns cached.
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        # Seed FK chain so apply_event! on entity_id=1 succeeds.
        DBInterface.execute(db,
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, label) VALUES (1, 'A1')")
        DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")

        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "uuid-4"], UInt8[])
        called = Ref(0)
        function multi_event_body()
            called[] += 1
            # Emit 3 apply_event! calls with the same client_op_id (extracted from req).
            # Use distinct (action, entity_id) tuples — the I2 partial unique
            # index on (client_op_id, action, entity_id) rejects exact duplicates.
            for i in 1:3
                HimalayaUI.apply_event!(db, req;
                    kind = "speculative_created_$i",
                    entity_type = "exposure",
                    entity_id = 1,  # assume seeded
                    payload = Dict(:n => i))
            end
            HTTP.Response(201; body = "{\"created\":3}")
        end
        r1 = with_idempotency(db, req) do; multi_event_body(); end
        r2 = with_idempotency(db, req) do; multi_event_body(); end
        @test called[] == 1
        # 3 user_actions rows from first call only.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE client_op_id = 'uuid-4'"))
        @test length(rows) == 3
        # Cached body matches.
        @test String(r1.body) == String(r2.body)
    end
end

@testset "I2: with_idempotency body throw rolls back event AND cache atomically" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        DBInterface.execute(db,
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, label) VALUES (1, 'A1')")
        res = DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")
        exp_id = Int(DBInterface.lastrowid(res))

        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "op-atomic"], UInt8[])
        try
            with_idempotency(db, req) do
                HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
                    kind="peak_added", entity_type="exposure", entity_id=exp_id,
                    payload=Dict(:q => 4.0))
                error("body explosion")
            end
        catch
        end
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE client_op_id = 'op-atomic'"))
        @test isempty(rows)
        crows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM idempotent_responses WHERE client_op_id = 'op-atomic'"))
        @test isempty(crows)
        SQLite.close(db)
    end
end

@testset "OP_LOCKS entry removed when body throws (issue #37 Bug 4)" begin
    # Without cleanup, throw-without-cache-row leaks the OP_LOCKS entry until
    # process restart — gc_idempotent_responses! only collects locks for ops
    # that have a corresponding idempotent_responses row, and the row is
    # never written when the body throws.
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        req = HTTP.Request("POST", "/",
                           ["X-Client-Op-Id" => "op-throw-leak"], UInt8[])
        try
            with_idempotency(db, req) do
                error("body explosion")
            end
        catch
        end
        @test !haskey(HimalayaUI.OP_LOCKS, "op-throw-leak")
        # No idempotent_responses row was written either, so the GC would never
        # have collected this entry on its own.
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT 1 FROM idempotent_responses WHERE client_op_id = 'op-throw-leak'"))
        @test isempty(rows)
        SQLite.close(db)
    end
end

@testset "I2: with_idempotency success commits event AND cache atomically" begin
    mktempdir() do tmp
        db = open_db(joinpath(tmp, "test.db"))
        DBInterface.execute(db,
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, label) VALUES (1, 'A1')")
        res = DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")
        exp_id = Int(DBInterface.lastrowid(res))

        req = HTTP.Request("POST", "/", ["X-Client-Op-Id" => "op-success"], UInt8[])
        with_idempotency(db, req) do
            HimalayaUI.apply_event!(HimalayaUI.InTransaction(), db, req;
                kind="peak_added", entity_type="exposure", entity_id=exp_id,
                payload=Dict(:q => 5.0))
            HTTP.Response(200; body = Vector{UInt8}("{\"ok\":true}"))
        end
        @test length(Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM user_actions WHERE client_op_id = 'op-success'"))) == 1
        @test length(Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM idempotent_responses WHERE client_op_id = 'op-success'"))) == 1
        SQLite.close(db)
    end
end
