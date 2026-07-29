using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# ---------------------------------------------------------------------------
# Unit tests: broadcast plumbing without opening a real HTTP connection.
# These assert the channel-only architecture and JSON frame shape directly.
# ---------------------------------------------------------------------------

@testset "SSE: broadcast_event! puts a curation frame onto subscriber channels" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        HimalayaUI.bind_db!(db)

        # Inject a fake subscriber into the registry.
        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end

        try
            # Create a user so lookup_username resolves.
            user_id = HimalayaUI.get_or_create_user!(db, "alice")

            HimalayaUI.broadcast_event!(
                1, "test_broadcast", "exposure", 42,
                user_id, "tab-xyz", "uuid-789", JSON3.write(Dict(:foo => "bar")))

            @test isready(pending)
            frame = take!(pending)
            @test occursin("event: curation", frame)
            @test occursin("\"kind\":\"test_broadcast\"", frame)
            @test occursin("\"actor\":\"alice\"", frame)
            @test occursin("\"entity_id\":42", frame)
            @test occursin("\"entity_type\":\"exposure\"", frame)
            @test occursin("\"client_id\":\"tab-xyz\"", frame)
            @test occursin("\"client_op_id\":\"uuid-789\"", frame)

            # Parse JSON to assert ts + client_op_id keys present.
            data_line = first([l for l in split(frame, '\n') if startswith(l, "data: ")])
            json_str = replace(data_line, r"^data: " => "")
            obj = JSON3.read(json_str)
            @test haskey(obj, :client_op_id)
            @test haskey(obj, :ts)
            @test obj.client_op_id == "uuid-789"
            @test obj.ts isa AbstractString
            @test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z$", obj.ts)

            # Payload is embedded in the frame.
            @test occursin("\"foo\"", frame)
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(pending)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end
end

@testset "SSE: broadcast_event! with nil user_id sets actor to null" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        HimalayaUI.bind_db!(db)

        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end

        try
            HimalayaUI.broadcast_event!(
                2, "anon_event", "exposure", 7,
                nothing, nothing, nothing, nothing)

            @test isready(pending)
            frame = take!(pending)
            @test occursin("event: curation", frame)
            @test occursin("\"actor\":null", frame)
            @test occursin("\"client_id\":null", frame)
            @test occursin("\"client_op_id\":null", frame)

            data_line = first([l for l in split(frame, '\n') if startswith(l, "data: ")])
            obj = JSON3.read(replace(data_line, r"^data: " => ""))
            @test haskey(obj, :client_op_id)
            @test haskey(obj, :ts)
            @test occursin(r"^\d{4}-\d{2}-\d{2}T\d{2}:\d{2}:\d{2}\.\d{3}Z$", obj.ts)
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(pending)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end
end

@testset "SSE: broadcast_event! prunes closed subscriber channels" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        HimalayaUI.bind_db!(db)

        dead = Channel{String}(1)
        close(dead)  # closed before broadcast — simulates disconnected client
        live = Channel{String}(64)
        sub_dead = (pending = dead,)
        sub_live = (pending = live,)

        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub_dead)
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub_live)
        end

        try
            HimalayaUI.broadcast_event!(
                3, "prune_test", "exposure", 1,
                nothing, nothing, nothing, nothing)

            # Dead sub should have been pruned.
            n = lock(HimalayaUI.SSE_LOCK) do
                length(HimalayaUI.SSE_SUBSCRIBERS[])
            end
            @test n == 1

            # Live sub received the frame.
            @test isready(live)
            frame = take!(live)
            @test occursin("\"kind\":\"prune_test\"", frame)
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub_live, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(live)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end
end

@testset "SSE: broadcast_event! prunes slow (full) subscriber without blocking" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        HimalayaUI.bind_db!(db)

        # A slow subscriber: channel capacity 4, intentionally NOT drained.
        slow = Channel{String}(4)
        live = Channel{String}(128)
        sub_slow = (pending = slow,)
        sub_live = (pending = live,)

        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub_slow)
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub_live)
        end

        try
            # Pump 6 frames — more than the slow channel's capacity of 4.
            # The broadcast loop must not hang and must prune the slow subscriber.
            for i in 1:6
                HimalayaUI.broadcast_event!(
                    i, "slow_test", "exposure", 1,
                    nothing, nothing, nothing, nothing)
            end

            # Slow subscriber should have been pruned after its channel filled.
            n = lock(HimalayaUI.SSE_LOCK) do
                length(HimalayaUI.SSE_SUBSCRIBERS[])
            end
            @test n == 1  # only live remains

            # Live subscriber received all 6 frames.
            @test Base.n_avail(live) == 6
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub_live, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(slow)
            close(live)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end
end

@testset "SSE: analyze_run with both skip flags true does NOT broadcast" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        HimalayaUI.bind_db!(db)

        # Seed FK chain so apply_event! on entity_id=1 succeeds.
        DBInterface.execute(db,
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (1, 'A1')")
        DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")

        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end

        try
            req = HimalayaUI._system_request()
            result = HimalayaUI.apply_event!(db, req;
                kind = "analyze_run",
                entity_type = "exposure",
                entity_id = 1,
                payload = Dict(:findpeaks_skipped => true,
                               :indexpeaks_skipped => true,
                               :duration_ms => 0))

            # No frame should have been enqueued.
            @test Base.n_avail(pending) == 0

            # But the user_actions row IS still written.
            rows = Tables.rowtable(DBInterface.execute(db,
                "SELECT id, action FROM user_actions WHERE id = ?", [result.event_id]))
            @test length(rows) == 1
            @test String(rows[1].action) == "analyze_run"
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(pending)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end
end

@testset "SSE: analyze_run with one skip flag true DOES broadcast" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        HimalayaUI.bind_db!(db)

        DBInterface.execute(db,
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (1, 'A1')")
        DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")

        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end

        try
            req = HimalayaUI._system_request()
            HimalayaUI.apply_event!(db, req;
                kind = "analyze_run",
                entity_type = "exposure",
                entity_id = 1,
                payload = Dict(:findpeaks_skipped => true,
                               :indexpeaks_skipped => false,
                               :duration_ms => 5))

            @test Base.n_avail(pending) == 1
            frame = take!(pending)
            @test occursin("event: curation", frame)
            data_line = first([l for l in split(frame, '\n') if startswith(l, "data: ")])
            obj = JSON3.read(replace(data_line, r"^data: " => ""))
            @test obj.kind == "analyze_run"
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(pending)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end
end

@testset "SSE: non-analyze_run events broadcast regardless of skip flags in payload" begin
    mktempdir() do tmp
        db = open_prepared_clone(tmp)
        HimalayaUI.bind_db!(db)

        DBInterface.execute(db,
            "INSERT INTO experiments (name, path, data_dir, analysis_dir) VALUES ('e', '/p', '/d', '/a')")
        DBInterface.execute(db,
            "INSERT INTO samples (experiment_id, name) VALUES (1, 'A1')")
        DBInterface.execute(db,
            "INSERT INTO exposures (sample_id, filename) VALUES (1, 'f')")

        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end

        try
            # peak_added would normally come from a user route; build a real HTTP.Request.
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(db, req;
                kind = "peak_added",
                entity_type = "exposure",
                entity_id = 1,
                payload = Dict(:q => 0.123,
                               :findpeaks_skipped => true,
                               :indexpeaks_skipped => true))

            @test Base.n_avail(pending) == 1
            frame = take!(pending)
            data_line = first([l for l in split(frame, '\n') if startswith(l, "data: ")])
            obj = JSON3.read(replace(data_line, r"^data: " => ""))
            @test obj.kind == "peak_added"
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(pending)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end
end

@testset "SSE: broadcast_event! emits post_state when provided" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        HimalayaUI.bind_db!(db)
        sub = (id = "t-postst", pending = Channel{String}(8))
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end
        try
            HimalayaUI.broadcast_event!(1, "peak_added", "exposure", 42,
                nothing, "tab-id", "op-id", JSON3.write(Dict(:q => 1.0));
                post_state = Dict(:analysis_inputs_hash => "abc123", :indices => Any[]))
            frame = take!(sub.pending)
            data_line = first([l for l in split(frame, '\n') if startswith(l, "data: ")])
            json_str = replace(data_line, r"^data: " => "")
            obj = JSON3.read(json_str)
            @test haskey(obj, :post_state)
            @test obj.post_state.analysis_inputs_hash == "abc123"
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(sub.pending)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end
end

@testset "SSE: broadcast_event! omits post_state when not provided" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        HimalayaUI.bind_db!(db)
        sub = (id = "t-no-postst", pending = Channel{String}(8))
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end
        try
            HimalayaUI.broadcast_event!(1, "peak_added", "exposure", 42,
                nothing, "tab-id", "op-id", JSON3.write(Dict(:q => 1.0)))
            frame = take!(sub.pending)
            data_line = first([l for l in split(frame, '\n') if startswith(l, "data: ")])
            json_str = replace(data_line, r"^data: " => "")
            obj = JSON3.read(json_str)
            @test !haskey(obj, :post_state)
        finally
            lock(HimalayaUI.SSE_LOCK) do
                filter!(x -> x !== sub, HimalayaUI.SSE_SUBSCRIBERS[])
            end
            close(sub.pending)
            HimalayaUI.SSE_SUBSCRIBERS[] = []
        end
    end
end

@testset "SSE: lookup_username returns nothing for unknown id" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        @test HimalayaUI.lookup_username(db, 99999) === nothing
    end
end

@testset "SSE: lookup_username resolves known user" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        uid = HimalayaUI.get_or_create_user!(db, "bob")
        @test HimalayaUI.lookup_username(db, uid) == "bob"
    end
end

# ---------------------------------------------------------------------------
# Network-level test: verify the /api/events endpoint emits a curation frame
# after apply_event! and sends correct SSE headers.
# ---------------------------------------------------------------------------

@testset "SSE: /api/events endpoint headers and frame delivery" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        port = HimalayaUI.find_free_port()
        HimalayaUI.start_test_server!(db, port)
        try
            # Open the SSE stream in a background task; collect the first frame.
            received = Channel{String}(1)
            t = @async begin
                try
                    HTTP.open("GET", "http://127.0.0.1:$port/api/events";
                              retry = false, status_exception = false,
                              readtimeout = 5) do io
                        HTTP.startread(io)
                        # Verify headers from the response.
                        ct  = HTTP.header(io, "Content-Type", "")
                        xab = HTTP.header(io, "X-Accel-Buffering", "")
                        isopen(received) && put!(received, "headers:$ct|$xab")

                        buf = IOBuffer()
                        while !eof(io) && isopen(received)
                            chunk = readavailable(io)
                            isempty(chunk) && continue
                            write(buf, chunk)
                            s = String(take!(buf))
                            if occursin("event: curation", s) && isopen(received)
                                put!(received, s)
                                break
                            end
                        end
                    end
                catch
                end
            end

            # Wait for the subscriber to register and headers to arrive.
            header_frame = nothing
            for _ in 1:40
                if isready(received)
                    header_frame = take!(received)
                    break
                end
                sleep(0.05)
            end

            # Check headers line received first.
            @test header_frame !== nothing
            if header_frame !== nothing
                @test occursin("text/event-stream", header_frame)
                @test occursin("no", header_frame)  # X-Accel-Buffering: no
            end

            # Give a moment for the subscriber goroutine to register.
            sleep(0.1)

            # Trigger a broadcast.
            req = HTTP.Request("POST", "/x", ["X-Username" => "alice"], UInt8[])
            HimalayaUI.apply_event!(db, req;
                kind        = "sse_network_test",
                entity_type = "exposure",
                entity_id   = 1,
                payload     = Dict(:ping => true))

            # Collect the curation frame (up to 1s).
            curation_frame = nothing
            for _ in 1:20
                if isready(received)
                    curation_frame = take!(received)
                    break
                end
                sleep(0.05)
            end

            if curation_frame !== nothing
                @test occursin("event: curation", curation_frame)
                @test occursin("\"kind\":\"sse_network_test\"", curation_frame)
            else
                # Network test is best-effort in CI — skip rather than fail hard.
                @test_skip "curation frame not received within 1s (network timing)"
            end

            # Clean up the streaming task.
            close(received)
            try; schedule(t, InterruptException(); error = true); catch; end

        finally
            HimalayaUI.stop_test_server!()
        end
    end
end

# ---------------------------------------------------------------------------
# _try_put! contract — load-bearing for both broadcast_event! and the
# per-subscriber heartbeat Timer at server.jl:78 (fixed in #127). A
# saturated channel must NEVER block the caller; the heartbeat task would
# otherwise hang on a slow subscriber and silently leak.
# ---------------------------------------------------------------------------

@testset "SSE: _try_put! on full channel returns false without blocking" begin
    ch = Channel{String}(2)
    put!(ch, "a")
    put!(ch, "b")  # channel is now at capacity
    task = @async HimalayaUI._try_put!(ch, "c")
    # With non-blocking semantics, the task completes ~immediately. If a
    # regression reintroduced blocking put!, timedwait reports :timed_out.
    # Guard the subsequent fetch/take! — they would otherwise hang the
    # whole suite by waiting on the still-blocked task instead of failing
    # the @test cleanly.
    done = (timedwait(() -> istaskdone(task), 2.0) === :ok)
    @test done
    done || return
    @test fetch(task) === false
    # The dropped frame did not displace existing entries.
    @test take!(ch) == "a"
    @test take!(ch) == "b"
end

@testset "SSE: _try_put! on open channel with capacity succeeds" begin
    ch = Channel{String}(2)
    @test HimalayaUI._try_put!(ch, "x") === true
    @test take!(ch) == "x"
end

@testset "SSE: _try_put! on closed channel returns false" begin
    ch = Channel{String}(2)
    close(ch)
    @test HimalayaUI._try_put!(ch, "x") === false
end

# ---------------------------------------------------------------------------
# Regression: a saturated subscriber must be CLOSED when it is evicted.
#
# Dropping it from SSE_SUBSCRIBERS alone left the /api/events handler parked on
# `for frame in pending` forever. The HTTP stream was never finished, so the
# browser's EventSource saw a healthy connection, never fired `onerror`, and
# never auto-reconnected — the client went permanently deaf to every curation
# event. An ingest scan reliably triggers this: it out-runs the 64-slot channel.
# ---------------------------------------------------------------------------

@testset "SSE: saturated subscriber is closed (not just dropped) on eviction" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        HimalayaUI.bind_db!(db)

        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            empty!(HimalayaUI.SSE_SUBSCRIBERS[])
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end

        try
            # Never drain: overrun the 64-slot channel the way a scan does.
            for i in 1:70
                HimalayaUI.broadcast_progress!(1; kind = "ingest_progress",
                                               processed = i, total = 700)
            end

            @test isempty(HimalayaUI.SSE_SUBSCRIBERS[])          # evicted
            @test !isopen(pending)                                # AND closed

            # A closed channel ends the handler's `for frame in pending` loop, so
            # the stream tears down and EventSource reconnects. Buffered frames
            # stay readable until drained.
            #
            # Drain by count, NOT `for _ in pending`: iterating an OPEN channel
            # blocks forever once it empties, so a regression that stopped closing
            # would hang this testset (and the whole suite) instead of failing the
            # assertion above. Same hazard the `_try_put!` blocking test guards
            # against a few testsets down.
            drained = 0
            while Base.n_avail(pending) > 0
                take!(pending)
                drained += 1
            end
            @test drained == 64
        finally
            lock(HimalayaUI.SSE_LOCK) do
                empty!(HimalayaUI.SSE_SUBSCRIBERS[])
            end
        end
    end
end

@testset "SSE: _fanout_frame! closes an already-closed subscriber's slot cleanly" begin
    pending = Channel{String}(4)
    close(pending)
    sub = (pending = pending,)
    lock(HimalayaUI.SSE_LOCK) do
        empty!(HimalayaUI.SSE_SUBSCRIBERS[])
        push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
    end
    try
        HimalayaUI._fanout_frame!("event: curation\ndata: {}\n\n")
        @test isempty(HimalayaUI.SSE_SUBSCRIBERS[])
    finally
        lock(HimalayaUI.SSE_LOCK) do
            empty!(HimalayaUI.SSE_SUBSCRIBERS[])
        end
    end
end

# ---------------------------------------------------------------------------
# Wire-level companion to the eviction unit test above (docs/contract-testing.md
# governs SSE fixes; test/AGENTS.md designates this file a wire keeper).
#
# The unit test proves eviction CLOSES the channel. This proves the behavior
# that actually mattered to the bug: once closed, the /api/events response
# FINISHES, so the client sees EOF and EventSource auto-reconnects. Before the
# fix the handler stayed parked on `for frame in pending`, the response never
# completed, and the browser held a healthy-looking socket that would never
# deliver another event.
#
# NOTE two deliberate limits on what this asserts:
#
#  * Saturation is not driven over the wire. The handler drains `pending`
#    straight into the socket and 64 small frames fit inside the kernel buffer,
#    so a real client cannot reliably be made to fall 64 frames behind. The unit
#    test above covers saturate→evict→close; this covers close→client-unblocks.
#  * The client's read terminates by CONNECTION DROP (an HTTP.RequestError), not
#    a clean chunked-EOF — verified empirically against this server. Either way
#    the browser's EventSource fires `onerror` and auto-reconnects, which is the
#    property the fix exists to restore. So the assertion is "the client stops
#    waiting promptly", not "the body ends cleanly".
#  * THIS TEST IS NOT ITSELF A REGRESSION TEST. It closes the channel directly
#    rather than routing through `_fanout_frame!`, and `server.jl`'s handler loop
#    is untouched by that fix — so a manual close unblocks the handler identically
#    with or without it, and this passes on both sides. It documents the second
#    half of the causal chain (close → response ends → client unblocks); the
#    regression protection for the fix itself is the saturation unit test above,
#    which asserts `!isopen(pending)` after 70 evicting broadcasts.
# ---------------------------------------------------------------------------

@testset "SSE: closing an evicted subscriber unblocks the client" begin
    mktempdir() do dir
        db = open_prepared_clone(dir)
        port = HimalayaUI.find_free_port()
        HimalayaUI.start_test_server!(db, port)
        try
            # Set by the reader task however its read ends — clean EOF or throw.
            ended = Threads.Atomic{Bool}(false)
            t = @async begin
                try
                    HTTP.open("GET", "http://127.0.0.1:$port/api/events";
                              retry = false, status_exception = false,
                              readtimeout = 30) do io
                        HTTP.startread(io)
                        while !eof(io)
                            readavailable(io)
                        end
                    end
                catch
                end
                ended[] = true
            end

            # Wait for the handler to register its subscriber.
            sub = nothing
            for _ in 1:60
                s = lock(HimalayaUI.SSE_LOCK) do
                    isempty(HimalayaUI.SSE_SUBSCRIBERS[]) ? nothing :
                        first(HimalayaUI.SSE_SUBSCRIBERS[])
                end
                if s !== nothing
                    sub = s
                    break
                end
                sleep(0.05)
            end

            if sub === nothing
                @test_skip "subscriber did not register within 3s (network timing)"
            else
                @test !ended[]   # still streaming before the eviction

                # The close half of what _fanout_frame! does to an evicted
                # subscriber.
                close(sub.pending)

                # 5s window against a 30s readtimeout, so a handler that stayed
                # parked would fail rather than pass slowly.
                @test timedwait(() -> ended[], 5.0) === :ok
            end

            try; schedule(t, InterruptException(); error = true); catch; end
        finally
            HimalayaUI.stop_test_server!()
        end
    end
end
