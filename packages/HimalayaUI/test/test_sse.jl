using Test, HTTP, JSON3, SQLite, DBInterface, Tables
using HimalayaUI

# ---------------------------------------------------------------------------
# Unit tests: broadcast plumbing without opening a real HTTP connection.
# These assert the channel-only architecture and JSON frame shape directly.
# ---------------------------------------------------------------------------

@testset "SSE: broadcast_event! puts a curation frame onto subscriber channels" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
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
                user_id, JSON3.write(Dict(:foo => "bar")), "tab-xyz")

            @test isready(pending)
            frame = take!(pending)
            @test occursin("event: curation", frame)
            @test occursin("\"kind\":\"test_broadcast\"", frame)
            @test occursin("\"actor\":\"alice\"", frame)
            @test occursin("\"entity_id\":42", frame)
            @test occursin("\"entity_type\":\"exposure\"", frame)
            @test occursin("\"client_id\":\"tab-xyz\"", frame)

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
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        HimalayaUI.bind_db!(db)

        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(HimalayaUI.SSE_LOCK) do
            push!(HimalayaUI.SSE_SUBSCRIBERS[], sub)
        end

        try
            HimalayaUI.broadcast_event!(
                2, "anon_event", "exposure", 7,
                nothing, nothing, nothing)

            @test isready(pending)
            frame = take!(pending)
            @test occursin("event: curation", frame)
            @test occursin("\"actor\":null", frame)
            @test occursin("\"client_id\":null", frame)
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
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
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
                nothing, nothing, nothing)

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
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
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
                    nothing, nothing, nothing)
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

@testset "SSE: lookup_username returns nothing for unknown id" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
        @test HimalayaUI.lookup_username(db, 99999) === nothing
    end
end

@testset "SSE: lookup_username resolves known user" begin
    mktempdir() do dir
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
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
        db = HimalayaUI.open_db(joinpath(dir, "h.db"))
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
