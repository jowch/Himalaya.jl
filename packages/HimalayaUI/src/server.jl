using Oxygen
using HTTP
using SQLite
using JSON3
import Sockets

const _DB_REF = Ref{Union{SQLite.DB, Nothing}}(nothing)

# SSE subscribers. Each entry has a `pending::Channel{String}` queue. The
# handler loop blocks on `take!(pending)`; broadcast_event! `put!`s frames
# directly onto every subscriber's queue (no shared Condition needed —
# eliminates the TOCTOU race that an isready/wait combo would introduce).
const SSE_SUBSCRIBERS = Ref{Vector{Any}}([])
const SSE_LOCK        = ReentrantLock()

function current_db()
    db = _DB_REF[]
    db === nothing && error("no DB bound; call serve(db; ...) or start_test_server!")
    db
end

function bind_db!(db::SQLite.DB)
    _DB_REF[] = db
    nothing
end

function register_routes!()
    # Static frontend assets — only mounted if dist directory exists with content.
    # HIMALAYA_FRONTEND_DIST overrides the package-bundled location, supporting
    # /opt-style deployments where the build artefacts ship separately from src.
    dist_dir = get(ENV, "HIMALAYA_FRONTEND_DIST",
                   joinpath(pkgdir(HimalayaUI), "frontend", "dist"))
    if isdir(dist_dir)
        Oxygen.dynamicfiles(dist_dir, "/")
    end

    @get "/api/health" function(req::HTTP.Request)
        Dict("status" => "ok")
    end

    # Plan A: Oxygen 1.10.x @stream macro — receives an HTTP.Stream directly.
    # This is the cleanest path: no forked HTTP.serve, no response-body callback
    # gymnastics. The handler runs in its own task; we block on take!(pending)
    # so there's no busy-poll. A per-subscriber heartbeat Timer pushes
    # ":heartbeat\n\n" every 15 s onto the same channel so reverse proxies
    # don't drop idle connections. On disconnect the write throws and we prune.
    @stream "/api/events" function(stream::HTTP.Stream)
        HTTP.setheader(stream, "Content-Type"      => "text/event-stream")
        HTTP.setheader(stream, "Cache-Control"     => "no-cache")
        HTTP.setheader(stream, "Connection"        => "keep-alive")
        HTTP.setheader(stream, "X-Accel-Buffering" => "no")
        startwrite(stream)

        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(SSE_LOCK) do
            push!(SSE_SUBSCRIBERS[], sub)
        end

        # Per-subscriber heartbeat: keeps the TCP connection alive through
        # reverse proxies that close idle connections after ~60 s.
        heartbeat = Timer(15; interval = 15) do _
            try put!(pending, ":heartbeat\n\n") catch end
        end

        try
            for frame in pending
                write(stream, frame)
            end
        finally
            close(heartbeat)
            lock(SSE_LOCK) do
                filter!(x -> x !== sub, SSE_SUBSCRIBERS[])
            end
        end
        nothing
    end

    register_users_routes!()
    register_experiments_routes!()
    register_samples_routes!()
    register_exposures_routes!()
    register_peaks_routes!()
    register_messages_routes!()
    register_trace_routes!()
    register_analysis_routes!()
    register_export_routes!()
end

function serve(db::SQLite.DB; host::String = "127.0.0.1", port::Int = 8080)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    Oxygen.serve(; host, port, show_banner = false, docs = false, metrics = false)
end

function start_test_server!(db::SQLite.DB, port::Int)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    Oxygen.serve(; host = "127.0.0.1", port, async = true, show_banner = false, docs = false, metrics = false)
    _wait_for_server(port)
end

function stop_test_server!()
    Oxygen.terminate()
    Oxygen.resetstate()
    _DB_REF[] = nothing
    SSE_SUBSCRIBERS[] = []
    nothing
end

function find_free_port()
    server = Sockets.listen(Sockets.IPv4(0), 0)
    port   = Sockets.getsockname(server)[2]
    close(server)
    Int(port)
end

function _wait_for_server(port::Int; timeout_s = 5.0)
    deadline = time() + timeout_s
    while time() < deadline
        try
            resp = HTTP.get("http://127.0.0.1:$port/api/health";
                            connect_timeout = 1, readtimeout = 1,
                            retry = false, status_exception = false)
            resp.status == 200 && return
        catch
        end
        sleep(0.05)
    end
    error("test server on port $port did not become ready within $(timeout_s)s")
end
