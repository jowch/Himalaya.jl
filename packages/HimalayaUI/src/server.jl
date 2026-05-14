using Oxygen
using HTTP
using SQLite
using JSON3
import Sockets

const _DB_REF = Ref{Union{SQLite.DB, Nothing}}(nothing)

# Issue #122. Serializes every `SQLite.transaction(db)` call that targets the
# singleton `_DB_REF` connection. Without this, two concurrent writers race on
# `SQLite.transaction`'s TOCTOU: both can pass `intransaction(db) == false`
# and attempt `BEGIN` (loud 500), or one passes false/the other true and
# silently nests a SAVEPOINT inside the other's tx (silent corruption —
# Failure mode B in #122). The lock also serializes `db.stmt_wrappers` Dict
# mutation between concurrent writers (Race 2 writer-vs-writer).
#
# ReentrantLock so nested write paths in the SAME task — e.g. a route body
# inside `with_idempotency` that calls `analyze_exposure!` which opens its
# own tx — don't self-deadlock.
#
# Lock-ordering invariant: `OP_LOCKS_MU` (idempotency.jl) is never held
# across a `_DB_WRITE_LOCK` acquisition. Today `with_idempotency` does
# OP_LOCKS_MU → release → _DB_WRITE_LOCK (`_op_lock` takes the mutex only to
# `get!` from the Dict, then releases before the lock body runs), while
# `gc_idempotent_responses!` does _DB_WRITE_LOCK → release → OP_LOCKS_MU
# (the DELETE commits before the OP_LOCKS prune runs). Future code added
# under OP_LOCKS_MU must not call into anything that may acquire
# _DB_WRITE_LOCK, or the orderings will diverge into a deadlock.
const _DB_WRITE_LOCK = ReentrantLock()

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

        # SPA catch-all (spec §3.2). Doublestar `/**` is HTTP.jl's multi-
        # segment wildcard (Handlers.jl:174,219–227); single-conditional
        # captures like `{rest:.*}` only match one segment. The `api/`
        # guard is critical — without it an unregistered API route would
        # fall through to index.html and mask 404s as 200s.
        @get "/**" function(req::HTTP.Request)
            path = HTTP.URI(req.target).path
            rest = lstrip(path, '/')
            (startswith(rest, "api/") || rest == "api") && return HTTP.Response(404, "Not found")
            return HTTP.Response(200,
                ["Content-Type" => "text/html; charset=utf-8",
                 "Cache-Control" => "no-store"],
                read(joinpath(dist_dir, "index.html")))
        end
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
    register_comparisons_routes!()
    register_picker_routes!()
    register_resolve_routes!()
end

const GC_TIMER = Ref{Union{Timer, Nothing}}(nothing)

"""
    start_gc_timer!(db; interval_seconds = 1800, ttl_seconds = 3600)

Schedule periodic invocations of `gc_idempotent_responses!`. The first
invocation fires immediately on the timer thread; subsequent invocations
fire every `interval_seconds`. Idempotent: calling twice closes the
previous timer before installing a new one.
"""
function start_gc_timer!(db::SQLite.DB; interval_seconds::Real = 30 * 60,
                                        ttl_seconds::Int = 3600)
    GC_TIMER[] !== nothing && (close(GC_TIMER[]); GC_TIMER[] = nothing)
    GC_TIMER[] = Timer(0.0; interval = interval_seconds) do _
        try
            gc_idempotent_responses!(db; ttl_seconds = ttl_seconds)
        catch err
            @warn "idempotent_responses GC sweep failed" exception = err
        end
    end
    nothing
end

"""
    stop_gc_timer!()

Close the active GC timer if any. No-op when no timer was started.
"""
function stop_gc_timer!()
    if GC_TIMER[] !== nothing
        close(GC_TIMER[])
        GC_TIMER[] = nothing
    end
    nothing
end

function serve(db::SQLite.DB; host::String = "127.0.0.1", port::Int = 8080)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    start_gc_timer!(db)
    # parallel = true wraps each request in `Threads.@spawn`, dispatching
    # handlers across the default thread pool. Without it, every handler
    # runs cooperatively on HTTP.jl's single interactive thread — even with
    # JULIA_NUM_THREADS > 1 — so concurrent requests serialize. See #115.
    #
    # All routes still share one `_DB_REF` connection. Writes serialize at
    # the Julia level via `_DB_WRITE_LOCK` (above) — every `SQLite.transaction`
    # site on the singleton (`with_idempotency`, default `apply_event!`,
    # `persist_analysis!`, `analyze_exposure!`, the `gc_idempotent_responses!`
    # DELETE) acquires it. This closes the `SQLite.transaction` TOCTOU race
    # (#122 Race 1: loud 500s + silent savepoint nesting) and the
    # writer-vs-writer `db.stmt_wrappers` Dict race (#122 Race 2).
    #
    # Residual: reader-vs-writer `stmt_wrappers` mutation remains possible
    # because reads don't take the lock (and shouldn't — WAL is the whole
    # point of #115). Empirically rare; surfaces as intermittent test flakes
    # or sporadic 500s. The real fix is per-request reader connections;
    # tracked as a follow-up to #122.
    Oxygen.serve(; host, port, show_banner = false, docs = false, metrics = false,
                 parallel = true)
end

function start_test_server!(db::SQLite.DB, port::Int)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    # Match production's threading model so the test suite exercises the
    # same scheduling as `serve`. Any test that depends on single-threaded
    # request ordering is a bug.
    Oxygen.serve(; host = "127.0.0.1", port, async = true,
                 show_banner = false, docs = false, metrics = false,
                 parallel = true)
    _wait_for_server(port)
end

function stop_test_server!()
    stop_gc_timer!()
    Oxygen.terminate()
    Oxygen.resetstate()
    _DB_REF[] = nothing
    SSE_SUBSCRIBERS[] = []
    empty!(OP_LOCKS)
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
