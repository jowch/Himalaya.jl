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
# across a `_DB_WRITE_LOCK` acquisition. Every non-GET request already holds
# _DB_WRITE_LOCK (db_connection_middleware) when `with_idempotency` briefly
# takes OP_LOCKS_MU to `get!` its `_op_lock`, so the per-op lock is redundant
# for HTTP callers. `gc_idempotent_responses!` does _DB_WRITE_LOCK → release
# → OP_LOCKS_MU. Future code added under OP_LOCKS_MU must not call into
# anything that may acquire _DB_WRITE_LOCK, or the orderings will diverge
# into a deadlock.
const _DB_WRITE_LOCK = ReentrantLock()

# SSE subscribers. Each entry has a `pending::Channel{String}` queue. The
# handler loop blocks on `take!(pending)`; broadcast_event! `put!`s frames
# directly onto every subscriber's queue (no shared Condition needed —
# eliminates the TOCTOU race that an isready/wait combo would introduce).
const SSE_SUBSCRIBERS = Ref{Vector{Any}}([])
const SSE_LOCK        = ReentrantLock()
const SSE_OPENING_FRAME   = "retry: 3000\n\n"
const SSE_HEARTBEAT_FRAME = "event: heartbeat\ndata: {}\n\n"

function current_db()
    reader = _scoped_reader()
    reader === nothing || return reader
    db = _DB_REF[]
    db === nothing && error("no DB bound; call serve(db; ...) or start_test_server!")
    db
end

function bind_db!(db::SQLite.DB)
    _DB_REF[] = db
    _reset_reader_pool!(db)
    nothing
end

function unbind_db!()
    _DB_REF[] = nothing
    _reset_reader_pool!(nothing)
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
    # so there's no busy-poll. On disconnect the write throws and we prune.
    @stream "/api/events" function(stream::HTTP.Stream)
        HTTP.setheader(stream, "Content-Type"      => "text/event-stream")
        HTTP.setheader(stream, "Cache-Control"     => "no-cache")
        HTTP.setheader(stream, "Connection"        => "keep-alive")
        HTTP.setheader(stream, "X-Accel-Buffering" => "no")
        startwrite(stream)
        write(stream, SSE_OPENING_FRAME)

        pending = Channel{String}(64)
        sub = (pending = pending,)
        lock(SSE_LOCK) do
            push!(SSE_SUBSCRIBERS[], sub)
        end

        # Per-subscriber heartbeat: keeps the TCP connection alive through
        # reverse proxies that close idle connections after ~60 s. Uses
        # _try_put! (same helper as broadcast_event!) so a saturated
        # subscriber drops the heartbeat frame instead of blocking the
        # Timer task — a slow client must not be able to leak heartbeat
        # tasks at one-per-15s.
        heartbeat = Timer(15; interval = 15) do _
            _try_put!(pending, SSE_HEARTBEAT_FRAME)
        end

        drain = Threads.@spawn :interactive for frame in pending
            write(stream, frame)
        end
        try
            wait(drain)
        catch err
            err isa TaskFailedException ? throw(err.task.exception) : rethrow()
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
    register_trace_routes!()
    register_analysis_routes!()
    register_export_routes!()
    register_picker_routes!()
    register_resolve_routes!()
    register_bridge_routes!()
    register_label_routes!()
    register_tag_color_routes!()
    register_fs_routes!()
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

const INGEST_LOCKS = Dict{Int, ReentrantLock}()
const INGEST_LOCKS_MU = ReentrantLock()

_ingest_lock(experiment_id::Int) = lock(INGEST_LOCKS_MU) do
    get!(INGEST_LOCKS, experiment_id, ReentrantLock())
end

ingest_running(experiment_id::Int) = islocked(_ingest_lock(experiment_id))

function _set_ingest_status!(db::SQLite.DB, experiment_id::Int, status::AbstractString;
        scanned::Bool = false, error::Union{String,Nothing} = nothing)
    lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            if scanned
                DBInterface.execute(db,
                    "UPDATE experiments SET ingest_status = ?, ingest_error = NULL, last_scanned_at = ? WHERE id = ?",
                    [status, format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"), experiment_id])
            elseif error !== nothing
                DBInterface.execute(db, "UPDATE experiments SET ingest_status = ?, ingest_error = ? WHERE id = ?",
                    [status, first(error, 2000), experiment_id])
            else
                DBInterface.execute(db, "UPDATE experiments SET ingest_status = ? WHERE id = ?",
                    [status, experiment_id])
            end
        end
    end
end

"""
    run_ingest!(db, experiment_id; phase=nothing) -> Symbol

Bring one experiment up to date with its folder: back-fill beam monitors, then
`scan_and_group!`, which inserts new exposures and analyzes every exposure that
is new, never analyzed, whose trace has since appeared, or whose analysis is from
an older `ANALYSIS_VERSION`. Every trigger (create, Rescan, opening the
experiment, boot) goes through here, so one experiment never has two passes at
once: a second caller waits, then runs its own pass, which is a cheap no-op when
nothing changed in between. Returns `:ok` or `:failed`. Streams the ingest SSE frames; `phase`
is `"rescan"` for every pass after the first.
"""
function run_ingest!(db::SQLite.DB, experiment_id::Int; phase::Union{String,Nothing} = nothing, force::Bool = false)
    lk = _ingest_lock(experiment_id)
    lock(lk)
    try
        result = with_background_db(db) do db
            extra = phase === nothing ? (;) : (; phase)
            _set_ingest_status!(db, experiment_id, "scanning")
            broadcast_progress!(experiment_id; kind = "ingest_started", processed = 0, total = 0, extra...)
            try
                try_backfill_monitors!(db, experiment_id)
                scan_and_group!(db, experiment_id; retry_failed = force,
                    on_progress = (p, t, stage) -> broadcast_progress!(experiment_id;
                        kind = "ingest_progress", processed = p, total = t, stage = stage, extra...))
                suggest_labels_quietly!(db, experiment_id)
                _set_ingest_status!(db, experiment_id, "complete"; scanned = true)
                broadcast_progress!(experiment_id; kind = "ingest_complete", processed = 0, total = 0)
                :ok
            catch err
                @warn "ingest failed" experiment_id exception = err
                message = sprint(showerror, err)
                _set_ingest_status!(db, experiment_id, "failed"; error = message)
                broadcast_progress!(experiment_id; kind = "ingest_failed",
                    processed = 0, total = 0, error = message)
                :failed
            end
        end
        result === :ok && prewarm_after_ingest!(db, experiment_id)
        result
    finally
        unlock(lk)
    end
end

const THUMB_PREWARMS = Dict{Int,Task}()
const _THUMB_PREWARMS_LOCK = ReentrantLock()
const _THUMB_PREWARM_GATE = Ref{Union{Nothing,Channel{Nothing}}}(nothing)

function prewarm_after_ingest!(db::SQLite.DB, experiment_id::Int)
    thumb_cache_dir(db) === nothing && return nothing
    gate = _THUMB_PREWARM_GATE[]
    lock(_THUMB_PREWARMS_LOCK) do
        previous = get(THUMB_PREWARMS, experiment_id, nothing)
        task = Threads.@spawn try
            previous === nothing || wait(previous)
            gate === nothing || take!(gate)
            with_background_db(db) do bg
                prewarm_thumbnails!(bg; threads = false, overwrite = false, experiment_id)
            end
        catch e
            @warn "thumbnail prewarm failed" experiment_id exception = e
        end
        THUMB_PREWARMS[experiment_id] = task
    end
end

function wait_prewarm(experiment_id::Int)
    task = lock(() -> get(THUMB_PREWARMS, experiment_id, nothing), _THUMB_PREWARMS_LOCK)
    task === nothing || wait(task)
    nothing
end

function _ingest_phase(db::SQLite.DB, experiment_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT last_scanned_at FROM experiments WHERE id = ?", [experiment_id]))
    (isempty(rows) || ismissing(first(rows).last_scanned_at)) ? nothing : "rescan"
end

"""
    resume_ingests!(db)

Finish whatever a restart or an upgrade left undone: every experiment still
marked `scanning`, or with work by `ingest_has_work` (new files, traces that
have appeared, exposures never analyzed or analyzed by an older version), gets
one `run_ingest!` pass. Called from `start_background!` on a spawned task.
"""
function resume_ingests!(db::SQLite.DB)
    rows = Tables.rowtable(DBInterface.execute(db, "SELECT id, ingest_status FROM experiments ORDER BY id"))
    for r in rows
        id = Int(r.id)
        try
            stuck = coalesce(r.ingest_status, "") == "scanning"
            (stuck || ingest_has_work(db, id)) || continue
            run_ingest!(db, id; phase = _ingest_phase(db, id))
        catch err
            @warn "resume_ingests!: experiment skipped" experiment_id = id exception = err
        end
    end
    nothing
end

function try_backfill_monitors!(db::SQLite.DB, experiment_id::Int)
    try
        backfill_monitors!(db, experiment_id)
    catch err
        @warn "beam monitor backfill failed" experiment_id = experiment_id exception = err
        0
    end
end

function backfill_all_monitors!(db::SQLite.DB)
    ids = Tables.rowtable(DBInterface.execute(db,
        "SELECT DISTINCT experiment_id FROM exposures
          WHERE prp_path IS NOT NULL AND monitors_read = 0 AND experiment_id IS NOT NULL"))
    for r in ids
        try_backfill_monitors!(db, Int(r.experiment_id))
    end
    nothing
end

function start_background!(db::SQLite.DB; resume::Bool = true)
    start_gc_timer!(db)
    Threads.@spawn try
        with_background_db(db) do bg
            backfill_all_monitors!(bg)
            resume && resume_ingests!(bg)
        end
    catch err
        @warn "boot background work failed" exception = err
    end
end

function boot!(db::SQLite.DB; resume::Bool = true)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    start_background!(db; resume)
end

function serve(db::SQLite.DB; host::String = "127.0.0.1", port::Int = 8080)
    boot!(db)
    # parallel = true wraps each request in `Threads.@spawn`, dispatching
    # handlers across the default thread pool. Without it, every handler
    # runs cooperatively on HTTP.jl's single interactive thread — even with
    # JULIA_NUM_THREADS > 1 — so concurrent requests serialize. See #115.
    #
    # `db_connection_middleware` (db_connections.jl) gives GET/HEAD a pooled
    # reader with one snapshot each, runs the event stream with neither, and
    # holds `_DB_WRITE_LOCK` around every other request, so only one task at
    # a time uses the `_DB_REF` singleton (#122 Races 1 and 2).
    Oxygen.serve(; host, port, show_banner = false, docs = false, metrics = false,
                 parallel = true, middleware = [db_connection_middleware])
end

function start_test_server!(db::SQLite.DB, port::Int)
    boot!(db; resume = false)
    # Match production's threading model so the test suite exercises the
    # same scheduling as `serve`. Any test that depends on single-threaded
    # request ordering is a bug.
    Oxygen.serve(; host = "127.0.0.1", port, async = true,
                 show_banner = false, docs = false, metrics = false,
                 parallel = true, middleware = [db_connection_middleware])
    _wait_for_server(port)
end

function stop_test_server!()
    stop_gc_timer!()
    Oxygen.terminate()
    Oxygen.resetstate()
    unbind_db!()
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

function _wait_for_server(port::Int; timeout_s = 60.0)
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
