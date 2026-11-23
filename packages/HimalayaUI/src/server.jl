using Oxygen
using HTTP
using SQLite
using JSON3
import Sockets

const _DB_REF = Ref{Union{SQLite.DB, Nothing}}(nothing)

const _DB_WRITE_LOCK = ReentrantLock()

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
    dist_dir = get(ENV, "HIMALAYA_FRONTEND_DIST",
                   joinpath(pkgdir(HimalayaUI), "frontend", "dist"))
    if isdir(dist_dir)
        Oxygen.dynamicfiles(dist_dir, "/")

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
    register_setups_routes!()
    register_samples_routes!()
    register_exposures_routes!()
    register_peaks_routes!()
    register_trace_routes!()
    register_analysis_routes!()
    register_export_routes!()
    register_bridge_routes!()
    register_compare_routes!()
    register_label_routes!()
    register_label_color_routes!()
    register_fs_routes!()
end

const GC_TIMER = Ref{Union{Timer, Nothing}}(nothing)

function start_gc_timer!(db::SQLite.DB; interval_seconds::Real = 30 * 60,
                                        ttl_seconds::Int = 3600)
    GC_TIMER[] !== nothing && (close(GC_TIMER[]); GC_TIMER[] = nothing)
    GC_TIMER[] = Timer(_ -> _gc_sweep(db, ttl_seconds), 0.0; interval = interval_seconds)
    nothing
end

function _gc_sweep(db::SQLite.DB, ttl_seconds::Int)
    try
        gc_idempotent_responses!(db; ttl_seconds = ttl_seconds)
    catch err
        @warn "idempotent_responses GC sweep failed" exception = err
    end
    nothing
end

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
                add_new_setups!(db, experiment_id)
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
                prewarm_thumbnails!(bg; threads = false, experiment_id)
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
    rows = sqlrows(DBInterface.execute(db,
        "SELECT last_scanned_at FROM experiments WHERE id = ?", [experiment_id]))
    (isempty(rows) || ismissing(first(rows).last_scanned_at)) ? nothing : "rescan"
end

function resume_ingests!(db::SQLite.DB)
    rows = sqlrows(DBInterface.execute(db, "SELECT id, ingest_status FROM experiments ORDER BY id"))
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
    ids = sqlrows(DBInterface.execute(db,
        "SELECT DISTINCT experiment_id FROM exposures
          WHERE prp_path IS NOT NULL AND monitors_read = 0 AND experiment_id IS NOT NULL"))
    for r in ids
        eid = Int(r.experiment_id)
        try_backfill_monitors!(db, eid) > 0 && suggest_labels_quietly!(db, eid)
    end
    nothing
end

function start_background!(db::SQLite.DB; resume::Bool = true)
    start_gc_timer!(db)
    Threads.@spawn _boot_background_work(db, resume)
end

function _boot_background_work(db::SQLite.DB, resume::Bool)
    try
        with_background_db(db) do bg
            backfill_all_monitors!(bg)
            resume && resume_ingests!(bg)
        end
    catch err
        @warn "boot background work failed" exception = err
    end
    nothing
end

function boot!(db::SQLite.DB; resume::Bool = true)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    start_background!(db; resume)
end

function warn_if_single_threaded(n::Int = Threads.nthreads(:default))
    n >= 2 && return false
    @warn "The server has one worker thread, so a scan or an analysis holds every other request until it finishes. Start Julia with --threads=4,1, or use bin/himalaya, which does." threads = n
    true
end

function serve(db::SQLite.DB; host::String = "127.0.0.1", port::Int = 8080)
    warn_if_single_threaded()
    boot!(db)
    Oxygen.serve(; host, port, show_banner = false, docs = false, metrics = false,
                 parallel = true, middleware = [db_connection_middleware])
end

function start_test_server!(db::SQLite.DB, port::Int)
    boot!(db; resume = false)
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
