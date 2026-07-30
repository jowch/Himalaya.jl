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
        # reverse proxies that close idle connections after ~60 s. Uses
        # _try_put! (same helper as broadcast_event!) so a saturated
        # subscriber drops the heartbeat frame instead of blocking the
        # Timer task — a slow client must not be able to leak heartbeat
        # tasks at one-per-15s.
        heartbeat = Timer(15; interval = 15) do _
            _try_put!(pending, ":heartbeat\n\n")
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
    register_series_routes!()
    register_picker_routes!()
    register_resolve_routes!()
    register_grouping_routes!()
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

# Per-experiment rescan timer registry (spec §9.4).
# Keyed by experiment_id (Int). Each entry is a Timer whose tick runs the
# cheap change-check + additive scan on a @spawn'd task so the libuv timer
# thread is never blocked by multi-minute scan work.
# Distinct from GC_TIMER (a single module Ref). Not persisted across restarts;
# on boot, the scheduler is re-armed by the first scan or manual rescan.
const RESCAN_TIMERS = Dict{Int, Timer}()
const RESCAN_TIMERS_MU = ReentrantLock()

# Per-experiment reentrancy guard (spec §9.4: protect with an explicit
# ReentrantLock, not the DB_WRITE_LOCK which carries write-ordering semantics).
# Keyed by experiment_id. A timer tick attempts trylock; if already locked, the
# tick is skipped (a scan is already in flight).
const RESCAN_LOCKS = Dict{Int, ReentrantLock}()
const RESCAN_LOCKS_MU = ReentrantLock()

"""
    _rescan_lock(experiment_id) -> ReentrantLock

Return the per-experiment reentrancy lock, creating it on first access.
"""
function _rescan_lock(experiment_id::Int)
    lock(RESCAN_LOCKS_MU) do
        get!(RESCAN_LOCKS, experiment_id, ReentrantLock())
    end
end

"""
    stop_rescan_scheduler!(experiment_id)

Stop and remove the rescan Timer for `experiment_id`. No-op if no timer is
registered. Called before `DELETE /api/experiments/{id}` so the timer callback
cannot fire against a non-existent row.
"""
function stop_rescan_scheduler!(experiment_id::Int)
    lock(RESCAN_TIMERS_MU) do
        t = get(RESCAN_TIMERS, experiment_id, nothing)
        if t !== nothing
            close(t)
            delete!(RESCAN_TIMERS, experiment_id)
        end
    end
    nothing
end

"""
    stop_all_rescan_timers!()

Close all per-experiment rescan timers. Called on server shutdown (mirror
`stop_gc_timer!`).
"""
function stop_all_rescan_timers!()
    lock(RESCAN_TIMERS_MU) do
        for (_, t) in RESCAN_TIMERS
            try; close(t); catch; end
        end
        empty!(RESCAN_TIMERS)
    end
    nothing
end

"""
    start_rescan_scheduler!(db, experiment_id; tick_interval_seconds, cheap_check_fn)

Arm a per-experiment rescan Timer. On each tick the libuv timer thread does the
MINIMUM possible work — it only `@spawn`s the tick body. EVERYTHING that touches
the per-experiment `ReentrantLock` (both `trylock` AND `unlock`) runs INSIDE the
spawned task:

1. `@spawn` the body so the libuv timer thread is never blocked.
2. The spawned body `trylock`s the per-experiment lock; if already held (a scan
   is in flight) it returns immediately, skipping this tick.
3. With the lock held, it calls `_rescan_tick!`, then `unlock`s in a `finally`.

CRITICAL (Julia `ReentrantLock` is task-owned): the `trylock` and the matching
`unlock` MUST execute on the SAME task. If `trylock` ran on the timer thread and
`unlock` ran on the `@spawn`'d task, the unlock would throw "unlock from wrong
thread", the lock would never release, and every later tick would be silently
skipped. That is why both calls are inside the `@spawn` block below.

`cheap_check_fn` is an optional change-detection seam that defaults to the real
`cheap_change_check` (Phase B). Task 8's backoff test injects a stub here to drive
the tier transitions deterministically without filesystem fixtures; the production
path passes nothing and the real cheap-check runs. There is no `scan_fn` seam —
the tick calls `scan_and_group!` directly (a no-op on a directory with no matching
triplets, so the test's empty temp dir needs no stub).

Idempotent: calling twice closes the previous timer before installing a new one.
"""
function start_rescan_scheduler!(db::SQLite.DB, experiment_id::Int;
                                  tick_interval_seconds::Real = 3600.0,
                                  cheap_check_fn = nothing)
    # Close any existing timer for this experiment.
    stop_rescan_scheduler!(experiment_id)

    timer = Timer(tick_interval_seconds; interval = tick_interval_seconds) do _
        # Do NOTHING that touches the per-experiment lock on the timer thread.
        # Spawn first, then trylock/unlock both on the spawned task (same-task
        # ownership — see docstring).
        Threads.@spawn begin
            lk = _rescan_lock(experiment_id)
            trylock(lk) || return  # scan already in flight; skip this tick
            try
                _rescan_tick!(db, experiment_id; cheap_check_fn = cheap_check_fn)
            catch err
                @warn "rescan tick failed" experiment_id = experiment_id exception = err
            finally
                unlock(lk)
            end
        end
    end

    lock(RESCAN_TIMERS_MU) do
        RESCAN_TIMERS[experiment_id] = timer
    end
    nothing
end

"""
    rearm_rescan_schedulers!(db)

Re-arm the per-experiment rescan timers on server boot (6c-2). The registry
(`RESCAN_TIMERS`) is in-memory, so a restart drops every scheduler; without this
an already-ingested experiment would never auto-rescan again until someone
triggers a manual scan. Arms a scheduler for each experiment whose ingest has
completed. Per-experiment failures are logged, never fatal (one bad row must not
abort boot).
"""
function rearm_rescan_schedulers!(db::SQLite.DB)
    rows = try
        Tables.rowtable(DBInterface.execute(db,
            "SELECT id FROM experiments WHERE ingest_status = 'complete'"))
    catch err
        @warn "rearm_rescan_schedulers!: query failed" exception = err
        return nothing
    end
    for r in rows
        try
            start_rescan_scheduler!(db, Int(r.id))
        catch err
            @warn "rearm_rescan_schedulers!: failed to arm" experiment_id = Int(r.id) exception = err
        end
    end
    nothing
end

"""
    _rescan_tick!(db, experiment_id; cheap_check_fn=nothing, kwargs...)

One tick of the per-experiment rescan scheduler. Called from the `@spawn`'d
body inside the Timer callback so the libuv timer thread is never blocked.

Change detection runs `cheap_check_fn(db, experiment_id)` when a seam is injected
(Task 8's backoff test), otherwise the real `cheap_change_check(db, experiment_id)`
(Phase B). On a detected change it runs `scan_and_group!(db, experiment_id)`. Both
Phase B functions resolve the experiment's `data_dir` from the row themselves, so
this tick threads no path.

Steps:
1. Read `(last_scan_tier, consecutive_empty_ticks)`; if the experiment is gone, return.
2. Determine `changed::Bool`.
3. If changed: run the scan, reset `consecutive_empty_ticks`, re-arm at the `fast` tier.
4. If not changed: increment the counter; advance fast→daily at `ticks_before_daily`,
   daily→stopped at `ticks_before_stop`, resetting the counter on each transition.
5. Persist `last_scan_tier` + `consecutive_empty_ticks` so restarts don't reset
   quiet experiments to the fast tier.

kwargs (tunable thresholds — defaults match spec §9.4; exposed for testing):
- `fast_interval`: seconds between fast-tier ticks (default 3600.0 = 1 h)
- `daily_interval`: seconds between daily-tier ticks (default 86400.0 = 24 h)
- `ticks_before_daily`: consecutive empty fast ticks before advancing to daily (default 6)
- `ticks_before_stop`: consecutive empty daily ticks before stopping (default 3)
"""
function _rescan_tick!(db::SQLite.DB, experiment_id::Int;
                        cheap_check_fn = nothing,
                        fast_interval::Real    = 3600.0,
                        daily_interval::Real   = 86400.0,
                        ticks_before_daily::Int = 6,
                        ticks_before_stop::Int  = 3)

    # Single existence + state read. Empty ⇒ experiment deleted between timer arm
    # and tick (the DELETE route stops the timer first, so this is belt-and-suspenders).
    row = Tables.rowtable(DBInterface.execute(db,
        "SELECT last_scan_tier, consecutive_empty_ticks FROM experiments WHERE id = ?",
        [experiment_id]))
    isempty(row) && return

    # Change detection: injected seam (test) or the real cheap_change_check. Both
    # take (db, experiment_id) and resolve data_dir themselves.
    changed = try
        cheap_check_fn !== nothing ? cheap_check_fn(db, experiment_id) :
                                     cheap_change_check(db, experiment_id)
    catch err
        @warn "cheap change-check failed" experiment_id = experiment_id exception = err
        false
    end

    if changed
        # Drive the corpus "analyzing" surface: phase="rescan" so the frontend maps
        # these frames to the inline ProgressBar (not the initial-scan GroupingReviewPage).
        # A terminal frame (complete/failed) MUST follow ingest_started or the surface
        # would stick — the frontend clears ingestInFlight only on a terminal frame.
        broadcast_progress!(experiment_id; kind = "ingest_started", processed = 0, total = 0, phase = "rescan")
        try
            # 3-arg on_progress: `stage` names the bar segment, `phase` selects the
            # surface. THIRD call site of this lambda shape (the two in
            # routes_experiments.jl are the others) -- keep them in step, an arity
            # mismatch here MethodErrors on the first discovered file and the catch
            # below buries it in a @warn + ingest_failed.
            scan_and_group!(db, experiment_id;
                on_progress = (p, t, stage) -> broadcast_progress!(experiment_id;
                    kind = "ingest_progress", processed = p, total = t,
                    phase = "rescan", stage = stage))
            broadcast_progress!(experiment_id; kind = "ingest_complete", processed = 0, total = 0)
        catch err
            @warn "scan failed during rescan tick" experiment_id = experiment_id exception = err
            broadcast_progress!(experiment_id; kind = "ingest_failed",
                processed = 0, total = 0, error = sprint(showerror, err))
        end
        # Re-arm at fast tier on a detected change.
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db, """
                    UPDATE experiments
                       SET last_scan_tier         = 'fast',
                           consecutive_empty_ticks = 0
                     WHERE id = ?
                """, [experiment_id])
            end
        end
        start_rescan_scheduler!(db, experiment_id;
            tick_interval_seconds = fast_interval)
        return
    end

    # No change — increment tick counter and check backoff thresholds.
    tier      = String(row[1].last_scan_tier)
    new_ticks = Int(row[1].consecutive_empty_ticks) + 1

    if tier == "fast" && new_ticks >= ticks_before_daily
        # Advance to daily
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db, """
                    UPDATE experiments
                       SET last_scan_tier         = 'daily',
                           consecutive_empty_ticks = 0
                     WHERE id = ?
                """, [experiment_id])
            end
        end
        start_rescan_scheduler!(db, experiment_id;
            tick_interval_seconds = daily_interval)
    elseif tier == "daily" && new_ticks >= ticks_before_stop
        # Stop the scheduler
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db, """
                    UPDATE experiments
                       SET last_scan_tier         = 'stopped',
                           consecutive_empty_ticks = 0
                     WHERE id = ?
                """, [experiment_id])
            end
        end
        stop_rescan_scheduler!(experiment_id)
    else
        # Stay in current tier, just increment the counter
        lock(_DB_WRITE_LOCK) do
            SQLite.transaction(db) do
                DBInterface.execute(db, """
                    UPDATE experiments
                       SET consecutive_empty_ticks = ?
                     WHERE id = ?
                """, [new_ticks, experiment_id])
            end
        end
    end
    nothing
end

function serve(db::SQLite.DB; host::String = "127.0.0.1", port::Int = 8080)
    Oxygen.resetstate()
    bind_db!(db)
    register_routes!()
    start_gc_timer!(db)
    rearm_rescan_schedulers!(db)   # 6c-2: timers are in-memory; re-arm on boot
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
    stop_all_rescan_timers!()
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
