using SQLite, DBInterface, HTTP, Tables

"""
    OP_LOCKS::Dict{String, ReentrantLock}
    OP_LOCKS_MU::ReentrantLock

Per-`client_op_id` lock registry used by `with_idempotency` to serialize
concurrent retries of the same mutation. `OP_LOCKS_MU` guards insertion into
`OP_LOCKS`; the per-op-id `ReentrantLock` itself guards the cache-check +
body-execution + cache-write sequence.

Entries accumulate per unique `client_op_id` for the lifetime of the process.
A TTL sweep is added in M0.7 alongside the `idempotent_responses` GC.

Lifecycle notes:

- **Single-process-safe only.** The deployment model is one-experiment-per-
  process (see CLAUDE.md). A multi-process or sharded deployment would need a
  different primitive — e.g. `INSERT OR IGNORE` on `idempotent_responses` plus
  a cached re-read — because in-process `ReentrantLock`s don't cross processes.
- **Sweep contract.** `gc_idempotent_responses!` only collects locks for
  ops that have a corresponding `idempotent_responses` row OLDER than the
  TTL — i.e. the body has executed and committed. Mid-body ops have no
  response row yet (the row is INSERTed inside the body's transaction), so
  the sweep cannot race-delete a lock another task is about to acquire
  for a never-completed op. Both the response delete and the lock delete
  happen under `OP_LOCKS_MU` so concurrent `_op_lock` callers either see
  the old lock or `get!` a fresh one — never observe a torn deletion.
"""
const OP_LOCKS    = Dict{String, ReentrantLock}()
const OP_LOCKS_MU = ReentrantLock()

function _op_lock(op_id::String)::ReentrantLock
    lock(OP_LOCKS_MU) do
        get!(OP_LOCKS, op_id, ReentrantLock())
    end
end

function _lookup_cached_response(db::SQLite.DB, op_id::String)::Union{HTTP.Response, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT status_code, body FROM idempotent_responses WHERE client_op_id = ?",
        [op_id]))
    isempty(rows) && return nothing
    row = rows[1]
    # Cached replay synthesizes a fresh Response with `Content-Type: application/json`
    # only — any other headers the original route emitted (e.g. CORS, custom
    # `Cache-Control`) are lost on retry. All M2 idempotent routes return JSON
    # only and don't add custom headers, so this is fine today; if a future
    # idempotent route needs richer headers, persist them in the cache row.
    return HTTP.Response(Int(row.status_code),
                         ["Content-Type" => "application/json"];
                         body = String(row.body))
end

"""
    with_idempotency(f, db, req) -> HTTP.Response

Stripe-style request-level idempotency keyed by the `X-Client-Op-Id` header.
Wraps a route body `f` so that:

- Without `X-Client-Op-Id`, `f` runs once and its response is returned (no
  caching, fast pass-through).
- With `X-Client-Op-Id`, the first successful response (status < 400) is
  cached in `idempotent_responses` keyed by the op-id. Subsequent calls
  with the same op-id return the cached response without re-executing `f`.
- Failures (status ≥ 400) are NOT cached; the next retry re-executes `f`.
- A per-op-id `ReentrantLock` (from `OP_LOCKS`) serializes concurrent
  retries. The cache is checked once outside the lock (lock-free fast path)
  and again inside the lock (double-check) so two racing tasks with the
  same op-id execute the body exactly once — the second sees the cached
  row written by the first.

Body is guaranteed to execute exactly once per successful op-id, even under
concurrent retry.

Single-process-safe only. The in-process `OP_LOCKS` registry doesn't cross
processes, AND the cache `INSERT` (line below) isn't `INSERT OR IGNORE` — a
multi-process retry where two processes both see no cached row and race to
INSERT will hit `UNIQUE constraint failed` and 500. To support multi-process
deployment, change the cache write to `INSERT OR IGNORE` and on `changes()=0`
re-read the cached row to return.
"""
function with_idempotency(f, db::SQLite.DB, req::HTTP.Request)
    op_id = get_client_op_id(req)
    if op_id === nothing
        # No-op-id path: still support post-commit broadcasts for routes that
        # use the M2 pattern without idempotency. We don't open a tx here, so
        # any enqueued broadcast is paired with whatever durable write the
        # body did itself; flush on success, clear on throw.
        try
            response = f()
            _flush_post_commit_broadcasts!()
            return response
        catch
            _clear_post_commit_broadcasts!()
            rethrow()
        end
    end

    # Fast path: lock-free cache check (outside any tx).
    cached = _lookup_cached_response(db, op_id)
    cached === nothing || return cached

    # Acquire per-op-id lock; re-check cache inside the lock + tx.
    return lock(_op_lock(op_id)) do
        # I2 fix: wrap the body's event writes AND the cache-row insert in a
        # single SQLite transaction so they commit or roll back together.
        # Closes the crash window where apply_event!'s event row would
        # commit but the cache row wouldn't, allowing duplicate event rows
        # on retry. Routes whose body emits events should use
        # `apply_event!(InTransaction(), db, req; ...)` to participate in
        # this transaction rather than opening a nested one.
        local response
        local replayed_cache::Bool = false
        try
            response = SQLite.transaction(db) do
                # Double-check the cache inside the lock + tx.
                cached2 = _lookup_cached_response(db, op_id)
                if cached2 !== nothing
                    replayed_cache = true
                    return cached2
                end

                resp = f()
                if resp.status < 400
                    DBInterface.execute(db,
                        "INSERT INTO idempotent_responses (client_op_id, status_code, body) VALUES (?, ?, ?)",
                        [op_id, Int(resp.status), String(copy(resp.body))])
                end
                return resp
            end
        catch
            # Rollback path: the tx threw and rolled back. Any queued
            # post-commit broadcast must be discarded — its underlying
            # writes never committed.
            _clear_post_commit_broadcasts!()
            # Issue #37 Bug 4: clean up the OP_LOCKS entry if no cache row
            # was written. gc_idempotent_responses! only collects locks
            # whose ops have an idempotent_responses row, so without this
            # cleanup, locks for permanently-failed ops would leak until
            # process restart. The cache check is defensive — under the
            # current control flow, a thrown tx body has already rolled
            # back any INSERT.
            lock(OP_LOCKS_MU) do
                if _lookup_cached_response(db, op_id) === nothing
                    delete!(OP_LOCKS, op_id)
                end
            end
            rethrow()
        end

        # On a cache replay, the body did NOT execute, so the queue should
        # already be empty. On a fresh successful body execution (status<400),
        # the queue is flushed AFTER the tx commits so subscribers can never
        # see uncommitted state. Failed-but-not-thrown bodies (status≥400)
        # have nothing cached and any speculative enqueues are dropped.
        if replayed_cache || response.status >= 400
            _clear_post_commit_broadcasts!()
        else
            _flush_post_commit_broadcasts!()
        end
        return response
    end
end

"""
    gc_idempotent_responses!(db; ttl_seconds = 3600)

Sweep `idempotent_responses` rows older than `ttl_seconds` and drop the
matching `OP_LOCKS` entries. Mid-body ops have no response row yet (the row
is INSERTed inside the body's tx), so they're invisible to the sweep — the
sweep can only collect locks for ops that have COMPLETED past their TTL.
This avoids the race where a naive sweep would delete a lock a fresh retry
is about to acquire for a never-completed op.

Holds `OP_LOCKS_MU` while collecting the to-prune set and deleting from
`OP_LOCKS`, so concurrent `_op_lock` callers either see the old lock or
`get!` a fresh one — never observe a torn deletion.

Recommended TTL: 1 hour (3600s). Long enough to cover any plausible client
retry window; short enough that long-running processes don't accumulate
unbounded state.
"""
function gc_idempotent_responses!(db::SQLite.DB; ttl_seconds::Int = 3600)
    # Single DELETE … RETURNING evaluates the cutoff once AND captures
    # exactly the set of expired client_op_ids in one shot. Avoids both:
    # (a) the SELECT/DELETE clock-skew race that an earlier two-query
    # version had (suggestion #8), and (b) the SQLite parameter limit
    # (`SQLITE_LIMIT_VARIABLE_NUMBER` is 999 on older builds) that an
    # `IN (?,?,...)`-based prune would hit on busy multiplayer sessions
    # (suggestion #10). RETURNING is supported in SQLite 3.35+, which
    # has been the system bundled version on all our deployment
    # targets since 2022.
    expired = Tables.rowtable(DBInterface.execute(db,
        """DELETE FROM idempotent_responses
           WHERE created_at < datetime('now', ?)
           RETURNING client_op_id""",
        ["-$(ttl_seconds) seconds"]))
    isempty(expired) && return nothing
    lock(OP_LOCKS_MU) do
        for r in expired
            delete!(OP_LOCKS, String(r.client_op_id))
        end
    end
    nothing
end
