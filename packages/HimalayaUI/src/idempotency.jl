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
- **Sweep contract for M0.7.** A sweeper cannot delete an entry from
  `OP_LOCKS` without holding `OP_LOCKS_MU` AND verifying no other task may
  still hold a reference to that lock. A naive delete races with `_op_lock`
  callers that already received the lock and are about to call `lock(it)`.
  M0.7 will likely sweep entries whose corresponding `idempotent_responses`
  row exists (i.e. the body has executed and won't be re-entered) under
  `OP_LOCKS_MU` only.
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

Single-process-safe; relies on the `idempotent_responses(client_op_id)` PK
constraint as defense-in-depth if a future multi-process deployment is
introduced (the in-process `OP_LOCKS` registry doesn't cross processes).
"""
function with_idempotency(f, db::SQLite.DB, req::HTTP.Request)
    op_id = get_client_op_id(req)
    op_id === nothing && return f()

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
        return SQLite.transaction(db) do
            # Double-check the cache inside the lock + tx.
            cached2 = _lookup_cached_response(db, op_id)
            cached2 === nothing || return cached2

            response = f()
            if response.status < 400
                DBInterface.execute(db,
                    "INSERT INTO idempotent_responses (client_op_id, status_code, body) VALUES (?, ?, ?)",
                    [op_id, Int(response.status), String(copy(response.body))])
            end
            return response
        end
    end
end

"""
    gc_idempotent_responses!(db; ttl_seconds = 3600)

Sweep `idempotent_responses` rows older than `ttl_seconds` and prune the
corresponding `OP_LOCKS` entries. Safe to call concurrently with
`with_idempotency`: the sweep holds `OP_LOCKS_MU` while reading the live set
and rebuilding the Dict, which serialises against `_op_lock`'s `get!`.

The combined pruning preserves the lock-free fast path: live ops keep their
locks; only ops whose response is gone are eligible for collection. An op
already in flight when the sweep runs will appear in the live set (because
its response has already been INSERTed by the time it returns) — so the
sweep cannot race-delete an active lock.

Recommended TTL: 1 hour (3600s). Long enough to cover any plausible client
retry window; short enough that long-running processes don't accumulate
unbounded state.
"""
function gc_idempotent_responses!(db::SQLite.DB; ttl_seconds::Int = 3600)
    DBInterface.execute(db,
        "DELETE FROM idempotent_responses WHERE created_at < datetime('now', ?)",
        ["-$(ttl_seconds) seconds"])
    lock(OP_LOCKS_MU) do
        live = Set{String}()
        rows = Tables.rowtable(DBInterface.execute(db,
            "SELECT client_op_id FROM idempotent_responses"))
        for r in rows
            push!(live, String(r.client_op_id))
        end
        for k in collect(keys(OP_LOCKS))
            k in live || delete!(OP_LOCKS, k)
        end
    end
    nothing
end
