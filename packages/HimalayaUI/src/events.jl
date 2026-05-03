using JSON3, SQLite, DBInterface, HTTP, Tables
using Dates: now, UTC, format, @dateformat_str

"""
    apply_event!(db, req; kind, entity_type, entity_id, payload, undoes_event_id=nothing)
      -> NamedTuple{(:event_id, :view_row_id), Tuple{Int, Union{Int, Nothing}}}

Atomic event-append + view-update. The log and the views must move together
or neither moves. Returns a named tuple with two fields:
- `event_id`: the newly-inserted event id in user_actions.
- `view_row_id`: the id of the view row inserted by the dispatcher, or
  `nothing` for non-insert dispatcher branches (DELETE, no-op, etc.).
  Callers that need the inserted row id (e.g. POST /peaks) use this directly
  instead of re-querying, eliminating a read-back race with concurrent writers.

`payload` is any JSON-serializable Dict / NamedTuple / nothing. If nothing,
the event is recorded but no view update fires (use sparingly — most actions
should carry a payload).
"""
function apply_event!(db::SQLite.DB, req;
                      kind::String,
                      entity_type::String,
                      entity_id::Integer,
                      payload = nothing,
                      undoes_event_id::Union{Int,Nothing} = nothing,
                      defer_broadcast::Bool = false,
                      post_state::Union{Dict, Nothing} = nothing)
    username = get_username(req)
    client_id = get_client_id(req)
    client_op_id = get_client_op_id(req)
    user_id  = username === nothing ? nothing : get_or_create_user!(db, username)
    payload_json = payload === nothing ? nothing : JSON3.write(payload)

    view_row_id_ref = Ref{Union{Int, Nothing}}(nothing)
    event_id = SQLite.transaction(db) do
        res = DBInterface.execute(db,
            """INSERT INTO user_actions
               (user_id, action, entity_type, entity_id, payload, undoes_event_id, client_id, client_op_id)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            [user_id, kind, entity_type, Int(entity_id), payload_json, undoes_event_id, client_id, client_op_id])
        eid = Int(DBInterface.lastrowid(res))

        # Canonicalize payload before dispatch — round-trip through JSON3 so
        # the live dispatcher and rebuild_views_from_log! see exactly the
        # same shape (JSON3.Object that supports both .field and [:field]
        # access, eliminating the Symbol-key vs String-key footgun).
        if payload_json !== nothing
            payload_canonical = JSON3.read(payload_json)
            view_row_id_ref[] = update_view_for_event!(db, kind, entity_id, payload_canonical, eid)
        end
        eid
    end

    # Best-effort SSE broadcast — fires AFTER the transaction commits, so a
    # subscriber never sees an event that was rolled back. If the process
    # dies between commit and broadcast, the event is durable in user_actions
    # but the frame is lost; clients reconcile on reconnect via TanStack
    # Query refetch (see R5a).
    # Defense-in-depth: broadcast is best-effort; if broadcast_event! is ever
    # detached (e.g. in a stripped-down deployment that doesn't ship the SSE
    # endpoint), the guard prevents an UndefVarError. Today broadcast_event!
    # is always defined alongside apply_event! — the try/catch below catches
    # runtime issues; this guard catches definition-time issues.
    if !defer_broadcast && isdefined(@__MODULE__, :broadcast_event!)
        # M0.4: suppress SSE broadcast for analyze_run no-ops (both skip flags true).
        # M2 wires synchronous reanalyze inside curation routes; without this guard
        # every curation event would also fan out an analyze_run frame even when
        # nothing changed — O(N) extra frames per session. The user_actions row is
        # still written; only the broadcast is suppressed. Strict `=== true` guards
        # against the JSON3.Object case where a missing key would return `nothing`.
        suppress = kind == "analyze_run" &&
                   payload !== nothing &&
                   get(payload, :findpeaks_skipped, false) === true &&
                   get(payload, :indexpeaks_skipped, false) === true
        if !suppress
            try
                broadcast_event!(event_id, kind, entity_type, Int(entity_id),
                                 user_id, client_id, client_op_id, payload_json;
                                 post_state = post_state)
            catch err
                @warn "broadcast_event! failed (event still durable in user_actions)" exception=err
            end
        end
    end
    (event_id = event_id, view_row_id = view_row_id_ref[])
end

"""
    _system_request() -> HTTP.Request

Synthetic request with no `X-Username` so the resulting event's `user_id` is
NULL — pipeline runs aren't attributed to any actor. Used by `analyze_run`
events emitted from `analyze_exposure!`.
"""
_system_request() = HTTP.Request("INTERNAL", "/_system", Pair{String,String}[], UInt8[])

"""
Dispatcher that updates materialized views based on event kind.
Currently a no-op for most events; populated as routes migrate to apply_event!.

**Payload contract:** `payload` is normalized to a `JSON3.Object` (or `nothing`)
before this is called. The live path in `apply_event!` round-trips through
`JSON3.write` → `JSON3.read` so dispatcher branches see the same shape they'd
see during `rebuild_views_from_log!`. This eliminates the Symbol-key vs
String-key footgun: `JSON3.Object` supports both `obj.q` and `obj[:q]` /
`obj["q"]`. Branch code accesses fields uniformly via `payload.q` style.
"""
function update_view_for_event!(db, kind, entity_id, payload, event_id)
    # R4.2 dispatcher branches — one per view-producing curation kind.
    # All writes happen inside the transaction opened by apply_event!.
    # INSERT branches return the lastrowid of their inserted row as Union{Int,Nothing};
    # non-INSERT branches (DELETE, no-op) return nothing.

    if kind == "peak_added"
        res = DBInterface.execute(db,
            """INSERT INTO peak_curations (exposure_id, kind, q, created_by)
               VALUES (?, 'add', ?, (SELECT user_id FROM user_actions WHERE id = ?))""",
            [Int(entity_id), Float64(payload.q), event_id])
        return Int(DBInterface.lastrowid(res))
    end

    if kind == "peak_excluded"
        res = DBInterface.execute(db,
            """INSERT INTO peak_curations (exposure_id, kind, q, created_by)
               VALUES (?, 'exclude', ?, (SELECT user_id FROM user_actions WHERE id = ?))""",
            [Int(entity_id), Float64(payload.q), event_id])
        return Int(DBInterface.lastrowid(res))
    end

    if kind == "peak_unexcluded"
        # payload.q is the auto peak's q; remove the matching exclude curation.
        # Tolerance shape mirrors effective_peaks.
        DBInterface.execute(db,
            """DELETE FROM peak_curations
               WHERE exposure_id = ? AND kind = 'exclude'
                 AND ABS(q - ?) <= MAX(1e-6, ABS(?) * 0.001)""",
            [Int(entity_id), Float64(payload.q), Float64(payload.q)])
        return nothing
    end

    if kind == "index_confirmed"
        DBInterface.execute(db,
            """INSERT OR IGNORE INTO index_group_members (group_id, index_id)
               VALUES (?, ?)""",
            [Int(payload.group_id), Int(payload.index_id)])
        return nothing
    end

    if kind == "index_unconfirmed"
        DBInterface.execute(db,
            """DELETE FROM index_group_members
               WHERE group_id = ? AND index_id = ?""",
            [Int(payload.group_id), Int(payload.index_id)])
        return nothing
    end

    # Scaffolding / legacy:
    kind == "noop_test" && return nothing
    # default: no view update (analyze_run and other instrumentation events land here)
    nothing
end

"""
    lookup_username(db, user_id) -> Union{String, Nothing}

Resolve a user_id to its username string. Returns nothing for NULL user_id.
Used by broadcast_event! to format the `actor` field in SSE frames so
clients can self-echo-filter their own edits.
"""
function lookup_username(db::SQLite.DB, user_id::Integer)::Union{String, Nothing}
    rows = Tables.rowtable(DBInterface.execute(db,
        "SELECT username FROM users WHERE id = ?", [Int(user_id)]))
    isempty(rows) || ismissing(rows[1].username) ? nothing : String(rows[1].username)
end

"""
    _try_put!(ch, value) -> Bool

Non-blocking put! for SSE subscriber channels. Returns true if the value was
enqueued, false if the channel is closed or full (slow subscriber).

A closed channel means the subscriber disconnected; a full channel means the
subscriber's take! loop has fallen behind. Both cases are treated as dead —
the SSE design is best-effort with EventSource auto-reconnect + TanStack Query
refetch providing reconciliation.

TOCTOU note: the heartbeat Timer for this subscriber could race between the
n_avail check and the put!, in theory filling the last slot and causing put!
to block. In practice the heartbeat fires at most once per 15 s; with cap=64
an idle subscriber would need 16+ minutes of heartbeats to fill. The race
window is negligible; Option A (document + ship) is the right tradeoff here.
"""
function _try_put!(ch::Channel{String}, value::String)::Bool
    isopen(ch) || return false
    # Channel is full → subscriber is too slow. Skip the put rather than block;
    # the SSE design is best-effort with reconnect-driven reconciliation.
    Base.n_avail(ch) >= ch.sz_max && return false
    put!(ch, value)
    return true
end

"""
    broadcast_event!(event_id, kind, entity_type, entity_id, user_id, client_id, client_op_id, payload_json)

Format a single SSE frame and enqueue it onto every subscriber's pending
channel. The frame carries `client_op_id` (the per-mutation idempotency key
echoed from the originating request's `X-Client-Op-Id` header) and `ts`
(the server-side broadcast timestamp, ISO-8601 UTC). Closed channels
(disconnected clients) and full channels (slow subscribers) are pruned —
the client will reconnect via EventSource auto-reconnect and refetch via
TanStack Query.

`client_id` is the per-tab SSE routing identity (from the `X-Client-Id`
request header) embedded in the frame so subscribers can self-echo-filter
events that originated in their own tab. `nothing` for system-emitted
events (no originating request).

Best-effort: this fires AFTER apply_event!'s transaction commits, so a
subscriber never sees an event that was rolled back. If the process dies
between commit and broadcast, the event is durable in user_actions but the
frame is lost; clients reconcile on reconnect via TanStack Query refetch.

SSE_SUBSCRIBERS and SSE_LOCK live in server.jl but are visible here because
both files are included into the same HimalayaUI module.
"""
function broadcast_event!(event_id::Integer, kind::String, entity_type::String,
                          entity_id::Integer, user_id::Union{Integer, Nothing},
                          client_id::Union{String, Nothing},
                          client_op_id::Union{String, Nothing},
                          payload_json::Union{String, Nothing};
                          post_state::Union{Dict, Nothing} = nothing)
    actor = user_id === nothing ? nothing : lookup_username(current_db(), user_id)
    fields = Dict{Symbol, Any}(
        :id           => Int(event_id),
        :kind         => kind,
        :entity_type  => entity_type,
        :entity_id    => Int(entity_id),
        :actor        => actor,
        :client_id    => client_id,
        :client_op_id => client_op_id,
        :ts           => format(now(UTC), dateformat"yyyy-mm-ddTHH:MM:SS.sssZ"),
        :payload      => payload_json === nothing ? nothing : JSON3.read(payload_json),
    )
    post_state === nothing || (fields[:post_state] = post_state)
    msg = JSON3.write(fields)
    frame = "event: curation\ndata: $msg\n\n"
    lock(SSE_LOCK) do
        to_drop = []
        for sub in SSE_SUBSCRIBERS[]
            _try_put!(sub.pending, frame) || push!(to_drop, sub)
        end
        for sub in to_drop
            filter!(x -> x !== sub, SSE_SUBSCRIBERS[])
        end
    end
    nothing
end

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

    # Fast path: lock-free cache check.
    cached = _lookup_cached_response(db, op_id)
    cached === nothing || return cached

    # Acquire per-op-id lock; re-check cache inside the lock.
    return lock(_op_lock(op_id)) do
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

"""
    rebuild_views_from_log!(db, exposure_id) -> Nothing

Re-fold every event for `exposure_id` from `user_actions` into the materialized
view tables (`peak_curations`, `index_group_members`). Used for migration,
disaster recovery, and replay testing.

The dispatcher (`update_view_for_event!`) is the sole writer to view tables;
this function exercises that contract by calling it for every historical
event in id-order. Tested via property: starting from empty views, applying
every event in order produces the same state as live `apply_event!` does.

For R4.1 this is a stub — R4.2 fills in dispatcher branches per kind, and
the property test asserts the round-trip.
"""
function rebuild_views_from_log!(db::SQLite.DB, exposure_id::Int)
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, action, entity_id, payload
           FROM user_actions
           WHERE entity_type = 'exposure' AND entity_id = ?
           ORDER BY id""", [exposure_id]))
    for r in rows
        kind = String(r.action)
        eid  = Int(r.id)
        ent  = Int(r.entity_id)
        ismissing(r.payload) && continue
        payload = JSON3.read(String(r.payload))
        update_view_for_event!(db, kind, ent, payload, eid)
    end
    nothing
end
