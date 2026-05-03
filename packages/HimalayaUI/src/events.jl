using JSON3, SQLite, DBInterface, HTTP, Tables
using Dates: now, UTC, format, @dateformat_str

"""
    InTransaction

Sentinel singleton type indicating the caller has already opened a SQLite
transaction. Use the `apply_event!(::InTransaction, db, req; ...)` method
to participate in that transaction rather than opening a nested one. This
exists so `with_idempotency` can wrap the entire event-write + cache-write
sequence in a single atomic transaction — closing the I2 crash window where
the event row would commit but the cache row wouldn't, allowing duplicate
events on retry.
"""
struct InTransaction end

"""
    apply_event!(db, req; kind, entity_type, entity_id, payload, undoes_event_id=nothing)
      -> NamedTuple{(:event_id, :view_row_id), Tuple{Int, Union{Int, Nothing}}}

Atomic event-append + view-update. The log and the views must move together
or neither moves. Returns a named tuple with two fields:
- `event_id`: the newly-inserted event id in user_actions (or, on idempotent
  retry, the id of the prior event row with the same (client_op_id, action,
  entity_id) tuple).
- `view_row_id`: the id of the view row inserted by the dispatcher, or
  `nothing` for non-insert dispatcher branches (DELETE, no-op, etc.) and on
  idempotent retry (where the dispatcher's prior insert is not re-derivable
  from the unique-index lookup). Callers that need the inserted row id
  (e.g. POST /peaks) use this directly instead of re-querying.

`payload` is any JSON-serializable Dict / NamedTuple / nothing. If nothing,
the event is recorded but no view update fires (use sparingly — most actions
should carry a payload).

This default method opens its own `SQLite.transaction`, then delegates to
`apply_event!(::InTransaction, db, req; ...)`. The SSE broadcast fires AFTER
the transaction commits so subscribers can never see uncommitted state.
Callers that need to participate in an outer transaction (e.g. routes
wrapped in `with_idempotency`) should use the `InTransaction` variant.
"""
function apply_event!(db::SQLite.DB, req;
                      kind::String,
                      entity_type::String,
                      entity_id::Integer,
                      payload = nothing,
                      undoes_event_id::Union{Int,Nothing} = nothing,
                      defer_broadcast::Bool = false,
                      post_state::Union{Dict, Nothing} = nothing)
    # Run the durable write inside a tx via the InTransaction variant, but
    # always defer the broadcast there — broadcast must wait until AFTER the
    # tx commits so subscribers can't see uncommitted state.
    result = SQLite.transaction(db) do
        apply_event!(InTransaction(), db, req;
                     kind = kind, entity_type = entity_type, entity_id = entity_id,
                     payload = payload, undoes_event_id = undoes_event_id,
                     defer_broadcast = true,
                     post_state = post_state)
    end

    # Now committed. Fire the broadcast unless the outer caller asked to defer
    # it themselves (e.g. coalesced batch broadcast).
    if !defer_broadcast
        _maybe_broadcast_event!(db, req, result, kind, entity_type, entity_id,
                                payload, post_state)
    end

    return (event_id = result.event_id, view_row_id = result.view_row_id)
end

"""
    apply_event!(::InTransaction, db, req; kwargs...)

In-transaction variant: the caller has already opened a `SQLite.transaction`.
Performs the INSERT into `user_actions` plus dispatcher view-update inside
that transaction. Does NOT broadcast (that happens after the outer tx
commits — caller is responsible).

Idempotent at the DB layer: when `client_op_id` is set and the partial
unique index on `(client_op_id, action, entity_id)` rejects the INSERT,
this looks up the existing event row and returns its `event_id`.

**On idempotent retry (UNIQUE constraint trip):** the returned `view_row_id`
is `nothing` because the dispatcher's prior INSERT isn't re-derivable from
the event row alone, and the dispatcher is NOT re-run (the prior application
already moved the views). Routes whose response shape depends on
`view_row_id` must be wrapped in `with_idempotency` so the cached HTTP
response — which carries the original `view_row_id` — is replayed on retry.
The current default-method-only callers (routes_peaks, routes_analysis)
don't yet wrap in `with_idempotency`, so a same-`X-Client-Op-Id` retry today
would 500. M2 routes will adopt `with_idempotency` and resolve this.

Returns a richer NamedTuple than the public default method — includes the
fields needed for a deferred post-commit broadcast.
"""
function apply_event!(::InTransaction, db::SQLite.DB, req;
                      kind::String,
                      entity_type::String,
                      entity_id::Integer,
                      payload = nothing,
                      undoes_event_id::Union{Int,Nothing} = nothing,
                      defer_broadcast::Bool = true,  # accepted for kw-symmetry; ignored here
                      post_state::Union{Dict, Nothing} = nothing)
    username     = get_username(req)
    client_id    = get_client_id(req)
    client_op_id = get_client_op_id(req)
    user_id      = username === nothing ? nothing : get_or_create_user!(db, username)
    payload_json = payload === nothing ? nothing : JSON3.write(payload)

    event_id::Int = 0
    view_row_id::Union{Int, Nothing} = nothing
    fresh_insert = true

    # Narrow scope: the try/catch wraps ONLY the user_actions INSERT (and
    # lastrowid extraction). The dispatcher runs *outside* the catch's reach
    # so a future view-INSERT that happens to trip its own UNIQUE constraint
    # can't be misclassified as an idempotent retry of the event-log INSERT.
    try
        res = DBInterface.execute(db,
            """INSERT INTO user_actions
               (user_id, action, entity_type, entity_id, payload, undoes_event_id, client_id, client_op_id)
               VALUES (?, ?, ?, ?, ?, ?, ?, ?)""",
            [user_id, kind, entity_type, Int(entity_id), payload_json, undoes_event_id, client_id, client_op_id])
        event_id = Int(DBInterface.lastrowid(res))
    catch err
        # Idempotent retry: the partial unique index on
        # (client_op_id, action, entity_id) rejected the INSERT because a
        # prior call already applied this op. SELECT the existing event_id
        # and return it. (SQLite raises a generic SQLiteException; match
        # by message text since there's no stable error code surface.)
        if client_op_id !== nothing && occursin("UNIQUE constraint failed", sprint(showerror, err))
            existing = Tables.rowtable(DBInterface.execute(db,
                """SELECT id FROM user_actions
                   WHERE client_op_id = ? AND action = ? AND entity_id = ?
                   LIMIT 1""",
                [client_op_id, kind, Int(entity_id)]))
            if !isempty(existing)
                event_id = Int(existing[1].id)
                fresh_insert = false
            else
                rethrow()
            end
        else
            rethrow()
        end
    end

    # Run the dispatcher only on a fresh INSERT — on retry the prior
    # application already moved the views, and re-running would double-apply.
    # Canonicalize payload before dispatch — round-trip through JSON3 so the
    # live dispatcher and rebuild_views_from_log! see exactly the same shape
    # (JSON3.Object supports both .field and [:field] access, eliminating
    # the Symbol-key vs String-key footgun).
    if fresh_insert && payload_json !== nothing
        payload_canonical = JSON3.read(payload_json)
        view_row_id = update_view_for_event!(db, kind, entity_id, payload_canonical, event_id)
    end

    return (event_id    = event_id,
            view_row_id = view_row_id,
            user_id     = user_id,
            client_id   = client_id,
            client_op_id = client_op_id,
            payload_json = payload_json)
end

"""
    _maybe_broadcast_event!(db, req, result, kind, entity_type, entity_id, payload, post_state)

Internal: called by the default `apply_event!` after the durable transaction
commits. Skips broadcast for analyze_run no-ops (the M0.4 suppression rule)
and tolerates a missing or failing `broadcast_event!`.
"""
function _maybe_broadcast_event!(db, req, result, kind, entity_type, entity_id, payload, post_state)
    isdefined(@__MODULE__, :broadcast_event!) || return nothing

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
    suppress && return nothing

    try
        broadcast_event!(result.event_id, kind, entity_type, Int(entity_id),
                         result.user_id, result.client_id, result.client_op_id,
                         result.payload_json;
                         post_state = post_state)
    catch err
        @warn "broadcast_event! failed (event still durable in user_actions)" exception=err
    end
    nothing
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

    # M2.1 trivial-route migrations: routes write to view tables directly,
    # so the dispatcher is a no-op for these kinds. Branches exist for
    # exhaustiveness so the rebuild_views_from_log! property test treats
    # them as known kinds rather than silently falling through.
    kind == "update_sample" && return nothing
    kind == "add_tag" && return nothing
    kind == "remove_tag" && return nothing
    kind == "post_message" && return nothing
    kind == "set_exposure_status" && return nothing
    kind == "select_exposure" && return nothing

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
