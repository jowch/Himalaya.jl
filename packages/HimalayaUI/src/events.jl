using JSON3, SQLite, DBInterface, HTTP, Tables

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
                      undoes_event_id::Union{Int,Nothing} = nothing)
    username = get_username(req)
    user_id  = username === nothing ? nothing : get_or_create_user!(db, username)
    payload_json = payload === nothing ? nothing : JSON3.write(payload)

    view_row_id_ref = Ref{Union{Int, Nothing}}(nothing)
    event_id = SQLite.transaction(db) do
        res = DBInterface.execute(db,
            """INSERT INTO user_actions
               (user_id, action, entity_type, entity_id, payload, undoes_event_id)
               VALUES (?, ?, ?, ?, ?, ?)""",
            [user_id, kind, entity_type, Int(entity_id), payload_json, undoes_event_id])
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
    if isdefined(@__MODULE__, :broadcast_event!)
        try
            broadcast_event!(event_id, kind, entity_type, Int(entity_id), user_id, payload_json)
        catch err
            @warn "broadcast_event! failed (event still durable in user_actions)" exception=err
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
    broadcast_event!(event_id, kind, entity_type, entity_id, user_id, payload_json)

Format a single SSE frame and enqueue it onto every subscriber's pending
channel. Closed channels (disconnected clients) and full channels (slow
subscribers) are pruned — the client will reconnect via EventSource
auto-reconnect and refetch via TanStack Query.

Best-effort: this fires AFTER apply_event!'s transaction commits, so a
subscriber never sees an event that was rolled back. If the process dies
between commit and broadcast, the event is durable in user_actions but the
frame is lost; clients reconcile on reconnect via TanStack Query refetch.

SSE_SUBSCRIBERS and SSE_LOCK live in server.jl but are visible here because
both files are included into the same HimalayaUI module.
"""
function broadcast_event!(event_id::Integer, kind::String, entity_type::String,
                          entity_id::Integer, user_id::Union{Integer, Nothing},
                          payload_json::Union{String, Nothing})
    actor = user_id === nothing ? nothing : lookup_username(current_db(), user_id)
    msg = JSON3.write(Dict(
        :id          => Int(event_id),
        :kind        => kind,
        :entity_type => entity_type,
        :entity_id   => Int(entity_id),
        :actor       => actor,
        :payload     => payload_json === nothing ? nothing : JSON3.read(payload_json),
    ))
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
