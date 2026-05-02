using JSON3, SQLite, DBInterface, HTTP

"""
    apply_event!(db, req; kind, entity_type, entity_id, payload, undoes_event_id=nothing)
      -> event_id::Int

Atomic event-append + view-update. The log and the views must move together
or neither moves. Returns the newly-inserted event id.

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
            update_view_for_event!(db, kind, entity_id, payload_canonical, eid)
        end
        eid
    end

    # Best-effort SSE broadcast — fires AFTER the transaction commits, so a
    # subscriber never sees an event that was rolled back. If the process
    # dies between commit and broadcast, the event is durable in user_actions
    # but the frame is lost; clients reconcile on reconnect via TanStack
    # Query refetch (see R5a). isdefined check keeps the helper a soft
    # dependency: R4 lands before R5a, and apply_event! must work without
    # broadcast wired up yet.
    if isdefined(@__MODULE__, :broadcast_event!)
        try
            broadcast_event!(event_id, kind, entity_type, Int(entity_id), user_id, payload_json)
        catch err
            @warn "broadcast_event! failed (event still durable in user_actions)" exception=err
        end
    end
    event_id
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
    # Implementations added per route in R4.2. Initial scaffolding only:
    kind == "noop_test" && return
    # default: no view update
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
