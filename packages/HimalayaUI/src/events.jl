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
    # _DB_WRITE_LOCK (#122) serializes the tx against any other singleton
    # writer; reentrant so a route already holding the lock (e.g. through
    # `with_idempotency`) re-enters cleanly.
    result = lock(_DB_WRITE_LOCK) do
        SQLite.transaction(db) do
            apply_event!(InTransaction(), db, req;
                         kind = kind, entity_type = entity_type, entity_id = entity_id,
                         payload = payload, undoes_event_id = undoes_event_id,
                         defer_broadcast = true,
                         post_state = post_state)
        end
    end

    # Now committed. Fire the broadcast unless the outer caller asked to defer
    # it themselves (e.g. coalesced batch broadcast).
    if !defer_broadcast
        _maybe_broadcast_event!(db, result, kind, entity_type, entity_id,
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
    _maybe_broadcast_event!(db, result, kind, entity_type, entity_id, payload, post_state)

Internal: called by the default `apply_event!` after the durable transaction
commits. Skips broadcast for analyze_run no-ops (the M0.4 suppression rule)
and tolerates a missing or failing `broadcast_event!`.
"""
function _maybe_broadcast_event!(db, result, kind, entity_type, entity_id, payload, post_state)
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
    _enqueue_post_commit_broadcast!(args...)

Queue a `broadcast_event!` invocation to fire AFTER the current
`with_idempotency` transaction commits. Stored in task-local storage so each
request handler gets its own queue. Cleared without firing on rollback.

Used by M2.2 peak/curation routes that need to emit a single enriched SSE
frame (carrying `post_state`) only after the outer with_idempotency tx
commits — broadcasting earlier would let subscribers see state that may roll
back.
"""
const POST_COMMIT_BROADCAST_KEY = :himalaya_post_commit_broadcasts

function _enqueue_post_commit_broadcast!(args...; kwargs...)
    queue = get!(task_local_storage(), POST_COMMIT_BROADCAST_KEY) do
        Vector{Tuple{Tuple, Base.Pairs}}()
    end
    push!(queue, (args, kwargs))
    nothing
end

"""
    _enqueue_broadcast_from_result!(result, kind, entity_type, entity_id;
                                    post_state = nothing)

Convenience: queue a post-commit SSE broadcast from the NamedTuple returned
by `apply_event!(InTransaction(), ...)`. Routes wrapped in `with_idempotency`
that call the InTransaction variant must explicitly enqueue their broadcast
(the InTransaction variant does NOT broadcast — see `apply_event!` docstring).
This helper centralizes the field-by-field unpack so callers don't drift.

Reuses `result.payload_json` (the canonical serialization frozen by
`apply_event!` after JSON3 round-trip) — no re-serialization needed.
"""
function _enqueue_broadcast_from_result!(result, kind::String,
                                         entity_type::String, entity_id::Integer;
                                         post_state::Union{Dict, Nothing} = nothing)
    _enqueue_post_commit_broadcast!(
        Int(result.event_id), kind, entity_type, Int(entity_id),
        result.user_id, result.client_id, result.client_op_id,
        result.payload_json;
        post_state = post_state)
end

"""
    _flush_post_commit_broadcasts!()

Fire every queued broadcast for the current task and clear the queue. Called
by `with_idempotency` after its `SQLite.transaction` commits successfully.
Failures in individual broadcasts are logged but do not abort the flush —
each frame is best-effort (matches `_maybe_broadcast_event!` semantics).
"""
function _flush_post_commit_broadcasts!()
    queue = get(task_local_storage(), POST_COMMIT_BROADCAST_KEY, nothing)
    queue === nothing && return nothing
    for (args, kwargs) in queue
        try
            broadcast_event!(args...; kwargs...)
        catch err
            @warn "post-commit broadcast failed (event still durable in user_actions)" exception=err
        end
    end
    delete!(task_local_storage(), POST_COMMIT_BROADCAST_KEY)
    nothing
end

"""
    _clear_post_commit_broadcasts!()

Discard any queued post-commit broadcasts for the current task without
firing them. Called by `with_idempotency` when the tx body throws (rollback)
so subscribers never see events whose underlying writes didn't durably
commit.
"""
function _clear_post_commit_broadcasts!()
    delete!(task_local_storage(), POST_COMMIT_BROADCAST_KEY)
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

    # peak_removed: the route handler deletes the peak_curations(kind='add') row
    # directly (it has the integer id from the URL), so the dispatcher is a
    # no-op. Branch exists for exhaustiveness — rebuild_views_from_log! treats
    # it as a known kind rather than silently falling through.
    kind == "peak_removed" && return nothing

    # index_confirmed / index_unconfirmed: RETIRED legacy group-membership kinds
    # (plotting redesign D-10). The /groups routes that emitted them are gone and
    # confirmed_index now sources from the durable assignment, not index_groups.
    # These remain as explicit no-op GUARDS — never delete them: historical events
    # still live in user_actions, and this keeps them recognized KNOWN kinds (the
    # branch is for exhaustiveness/intent, matching the other retired-kind guards;
    # the dispatcher's default is itself a silent no-op, so this isn't strictly
    # throw-prevention). No-op is replay-consistent: live apply and replay both
    # write nothing, so the per-event round-trip property holds.
    #
    # RECOVERY CAVEAT (see migrate_assignments!): the *state* a pre-Plan-A
    # confirmation produced is NOT reproducible from the log alone — those
    # exposures have only index_confirmed events (now no-ops) and no paired
    # assignment_add, so their durable assignment exists solely because the
    # sentinel-gated migrate_assignments! backfilled it from index_groups. A
    # from-empty rebuild_views_from_log! (drop assignment tables + re-fold) will
    # NOT reconstruct them; it must be preceded by clearing the assignments_v1
    # sentinel and re-running migrate_assignments!. Post-Plan-A confirmations
    # ride assignment_add and DO round-trip through the log.
    (kind == "index_confirmed" || kind == "index_unconfirmed") && return nothing

    # Plotting redesign Plan A: durable per-exposure assignment kinds. Sole
    # writer to assignments/assignment_members. entity_id is the exposure id.
    # Replay-idempotent (UPSERT / INSERT OR IGNORE / DELETE) so the fold from
    # an empty view reproduces live state.
    if kind == "assignment_add"
        DBInterface.execute(db,
            """INSERT INTO assignments (exposure_id, state) VALUES (?, 'indexed')
               ON CONFLICT(exposure_id) DO UPDATE SET state = 'indexed'""",
            [Int(entity_id)])
        DBInterface.execute(db,
            """INSERT OR IGNORE INTO assignment_members (exposure_id, index_id)
               VALUES (?, ?)""",
            [Int(entity_id), Int(payload.index_id)])
        return nothing
    end

    if kind == "assignment_remove"
        DBInterface.execute(db,
            "DELETE FROM assignment_members WHERE exposure_id = ? AND index_id = ?",
            [Int(entity_id), Int(payload.index_id)])
        return nothing
    end

    if kind == "assignment_set_state"
        state = String(payload.state)
        DBInterface.execute(db,
            """INSERT INTO assignments (exposure_id, state) VALUES (?, ?)
               ON CONFLICT(exposure_id) DO UPDATE SET state = excluded.state""",
            [Int(entity_id), state])
        if state != "indexed"
            DBInterface.execute(db,
                "DELETE FROM assignment_members WHERE exposure_id = ?", [Int(entity_id)])
        end
        return nothing
    end

    # M2.1 trivial-route migrations: routes write to view tables directly,
    # so the dispatcher is a no-op for these kinds. Branches exist for
    # exhaustiveness so the rebuild_views_from_log! property test treats
    # them as known kinds rather than silently falling through.
    kind == "update_sample" && return nothing
    kind == "add_tag" && return nothing
    kind == "remove_tag" && return nothing
    kind == "edit_tag" && return nothing
    kind == "post_message" && return nothing
    kind == "set_exposure_status" && return nothing
    kind == "select_exposure" && return nothing

    # Speculative index lifecycle: route handlers insert/delete the indices
    # row directly (the speculative create/delete paths in routes_analysis.jl
    # need the auto-generated id immediately to populate the response body).
    # Dispatcher branches exist for exhaustiveness only.
    kind == "speculative_created" && return nothing
    kind == "speculative_deleted" && return nothing

    # analyze_run: pure observability event — no view writes. The synchronous
    # reanalyze inside curation routes mutates indices/auto_peaks via
    # persist_analysis!, not via this dispatcher. Branch exists so the
    # rebuild_views_from_log! property test treats it as a known kind.
    kind == "analyze_run" && return nothing

    # Compare page (Plan §Phase 1, Task 1.2): three view-producing kinds.
    # (compare retired; see git history)
    if kind == "comparison_created"
        return _update_view_for_comparison_created!(db, entity_id, payload, event_id)
    end
    if kind == "comparison_submitted"
        return _update_view_for_comparison_submitted!(db, entity_id, payload, event_id)
    end
    if kind == "comparison_deleted"
        DBInterface.execute(db, "DELETE FROM comparisons WHERE id = ?", [Int(entity_id)])
        # The schema's ON DELETE CASCADE clauses drop comparison_members and
        # comparison_messages automatically.
        return nothing
    end

    # Per-user pin/unpin (Plan §Phase 13 follow-up): pin/unpin events live on
    # `entity_type='user'` so the durable history is queryable per-user via
    # the existing `user_actions` indexes. The view-table is `comparison_pins`
    # (composite PK: user_id, comparison_id). Payload carries `comparison_id`
    # plus `pinned_at` for the pin variant. Dispatcher derives user_id by
    # joining `user_actions WHERE id = event_id` — same pattern peak_added
    # uses to extract created_by.
    if kind == "comparison_pinned"
        DBInterface.execute(db,
            """INSERT OR REPLACE INTO comparison_pins (user_id, comparison_id, pinned_at)
               VALUES ((SELECT user_id FROM user_actions WHERE id = ?),
                       ?, CURRENT_TIMESTAMP)""",
            [event_id, Int(payload.comparison_id)])
        return nothing
    end

    if kind == "comparison_unpinned"
        DBInterface.execute(db,
            """DELETE FROM comparison_pins
               WHERE user_id = (SELECT user_id FROM user_actions WHERE id = ?)
                 AND comparison_id = ?""",
            [event_id, Int(payload.comparison_id)])
        return nothing
    end

    # Series event kinds (#166 / I2.3). The view-producing branches are
    # pure-replace: parent upsert + child DELETE-by-series_id + INSERT-all from
    # a full-snapshot payload — explicitly NOT comparison_submitted-shaped, so
    # every fold is idempotent and order-independent (folds from an empty view).
    if kind == "series_created"
        return _update_view_for_series_created!(db, entity_id, payload, event_id)
    end
    if kind == "series_recipe_updated"
        return _update_view_for_series_recipe_updated!(db, entity_id, payload, event_id)
    end
    if kind == "series_deleted"
        DBInterface.execute(db, "DELETE FROM series WHERE id = ?", [Int(entity_id)])
        # The schema's four ON DELETE CASCADE clauses drop series_members,
        # series_samples, series_messages and series_pins automatically.
        return nothing
    end
    if kind == "series_plate_committed"
        return _update_view_for_series_plate_committed!(db, entity_id, payload, event_id)
    end

    # Per-user series pins (#168 / I2.5). Stored with entity_type='user',
    # entity_id=user_id (the comparison_pinned precedent) — the affected series
    # rides in the payload as `series_id`. Five-layer: no post_state, no
    # mutator. The dispatcher derives user_id by joining user_actions on
    # event_id, exactly as the comparison pin branches do.
    if kind == "series_pinned"
        # Guard the pin on series existence: a rebuild that folds this event
        # after the series was deleted would otherwise FK-throw and abort the
        # whole rebuild. INSERT … SELECT … WHERE EXISTS makes the fold a no-op
        # in that case. (Deliberately stricter than the unguarded
        # comparison_pinned branch, which carries the same latent issue.)
        DBInterface.execute(db,
            """INSERT OR REPLACE INTO series_pins (user_id, series_id, pinned_at)
               SELECT (SELECT user_id FROM user_actions WHERE id = ?),
                      ?, CURRENT_TIMESTAMP
               WHERE EXISTS (SELECT 1 FROM series WHERE id = ?)""",
            [event_id, Int(payload.series_id), Int(payload.series_id)])
        return nothing
    end
    if kind == "series_unpinned"
        DBInterface.execute(db,
            """DELETE FROM series_pins
               WHERE user_id = (SELECT user_id FROM user_actions WHERE id = ?)
                 AND series_id = ?""",
            [event_id, Int(payload.series_id)])
        return nothing
    end

    # Scaffolding / legacy:
    kind == "noop_test" && return nothing
    # default: no view update
    nothing
end

"""
    _update_view_for_comparison_created!(db, entity_id, payload, event_id)

Insert the `comparisons` row and any payload-supplied members. `entity_id`
is the comparison id (always pre-allocated by the route handler — the
event-driven dispatcher cannot assign autoincrement ids reliably across
idempotent retry, so the route mints the id with an explicit INSERT and
passes it as `entity_id`). Returns the comparison id for `view_row_id`.

The route is responsible for the route-side `INSERT INTO comparisons`
that mints the id; this dispatcher branch instead **upserts** so a fresh
INSERT plus a replay land on the same row. In practice the route+dispatcher
sequence is:
  1. route INSERTs a placeholder comparison row, captures its autoincrement id.
  2. apply_event!(InTransaction, …, entity_id=that_id).
  3. dispatcher updates title/description/content_hash on the placeholder
     and inserts members.
"""
function _update_view_for_comparison_created!(db, entity_id, payload, event_id)
    user_id = user_id_for_event(db, event_id)
    now_str = comparison_now_iso()
    title = String(payload.title)
    description = haskey(payload, :description) && payload.description !== nothing ?
                  String(payload.description) : nothing
    forked_from_id = haskey(payload, :forked_from_id) && payload.forked_from_id !== nothing ?
                     Int(payload.forked_from_id) : nothing
    forked_at_hash = haskey(payload, :forked_at_hash) && payload.forked_at_hash !== nothing ?
                     String(payload.forked_at_hash) : nothing
    # View-choice columns (spec §6.4). Bare write (no COALESCE) — clients can
    # null these to reset to per-tab defaults. `String()` coercion normalises
    # JSON3 string subtypes before write.
    vgm  = haskey(payload, :view_grouping_mode) && payload.view_grouping_mode !== nothing ?
           String(payload.view_grouping_mode) : nothing
    vspt = haskey(payload, :view_show_peak_ticks)  ? payload.view_show_peak_ticks  : nothing
    vspl = haskey(payload, :view_show_peak_labels) ? payload.view_show_peak_labels : nothing

    # Insert/upsert the comparisons row at this id. The route may have minted
    # the row already (to capture the AUTOINCREMENT id); if so this UPDATE
    # supersedes that placeholder. Otherwise INSERT — the route can choose to
    # rely on dispatcher-driven creation.
    existing = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM comparisons WHERE id = ?", [Int(entity_id)]))
    if isempty(existing)
        # Use INSERT with explicit id; sqlite's AUTOINCREMENT counter will
        # advance past entity_id automatically. (Suitable for replay paths
        # where the original row was deleted but the event is being re-folded.)
        # `content_hash` lands NULL here and is overwritten unconditionally by
        # compute_content_hash below; the column is nullable post-#67.
        DBInterface.execute(db,
            """INSERT INTO comparisons
               (id, title, description, content_hash, created_by, created_at,
                updated_at, forked_from_id, forked_at_hash,
                view_grouping_mode, view_show_peak_ticks, view_show_peak_labels)
               VALUES (?, ?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, ?)""",
            [Int(entity_id), title, description, user_id,
             now_str, now_str, forked_from_id, forked_at_hash,
             vgm, vspt, vspl])
    else
        # The route's mint-the-id INSERT seeds NULL placeholders (#67), so
        # plain COALESCE stamps `created_at` on first fold. A real timestamp
        # from a prior fold survives untouched (replay path). Pre-#67 this
        # was `COALESCE(NULLIF(created_at, ''), ?)` — the NULLIF wrapper was
        # the #66 patch for the placeholder being `''` rather than NULL.
        DBInterface.execute(db,
            """UPDATE comparisons
               SET title = ?, description = ?, created_by = ?,
                   created_at = COALESCE(created_at, ?), updated_at = ?,
                   forked_from_id = ?, forked_at_hash = ?,
                   view_grouping_mode = ?, view_show_peak_ticks = ?,
                   view_show_peak_labels = ?
               WHERE id = ?""",
            [title, description, user_id, now_str, now_str,
             forked_from_id, forked_at_hash,
             vgm, vspt, vspl, Int(entity_id)])
    end

    # Insert members (all members in a comparison_created payload are NEW —
    # ids in the payload are ignored; the dispatcher mints fresh PKs).
    members = haskey(payload, :members) ? payload.members : []
    for m in members
        _insert_comparison_member!(db, Int(entity_id), m, user_id, now_str)
    end

    # Recompute the canonical content_hash from the just-written state and
    # bump the comparisons row.
    new_hash = compute_content_hash(db, Int(entity_id))
    DBInterface.execute(db,
        "UPDATE comparisons SET content_hash = ? WHERE id = ?",
        [new_hash, Int(entity_id)])

    return Int(entity_id)
end

"""
    _update_view_for_comparison_submitted!(db, entity_id, payload, event_id)

Diff the payload's member list against the current `comparison_members`
rows for this comparison. Three dispositions:
- **DELETE** rows present in the DB but absent from the payload.
- **UPDATE** rows present in both (`payload.id isa Int`); the snapshot is
  written unconditionally — see spec §Submission diff: "Every UPDATE
  recomputes the snapshot from the payload."
- **INSERT** rows in the payload with `id === nothing`.

Errors on a zero-row UPDATE (payload references a member id that doesn't
belong to this comparison). After the diff, recomputes `content_hash` and
bumps `updated_at`; may also update title/description.
"""
function _update_view_for_comparison_submitted!(db, entity_id, payload, event_id)
    members = haskey(payload, :members) ? payload.members : []
    payload_existing = Dict{Int, Any}()
    payload_new = Any[]
    for m in members
        # `m.id` arrives as JSON3.Object's accessor; nothing/null → INSERT,
        # integer → UPDATE. JSON3 decodes JSON null as `nothing` in this path.
        mid = haskey(m, :id) ? m.id : nothing
        if mid === nothing
            push!(payload_new, m)
        else
            payload_existing[Int(mid)] = m
        end
    end

    current_ids = member_ids_for_comparison(db, entity_id)

    # DELETE: in DB, not in payload.
    for id in setdiff(current_ids, keys(payload_existing))
        DBInterface.execute(db,
            "DELETE FROM comparison_members WHERE id = ? AND comparison_id = ?",
            [Int(id), Int(entity_id)])
    end

    # UPDATE: in both.
    for (id, m) in payload_existing
        DBInterface.execute(db,
            """UPDATE comparison_members
               SET exposure_id = ?, display_order = ?, band_height = ?,
                   y_offset = ?, normalization = ?, color_override = ?,
                   label_override = ?, q_window_min = ?, q_window_max = ?,
                   peak_display = ?, snapshot = ?
               WHERE id = ? AND comparison_id = ?""",
            [_member_field(m, :exposure_id),
             Int(_member_field(m, :display_order)),
             Float64(_member_field(m, :band_height; default=1.0)),
             Float64(_member_field(m, :y_offset; default=0.0)),
             String(_member_field(m, :normalization; default="none")),
             _member_field(m, :color_override),
             _member_field(m, :label_override),
             _member_field(m, :q_window_min),
             _member_field(m, :q_window_max),
             _member_json(m, :peak_display),
             _member_json(m, :snapshot)::String,  # NOT NULL — JSON3.write of payload
             Int(id), Int(entity_id)])
        # SQLite's changes() returns the row-count of the most recent statement.
        # Used here to reject UPDATEs whose `id` doesn't belong to this
        # comparison — payload references a stale or cross-comparison id.
        n_changed = Int(first(Tables.rowtable(DBInterface.execute(db,
            "SELECT changes() AS n"))).n)
        n_changed == 0 &&
            error("comparison_submitted: member id=$id not found for comparison $entity_id")
    end

    # INSERT: new members (id === nothing in payload).
    user_id = user_id_for_event(db, event_id)
    now_str = comparison_now_iso()
    for m in payload_new
        _insert_comparison_member!(db, Int(entity_id), m, user_id, now_str)
    end

    # Title / description / view-choice update — must happen BEFORE the hash
    # recompute, otherwise compute_content_hash sees the pre-update title and
    # the content_hash never moves on a rename-only submit (caught by the
    # routes' "rename + same membership → hash changes" regression test).
    # `String()` coercion normalises JSON3 string subtypes before write.
    # title / description use COALESCE (omit keeps current); the three
    # view-choice columns use a bare write so clients can null them to reset
    # to per-tab defaults (spec §6.4). The content_hash does NOT fold the
    # view-choice columns — no change to hash.jl.
    title_val = haskey(payload, :title) && payload.title !== nothing ?
                String(payload.title) : nothing
    desc_val  = haskey(payload, :description) && payload.description !== nothing ?
                String(payload.description) : nothing
    vgm  = haskey(payload, :view_grouping_mode) && payload.view_grouping_mode !== nothing ?
           String(payload.view_grouping_mode) : nothing
    vspt = haskey(payload, :view_show_peak_ticks)  ? payload.view_show_peak_ticks  : nothing
    vspl = haskey(payload, :view_show_peak_labels) ? payload.view_show_peak_labels : nothing
    DBInterface.execute(db,
        """UPDATE comparisons SET
              title = COALESCE(?, title),
              description = COALESCE(?, description),
              view_grouping_mode = ?,
              view_show_peak_ticks = ?,
              view_show_peak_labels = ?
           WHERE id = ?""",
        [title_val, desc_val, vgm, vspt, vspl, Int(entity_id)])
    new_hash = compute_content_hash(db, Int(entity_id))
    DBInterface.execute(db,
        "UPDATE comparisons SET content_hash = ?, updated_at = ? WHERE id = ?",
        [new_hash, now_str, Int(entity_id)])
    return nothing
end

# Helper: read a payload member field tolerating JSON3.Object's both-keys
# (`m.field` and `m[:field]`) interface and a missing key. Returns `nothing`
# for JSON null or absent key. Numeric defaults are applied via the kwarg.
function _member_field(m, key::Symbol; default=nothing)
    haskey(m, key) || return default
    v = getproperty(m, key)
    v === nothing ? default : v
end

# Helper: serialize a (possibly nil) payload sub-object back to JSON for
# the snapshot/peak_display TEXT columns. Returns `nothing` (→ NULL) when
# the field is absent or null. The payload value is already a JSON3.Object
# or nested Dict — re-serializing is the canonical way to land it in the
# TEXT column.
function _member_json(m, key::Symbol)
    haskey(m, key) || return nothing
    v = getproperty(m, key)
    v === nothing && return nothing
    JSON3.write(v)
end

# Helper: INSERT one comparison_members row, validating that the snapshot
# field is present (NOT NULL CHECK in the schema; the route is responsible
# for constructing it on the client side).
function _insert_comparison_member!(db, comparison_id, m, user_id, now_str)
    snap = _member_json(m, :snapshot)
    snap === nothing &&
        error("comparison member missing required `snapshot` field")
    DBInterface.execute(db,
        """INSERT INTO comparison_members
             (comparison_id, exposure_id, display_order, band_height, y_offset,
              normalization, color_override, label_override,
              q_window_min, q_window_max, peak_display, snapshot,
              created_by, created_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        [Int(comparison_id),
         _member_field(m, :exposure_id),
         Int(_member_field(m, :display_order)),
         Float64(_member_field(m, :band_height; default=1.0)),
         Float64(_member_field(m, :y_offset; default=0.0)),
         String(_member_field(m, :normalization; default="none")),
         _member_field(m, :color_override),
         _member_field(m, :label_override),
         _member_field(m, :q_window_min),
         _member_field(m, :q_window_max),
         _member_json(m, :peak_display),
         snap,
         user_id, now_str])
    nothing
end

# Helper: INSERT one series_samples row from a recipe payload entry. The
# payload entry is a JSON3.Object from `_series_sample_payload` —
# `{sample_id, position, pinned, excluded}`. `pinned`/`excluded` are coerced
# to 0/1 for the CHECK (… IN (0,1)) columns.
function _insert_series_sample!(db, series_id, s)
    pinned   = (haskey(s, :pinned)   && s.pinned   == true) ? 1 : 0
    excluded = (haskey(s, :excluded) && s.excluded == true) ? 1 : 0
    DBInterface.execute(db,
        """INSERT INTO series_samples
             (series_id, sample_id, position, pinned, excluded)
           VALUES (?, ?, ?, ?, ?)""",
        [Int(series_id), Int(s.sample_id), Int(s.position), pinned, excluded])
    nothing
end

"""
    _update_view_for_series_created!(db, entity_id, payload, event_id)

`series_created` dispatcher (#166). Upserts the `series` row at `entity_id`
(the route mints it with `INSERT INTO series DEFAULT VALUES` to capture the
AUTOINCREMENT id; a plain INSERT would collide on the live path — so this
SELECTs and UPDATEs an existing row, else INSERTs with an explicit id for the
replay-from-empty path). Sets `state='draft'`; `content_hash` stays NULL (a
draft has no committed plate — master plan §5.1). Then pure-replaces
`series_samples` from the full payload snapshot and resolves the plate
(`series_members`) from that recipe via `_resolve_series_plate!` (decision
2026-06: a created draft lands renderable — `content_hash` stays NULL since the
plate is not yet committed).
"""
function _update_view_for_series_created!(db, entity_id, payload, event_id)
    sid     = Int(entity_id)
    user_id = user_id_for_event(db, event_id)
    now_str = comparison_now_iso()

    title          = haskey(payload, :title) && payload.title !== nothing ?
                     String(payload.title) : nothing
    description    = haskey(payload, :description) && payload.description !== nothing ?
                     String(payload.description) : nothing
    forked_from_id = haskey(payload, :forked_from_id) && payload.forked_from_id !== nothing ?
                     Int(payload.forked_from_id) : nothing
    forked_at_hash = haskey(payload, :forked_at_hash) && payload.forked_at_hash !== nothing ?
                     String(payload.forked_at_hash) : nothing
    ordering_var   = haskey(payload, :ordering_variable) && payload.ordering_variable !== nothing ?
                     String(payload.ordering_variable) : nothing
    order_rule     = haskey(payload, :order_rule) && payload.order_rule !== nothing ?
                     String(payload.order_rule) : "manual"
    vgm  = haskey(payload, :view_grouping_mode) && payload.view_grouping_mode !== nothing ?
           String(payload.view_grouping_mode) : nothing
    vspt = haskey(payload, :view_show_peak_ticks)  ? payload.view_show_peak_ticks  : nothing
    vspl = haskey(payload, :view_show_peak_labels) ? payload.view_show_peak_labels : nothing

    existing = Tables.rowtable(DBInterface.execute(db,
        "SELECT id FROM series WHERE id = ?", [sid]))
    if isempty(existing)
        # Replay path: the original row was deleted. INSERT with an explicit id;
        # SQLite's AUTOINCREMENT counter advances past it automatically.
        DBInterface.execute(db,
            """INSERT INTO series
               (id, title, description, content_hash, created_by, created_at,
                updated_at, forked_from_id, forked_at_hash,
                view_grouping_mode, view_show_peak_ticks, view_show_peak_labels,
                ordering_variable, order_rule, state)
               VALUES (?, ?, ?, NULL, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, 'draft')""",
            [sid, title, description, user_id, now_str, now_str,
             forked_from_id, forked_at_hash, vgm, vspt, vspl,
             ordering_var, order_rule])
    else
        # Live path: the route's `INSERT … DEFAULT VALUES` placeholder exists.
        # COALESCE stamps created_at on first fold; a prior fold's value
        # survives. `content_hash = NULL` is set explicitly so the draft
        # invariant holds even on a defensive re-fold onto an already-committed
        # row.
        DBInterface.execute(db,
            """UPDATE series
               SET title = ?, description = ?, content_hash = NULL,
                   created_by = ?,
                   created_at = COALESCE(created_at, ?), updated_at = ?,
                   forked_from_id = ?, forked_at_hash = ?,
                   view_grouping_mode = ?, view_show_peak_ticks = ?,
                   view_show_peak_labels = ?,
                   ordering_variable = ?, order_rule = ?, state = 'draft'
               WHERE id = ?""",
            [title, description, user_id, now_str, now_str,
             forked_from_id, forked_at_hash, vgm, vspt, vspl,
             ordering_var, order_rule, sid])
    end

    # Pure-replace the recipe rows from the full-snapshot payload.
    DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
    samples = haskey(payload, :samples) ? payload.samples : []
    for s in samples
        _insert_series_sample!(db, sid, s)
    end

    # Resolve the plate (series_members) from the just-written recipe so the
    # created draft renders its waterfall immediately (decision 2026-06: a
    # created series lands in the builder already showing its traces). The plate
    # is a draft convenience here — content_hash stays NULL until commit.
    _resolve_series_plate!(db, sid, user_id, now_str)
    return sid
end

"""
    _update_view_for_series_recipe_updated!(db, entity_id, payload, event_id)

`series_recipe_updated` dispatcher (#166). Pure-replace: updates the recipe
scalars (`ordering_variable`, `order_rule`) on the `series` row, then
`DELETE`s every `series_samples` row for the series and re-`INSERT`s the full
payload snapshot. Never touches `content_hash` (recipe is excluded from the
hash — master plan §5.1), `state`, or `series_members`. `ordering_variable`
is a bare write (a PATCH is a full recipe replace, so an omitted field nulls
it); `order_rule` is `COALESCE`d because the column is `NOT NULL` and
`CHECK`-constrained.
"""
function _update_view_for_series_recipe_updated!(db, entity_id, payload, event_id)
    sid     = Int(entity_id)
    now_str = comparison_now_iso()
    ordering_var = haskey(payload, :ordering_variable) && payload.ordering_variable !== nothing ?
                   String(payload.ordering_variable) : nothing
    order_rule   = haskey(payload, :order_rule) && payload.order_rule !== nothing ?
                   String(payload.order_rule) : nothing

    # Bare UPDATE — assumes `series_created` has already folded this row.
    # Guaranteed by event-log ordering (a recipe edit cannot precede its
    # create); silently no-ops if folded standalone, which is acceptable
    # because `rebuild_views_from_log!` always replays the create first.
    DBInterface.execute(db,
        """UPDATE series
           SET ordering_variable = ?,
               order_rule = COALESCE(?, order_rule),
               updated_at = ?
           WHERE id = ?""",
        [ordering_var, order_rule, now_str, sid])

    DBInterface.execute(db, "DELETE FROM series_samples WHERE series_id = ?", [sid])
    samples = haskey(payload, :samples) ? payload.samples : []
    for s in samples
        _insert_series_sample!(db, sid, s)
    end
    return nothing
end

# Helper: INSERT one series_members row. Mirrors `_insert_comparison_member!`
# (the series_members and comparison_members column shapes are identical) —
# reuses `_member_field` / `_member_json`. Members in a series_plate_committed
# payload carry NO ids; the dispatcher mints fresh PKs.
function _insert_series_member!(db, series_id, m, user_id, now_str)
    snap = _member_json(m, :snapshot)
    snap === nothing &&
        error("series member missing required `snapshot` field")
    DBInterface.execute(db,
        """INSERT INTO series_members
             (series_id, exposure_id, display_order, band_height, y_offset,
              normalization, color_override, label_override,
              q_window_min, q_window_max, peak_display, snapshot,
              created_by, created_at)
           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)""",
        [Int(series_id),
         _member_field(m, :exposure_id),
         Int(_member_field(m, :display_order)),
         Float64(_member_field(m, :band_height; default=1.0)),
         Float64(_member_field(m, :y_offset; default=0.0)),
         String(_member_field(m, :normalization; default="none")),
         _member_field(m, :color_override),
         _member_field(m, :label_override),
         _member_field(m, :q_window_min),
         _member_field(m, :q_window_max),
         _member_json(m, :peak_display),
         snap,
         user_id, now_str])
    nothing
end

"""
    _resolve_series_plate!(db, series_id, user_id, now_str) -> Int

Resolve the **plate** (`series_members`) from the **recipe** (`series_samples`)
for a draft series. Pure-replace: `DELETE`s every existing `series_members` row,
then for each non-excluded recipe sample (in `position` order) resolves its
representative exposure and `INSERT`s one plate member with a freshly-computed
snapshot. Returns the number of members written.

Representative-exposure rule (the `_corpus_with_exposures` precedent,
`comparisons.jl`): highest-id `selected=1` exposure wins; else highest-id
exposure overall; else the sample resolves to NO exposure and is SKIPPED (a
recipe sample with no exposure has nothing renderable to plate). `display_order`
is sequential over the resolved members, so skipped rows leave no gap.

Used by the `series_created` dispatcher so a just-created draft lands with its
plate already resolved (the builder renders the waterfall immediately). Does
NOT touch `content_hash` or `state` — resolution is a draft convenience; the
plate is only frozen + hashed by `series_plate_committed` (the commit path).
"""
function _resolve_series_plate!(db, series_id, user_id, now_str)
    sid = Int(series_id)
    DBInterface.execute(db, "DELETE FROM series_members WHERE series_id = ?", [sid])

    recipe = Tables.rowtable(DBInterface.execute(db,
        """SELECT sample_id FROM series_samples
           WHERE series_id = ? AND excluded = 0
           ORDER BY position ASC, id ASC""", [sid]))

    display_order = 0
    for r in recipe
        sample_id = Int(r.sample_id)
        # Representative exposure: highest-id selected wins; else highest-id
        # overall; else nothing (skip — no renderable trace for this sample).
        exps = Tables.rowtable(DBInterface.execute(db,
            "SELECT id, selected FROM exposures WHERE sample_id = ? ORDER BY id ASC",
            [sample_id]))
        eid = nothing
        for e in Iterators.reverse(exps)
            if e.selected != 0
                eid = Int(e.id); break
            end
        end
        if eid === nothing && !isempty(exps)
            eid = Int(last(exps).id)
        end
        eid === nothing && continue

        # `_insert_series_member!` reads members via `getproperty` (the commit
        # path feeds it JSON3.Objects); round-trip the constructed dict through
        # JSON3 so property access resolves. Defaults (band_height, y_offset,
        # normalization, …) are supplied by `_member_field`'s `default=`.
        member = JSON3.read(JSON3.write(Dict{Symbol, Any}(
            :exposure_id   => eid,
            :display_order => display_order,
            :snapshot      => compute_member_snapshot(db, eid),
        )))
        _insert_series_member!(db, sid, member, user_id, now_str)
        display_order += 1
    end
    return display_order
end

"""
    _update_view_for_series_plate_committed!(db, entity_id, payload, event_id)

`series_plate_committed` dispatcher (#167). Pure-replace the plate:
`DELETE` every `series_members` row, then `INSERT` the full payload member
list (members carry no ids — `_insert_series_member!` mints them). Sets
`state='committed'` and computes `content_hash` from the plate via
`compute_series_content_hash` (which excludes the recipe — master plan §5.1).
Never touches `series_samples`.
"""
function _update_view_for_series_plate_committed!(db, entity_id, payload, event_id)
    sid     = Int(entity_id)
    user_id = user_id_for_event(db, event_id)
    now_str = comparison_now_iso()

    DBInterface.execute(db, "DELETE FROM series_members WHERE series_id = ?", [sid])
    members = haskey(payload, :members) ? payload.members : []
    for m in members
        _insert_series_member!(db, sid, m, user_id, now_str)
    end

    new_hash = compute_series_content_hash(db, sid)
    DBInterface.execute(db,
        "UPDATE series SET content_hash = ?, state = 'committed', updated_at = ? WHERE id = ?",
        [new_hash, now_str, sid])
    return nothing
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
function rebuild_views_from_log!(db::SQLite.DB, entity_id::Integer;
                                  entity_type::String = "exposure")
    rows = Tables.rowtable(DBInterface.execute(db,
        """SELECT id, action, entity_id, payload
           FROM user_actions
           WHERE entity_type = ? AND entity_id = ?
           ORDER BY id""", [entity_type, Int(entity_id)]))
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
