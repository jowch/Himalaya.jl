---
name: new-event-kind
description: Scaffold a new event kind in HimalayaUI's apply_event! pipeline — dispatcher branch (if view-producing), route call site, rebuild_views_from_log! round-trip test (if view-producing), and frontend sseSubscriber wiring (if view-producing). Reads docs/event-log.md and existing call sites as live reference. Usage: /new-event-kind <kind> [--no-view]
---

# new-event-kind

Scaffolds the wiring for a new event kind that flows through `apply_event!`. The dispatcher contract — *only `update_view_for_event!` writes to view tables* — is easy to violate silently, so this skill walks every step instead of assuming.

## Args

- `<kind>` — the event kind in snake_case (e.g. `peak_pinned`, `index_starred`). Use past-tense verbs to match existing kinds (`peak_added`, `index_confirmed`).
- `--no-view` — pure log event, no materialized view update. Skips dispatcher and rebuild test steps. Use for instrumentation events like `analyze_run` or for state changes that already have their own table (e.g. `set_status`).
- `--undoes <prior-kind>` — this is an undo of a prior event kind. Wires `undoes_event_id` resolution.

## Procedure

### 1. Read the current reference implementations

Before writing any code, read these files to understand the *current* shape — the dispatcher contract evolves and stale templates are worse than no template:

```
docs/event-log.md                                  ← THE CONTRACT — read §1, §2, §5 in full
packages/HimalayaUI/src/events.jl                  ← apply_event! + dispatcher + broadcast
packages/HimalayaUI/src/routes_peaks.jl            ← canonical apply_event! call sites:
                                                   #   peak_added (INSERT branch, uses view_row_id)
                                                   #   peak_excluded (INSERT branch, dedup)
                                                   #   peak_unexcluded (DELETE branch, undoes_event_id)
                                                   #   peak_removed (--no-view example)
packages/HimalayaUI/test/test_events.jl            ← rebuild_views_from_log! property test
packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts  ← handleCurationEvent + invalidation
packages/HimalayaUI/frontend/test/sse.test.tsx     ← subscriber test pattern
CLAUDE.md                                          ← Plan 7 / event-log gotchas
```

If the new kind is view-producing, `peak_added` is the closest template (INSERT into a view table, capture `view_row_id` for the route's response). If it's a DELETE, `peak_unexcluded` is the template.

### 2. Decide the entity

Every `apply_event!` call needs `entity_type` and `entity_id`. Use `"exposure"` whenever the event affects exposure-scoped state — this is what `idx_events_by_exposure` indexes, and it's what `rebuild_views_from_log!` folds over. Other entity types are valid (e.g. `"sample"`, `"experiment"`) but break the per-exposure rebuild contract; only use them if the event truly isn't exposure-scoped.

The fixup commit in PR #29 explicitly migrated `peak_removed` / `speculative_created` / `speculative_deleted` from `entity_type='peak'` to `entity_type='exposure'` for exactly this reason. Stay on the rails.

### 3. Design the payload

The payload becomes a `JSON3.Object` inside the dispatcher (`apply_event!` round-trips through `JSON3.write` → `JSON3.read` so live and replay see the same shape). Branch code accesses fields uniformly via `payload.q`-style.

Rules:
- **Use symbol keys in the route**: `Dict(:q => q, :group_id => gid)` — JSON3 emits string keys; the dispatcher reads back as `JSON3.Object` which accepts both.
- **Include enough context for replay.** `rebuild_views_from_log!` re-folds payloads in id-order — the dispatcher must be able to reproduce the view write from payload alone. If the dispatcher needs to query `auto_peaks` to resolve a reference, fine; just don't depend on the *current* state of a view it's about to write to.
- **Never include actor/user data in payload.** That's already on the `user_actions` row (`user_id`, `created_at`). Payload is what the dispatcher needs.

### 4. Add the dispatcher branch (skip if `--no-view`)

In `packages/HimalayaUI/src/events.jl`, add a branch to `update_view_for_event!`:

```julia
if kind == "<your-kind>"
    res = DBInterface.execute(db,
        """INSERT INTO <view_table> (...)
           VALUES (...)""",
        [...])
    return Int(DBInterface.lastrowid(res))
end
```

For DELETE branches return `nothing`. For INSERT branches return `Int(DBInterface.lastrowid(res))` — the route uses this via `view_row_id` to avoid a re-query race (see `docs/event-log.md` §1 / Issue 2 in the PR #29 review followup).

### 5. Add the route call site

In whichever `routes_*.jl` file owns the entity, replace direct view writes + `log_action!` with a single `apply_event!`:

```julia
result = apply_event!(db, req;
    kind        = "<your-kind>",
    entity_type = "exposure",
    entity_id   = exposure_id,
    payload     = Dict(:field => value))

# If view-producing INSERT and you need the row id:
new_id = result.view_row_id  # never re-query by content
```

If the route does multiple writes that must be atomic (e.g. speculative create + delete), wrap the whole block in `SQLite.transaction(db) do ... end` — `apply_event!`'s inner transaction nests safely under SQLite.jl's savepoint semantics.

**`--undoes` wiring.** If this kind undoes a prior event (e.g. `peak_unexcluded` undoes `peak_excluded`), resolve the prior event's id at the route layer and pass it as `undoes_event_id`. Pattern from `routes_peaks.jl:158`:

```julia
# Find the matching prior event (e.g. by exposure + payload q within tolerance)
undoes = find_prior_event_id(db, exposure_id, prior_kind, q)

apply_event!(db, req;
    kind            = "<your-kind>",
    entity_type     = "exposure",
    entity_id       = exposure_id,
    payload         = Dict(:q => q),
    undoes_event_id = undoes)  # may be nothing if not resolvable
```

The dispatcher branch for an undo-kind is typically a DELETE — see `peak_unexcluded` in `events.jl` for the canonical shape.

### 6. Add a `rebuild_views_from_log!` round-trip test (skip if `--no-view`)

In `packages/HimalayaUI/test/test_events.jl`, mirror the existing `peak_added` + `peak_excluded` test pattern. The property: starting from empty views, applying every event in order produces the same state as live `apply_event!` did.

```julia
@testset "<your-kind> round-trips through rebuild_views_from_log!" begin
    db = open_test_db()
    # ... set up exposure, fire apply_event! a few times with this kind ...
    live_state = Tables.rowtable(DBInterface.execute(db, "SELECT ... FROM <view_table> WHERE ..."))

    # Wipe the view, rebuild from the log
    DBInterface.execute(db, "DELETE FROM <view_table> WHERE <scope>")
    HimalayaUI.rebuild_views_from_log!(db, exposure_id)
    rebuilt_state = Tables.rowtable(DBInterface.execute(db, "SELECT ... FROM <view_table> WHERE ..."))

    @test live_state == rebuilt_state
end
```

If the test fails, the dispatcher is reading from state instead of payload — fix the dispatcher, not the test.

### 7. Add the route test

In `packages/HimalayaUI/test/test_routes_<resource>.jl`:
- Verify the route returns the expected status + body
- Verify the view row was created with the right fields
- Verify a `user_actions` row was created with `kind`, `entity_type='exposure'`, `entity_id=...`, `payload` JSON-decodable to the expected shape
- For INSERT branches: verify the response id matches the actual `view_row_id` (not re-queryable by content)

### 8. Frontend: extend `handleCurationEvent` (skip if `--no-view`)

In `packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts`, the existing handler invalidates `peaks(id)`, `indices(id)`, `groups(id)`, `exposure(id)` on every `entity_type === "exposure"` event. If the new kind affects a query key not in this list (e.g. a brand-new `tags(id)` resource), add it to the invalidation set.

If the new kind is on a different entity type (`"sample"`, `"experiment"`), `handleCurationEvent` currently early-returns — extend it with a parallel branch and add to `test/sse.test.tsx`.

### 9. Verify

```bash
# Backend
julia --project=packages/HimalayaUI -e 'using Pkg; Pkg.test("HimalayaUI")'
# Or just the affected files (with sysimage):
julia --sysimage build/himalaya.so --project=packages/HimalayaUI -e 'using Test, HimalayaUI; include("packages/HimalayaUI/test/test_events.jl"); include("packages/HimalayaUI/test/test_routes_<resource>.jl")'

# Frontend
(cd packages/HimalayaUI/frontend && npm test -- sse.test.tsx)
```

The PostToolUse hooks in `.claude/settings.json` will auto-fire matching tests on each edit if a sysimage is present.

## Anti-patterns to avoid

- **Don't write directly to `peak_curations` or `index_group_members`** anywhere outside `update_view_for_event!`. The whole point of the contract is that it's a single sole-writer; `rebuild_views_from_log!` will silently produce wrong state otherwise.
- **Don't re-query by content after an INSERT** to find the row id. Use `result.view_row_id`. Concurrent writers with the same content will both see the larger id (the bug PR #29's fix-3 commit explicitly fixed).
- **Don't broadcast manually.** `apply_event!` calls `broadcast_event!` for you, after the transaction commits. Calling it twice will double-fire.
- **Don't put dispatcher writes outside the `apply_event!` transaction.** The dispatcher runs inside the SQLite transaction opened by `apply_event!` so the log row and view row commit together. Writes outside this transaction can leave the log and view inconsistent on failure.
- **Don't use `log_action!` for view-producing events.** It bypasses the dispatcher and skips the broadcast. Use `apply_event!`. Reserve `log_action!` for the legacy non-view kinds explicitly listed in `docs/event-log.md` §1.

## When NOT to use this skill

- Adding a new column to `user_actions` itself — that's a schema migration, use the SQLite migration pattern in `db.jl`.
- Adding a new view table from scratch — design the schema and dispatcher contract together; this skill assumes the view table exists.
- Pure read endpoints — they don't need events. Use `/new-route` instead.
