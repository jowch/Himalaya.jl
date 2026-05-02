# Event log, hash memoization, and SSE multiplayer

This is the load-bearing reference for the Plan 7 architecture in
`packages/HimalayaUI`. Three coupled concerns:

1. **Event log + dispatcher** — `user_actions` is the source of truth for
   curation; view tables (`peak_curations`, `index_group_members`) are
   materialized exclusively by one dispatcher function.
2. **Hash memoization** — `findpeaks` and `indexpeaks` skip when their
   inputs are bit-identical to the last successful run.
3. **SSE multiplayer** — every successful event broadcasts to every
   connected client; clients invalidate scoped TanStack Query keys for
   the affected exposure.

Read this before touching `events.jl`, `hash.jl`, the `apply_event!` call
sites in `routes_*.jl`, the SSE handler in `server.jl`, or the
`StaleIndicesBanner` gating logic in the frontend.

---

## 1. The dispatcher contract

> **`update_view_for_event!` is the sole writer to view tables.**

Concretely: no route, migration, or pipeline function may `INSERT` /
`DELETE` into `peak_curations` or `index_group_members` except by going
through `apply_event!`. The dispatcher is the only function with
SQL-against-view-tables in `events.jl`, and the rest of the codebase is
expected to honour that.

Why it matters: `rebuild_views_from_log!` re-folds every event for an
exposure from `user_actions` and asserts the resulting view state matches
what live writes produced. If a route writes directly to
`peak_curations`, the rebuild reproduces the wrong state and the property
test catches it — but only for the event kinds the test exercises. New
view-writing event kinds need a corresponding dispatcher branch *and* a
round-trip test in `test_events.jl`.

### Event kinds today

View-producing (route through `apply_event!`):

| Kind | Dispatcher action |
|---|---|
| `peak_added` | `INSERT INTO peak_curations(kind='add', q=…)` |
| `peak_excluded` | `INSERT INTO peak_curations(kind='exclude', q=…)` |
| `peak_unexcluded` | `DELETE FROM peak_curations WHERE kind='exclude' AND q≈payload.q` (with `undoes_event_id` set when resolvable) |
| `index_confirmed` | `INSERT OR IGNORE INTO index_group_members(group_id, index_id)` |
| `index_unconfirmed` | `DELETE FROM index_group_members WHERE group_id=… AND index_id=…` |

Logged but no view update (use `log_action!`):

`set_status`, `add_tag`, `remove_tag`, `add_message`, `update_sample`,
`update_experiment`, `analyze`, `reingest`, `analyze_run`,
`peak_removed`, `speculative_created`, `speculative_deleted`.

`peak_removed` / `speculative_*` go through `apply_event!` for entity-type
discipline (so `idx_events_by_exposure` works) but produce no view write
— the dispatcher returns `nothing` for these.

### Payload contract

The dispatcher sees `payload` as a `JSON3.Object` (or `nothing`),
re-parsed from `payload_json` *inside* `apply_event!` so the live path
and `rebuild_views_from_log!` see exactly the same shape. Branch code
accesses fields uniformly via `payload.q` style — `JSON3.Object` supports
both symbol and string keys so the symbol-vs-string footgun is
neutralised at the boundary.

### Atomicity

`apply_event!` opens a single `SQLite.transaction` covering both the
`user_actions` insert and the dispatcher's view writes. Either both
commit or neither does. Routes that wrap multiple `apply_event!` calls
(e.g. speculative create + delete) must open an outer
`SQLite.transaction` themselves so the entire route is atomic; see
`routes_peaks.jl` and `routes_analysis.jl` for the pattern.

### Return shape

```julia
(event_id, view_row_id) = apply_event!(db, req; kind, entity_type, …)
```

`view_row_id` is the lastrowid of any INSERT the dispatcher performed,
or `nothing` for DELETE / no-op branches. Routes that need the inserted
row id (e.g. `POST /api/peaks` returning the new curation id) use this
field directly instead of re-querying by `q` — the re-query path was
prone to a concurrent-add collision where two same-q POSTs both saw the
larger id.

---

## 2. Hash memoization

Three SHA-256 hashes drive skip-when-unchanged:

| Column | Source | Skips when matched |
|---|---|---|
| `exposures.trace_hash` | `hash_trace_file(.dat path)` — bytes of the integration file | `findpeaks` |
| `exposures.analysis_inputs_hash` | `hash_peak_set(eff)` — sorted `(q, sharpness)` Float64 tuples | `indexpeaks` |
| `indices.inputs_hash` | snapshot of the `analysis_inputs_hash` that produced this index | drives `StaleIndicesBanner` |

**Skip predicate (both stages):**

```
skip = stored_hash == new_hash AND row_count > 0
```

The `row_count > 0` clause is load-bearing. A bare `stored_hash == new_hash`
check would short-circuit on a fresh DB whose hashes happened to match
defaults (or on a DB where the rows had been deleted out-of-band) without
actually having anything to read back. The `analyze_run` event payload
records the full skip predicate's outcome — instrumentation must reflect
the real branch taken, not the hash comparison alone.

**Staleness banner:** `StaleIndicesBanner` renders when *any* index's
`inputs_hash` differs from its exposure's current `analysis_inputs_hash`.
This replaces the retired `status='stale'` enum — there is no longer any
write-side "mark all indices stale" call. Staleness is purely derived.

**`hash_peak_set` is order-independent.** Inputs are sorted by `q` before
hashing, so an exclude-then-add or add-then-exclude that lands on the
same final peak set produces the same hash. This is the contract that
makes "no-op rerun → skip indexpeaks" hold across curation orderings.
There is a regression test for this in `test_pipeline.jl`.

**Hash migration.** Pre-Plan-7 DBs have NULL hashes; the first
`analyze_exposure!` on each exposure populates them. No backfill is
needed.

---

## 3. SSE multiplayer

### Server side (`server.jl`)

```
GET /api/events
  Content-Type: text/event-stream
  Cache-Control: no-cache
  X-Accel-Buffering: no       ← required for nginx
```

Architecture: Oxygen.jl 1.10.x `@stream` macro receives an `HTTP.Stream`
directly. Each subscriber gets a per-task `Channel{String}` for pending
frames plus a per-subscriber `Timer` that pushes `:heartbeat\n\n` every
15 s. The handler loop is just `for frame in pending` — race-free, no
busy-poll.

Both headers matter for deployment:

- `Cache-Control: no-cache` keeps proxies from caching the stream.
- `X-Accel-Buffering: no` defeats nginx's response buffering. Without
  this, nginx batches frames until the buffer fills and the multiplayer
  experience feels broken even though the backend is fine.

### Broadcast semantics

`broadcast_event!` fires from `apply_event!` *after* the transaction
commits. Implications:

- Subscribers never see an event that was rolled back.
- Process death between commit and broadcast loses the frame, but the
  event is durable in `user_actions`. Clients reconcile on EventSource
  reconnect via TanStack Query refetch — `App.tsx` invalidates the
  exposure's query keys on the first reconnect.
- Slow subscribers (channel full) and disconnected subscribers (channel
  closed) are pruned via `_try_put!` rather than blocking the broadcast
  loop. A pruned subscriber reconnects on the EventSource side and gets
  a fresh subscription.

### Client side (`src/lib/sseSubscriber.ts`)

`handleCurationEvent(data, { username, qc })`:

1. Parse the JSON frame; ignore on parse error or missing `entity_id`.
2. **Self-echo filter.** Skip if `event.actor === username`. This means
   two browser tabs sharing the same `X-Username` (same lab user) both
   filter their own echoes — by design; lab users with two tabs are
   collaborating with themselves, not racing.
3. Skip if `entity_type !== "exposure"` (defensive — only exposure
   events update view caches today).
4. Invalidate `peaks(id)`, `indices(id)`, `groups(id)`, `exposure(id)`
   for the affected exposure id. TanStack Query refetches what's
   currently mounted; nothing happens for queries the user can't see.

The EventSource connection is bound to `App.tsx`'s mount/unmount only;
`username` is read through a `useRef` inside the listener so onboarding
(`undefined → "alice"`) doesn't recycle the connection.

### Conflict resolution (deferred)

Optimistic concurrency via `If-Match` headers + 409-retry on the
frontend (R5b in Plan 7) is **deferred**. The gate is R4 instrumentation
showing real contention: ≥2% delta-event collision rate over ≥4 weeks /
≥500 events. Until then, multiplayer is last-write-wins, which the
event log makes auditable after the fact.

`exposures.selected` is intentionally LWW even after R5b ships — see the
note in `CLAUDE.md`.

---

## 4. Operational notes

- **Reverse proxy.** Set `proxy_buffering off` on the `/api/events`
  location in nginx, or rely on the `X-Accel-Buffering: no` header. SSE
  needs flushed-per-frame delivery.
- **Long-lived connections.** Each subscriber holds an open HTTP
  connection. Plan for connection limits (default Oxygen / HTTP.jl
  config is generous; nothing in HimalayaUI imposes a tighter cap).
- **Event log growth.** `user_actions` grows monotonically. No retention
  policy today. At lab scale (single-digit events per exposure, ~hundreds
  of exposures per experiment) this is a non-issue for years; revisit
  when an experiment crosses ~1M events.
- **Disaster recovery.** `rebuild_views_from_log!(db, exposure_id)` will
  rebuild `peak_curations` and `index_group_members` for one exposure
  from the event log. Useful if a view table is ever corrupted.

---

## 5. Adding a new view-producing event kind

1. Add a dispatcher branch in `update_view_for_event!`. Use
   `JSON3.Object` payload access (`payload.field`).
2. Make the originating route call `apply_event!` (not `log_action!`)
   with `entity_type='exposure'` so `idx_events_by_exposure` indexes it.
3. Add the kind to `rebuild_views_from_log!`'s round-trip test in
   `test_events.jl`. The property is: live `apply_event!` writes vs
   "wipe view table → fold every event" produce the same state.
4. If the dispatcher's INSERT exposes a row id the route needs, capture
   it via the `view_row_id` field on `apply_event!`'s return tuple
   instead of re-querying.
5. Frontend: if the new kind affects a query key not already in
   `handleCurationEvent`, add it there.

---

## See also

- [`docs/himalayaui-design.md`](himalayaui-design.md) §2.6, §2.9, §2.10 —
  design principles and the active-set preservation story.
- [`docs/superpowers/specs/2026-05-01-multiplayer-instrumentation-design.md`](superpowers/specs/2026-05-01-multiplayer-instrumentation-design.md) — original design spec.
- [`docs/superpowers/plans/2026-05-01-multiplayer-instrumentation.md`](superpowers/plans/2026-05-01-multiplayer-instrumentation.md) — implementation plan with R0–R5a phases.
- `packages/HimalayaUI/src/{events,hash,server,pipeline}.jl` and
  `packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts` — the code.
