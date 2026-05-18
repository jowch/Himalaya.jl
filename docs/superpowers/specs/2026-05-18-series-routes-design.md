# Series REST routes + business logic — design (issue #165, I2.2)

**Status:** approved 2026-05-18
**Issue:** #165 — "Add the series REST routes and business logic (I2.2)"
**Depends on:** #164 (series schema, merged as I2.1 — `migrate_series!` in `db.jl`)
**Blocks:** #166–#168 (series event kinds), #170 (native round-trip test), #172, #173, #175

## 1. Purpose

Add the series REST API and business logic so the series data model has a full
route surface before any UI consumes it. Ships invisibly: empty `series` tables,
routes unconsumed by the frontend, zero user-facing change.

Two new files, both adapted from the Compare-page equivalents:

- `packages/HimalayaUI/src/series.jl` — business logic, adapted from `comparisons.jl`.
- `packages/HimalayaUI/src/routes_series.jl` — REST routes, mirrors `routes_comparisons.jl`.

Two deliberate departures from the Compare-page originals:

1. **The `is_author` editing gate is dropped** (architecture decision 3 — the
   recipe/plate split and the event log already do gating's job). It lives only
   at `routes_comparisons.jl:255,312`; it has no counterpart in the series routes.
2. **The `last_event_at` mixed-timestamp sort bug (#76) is fixed, not copied.**

**This route contract is the frozen Phase 2 deliverable.** Issue #175 re-audits
it; a post-merge drift finding spawns a fast-follow patch issue owned by this
issue's author.

## 2. Background — the two constraints that shape the work

### 2.1 Mutating routes are written in full but non-functional until #166–#168

The series mutating routes emit events via `apply_event!` with kinds
`series_created`, `series_recipe_updated`, `series_plate_committed`,
`series_deleted`, `series_pinned`, `series_unpinned`. The `update_view_for_event!`
dispatcher branches for those kinds do not land until #166–#168.

`update_view_for_event!` (`events.jl:305`) falls through to a **silent no-op**
for an unknown kind (`events.jl:426-427` — "default: no view update"). It does
**not** `error`. Consequently, a mutating series route merged before #166–#168:

- writes the durable `user_actions` log row (the event is recorded), and
- writes **no view rows** — the dispatcher no-ops the unknown kind.

The issue-map cards for #166–#168 list their files as `events.jl`,
`applyRemoteToCache.ts`, `mutatorRegistry.ts`, the new mutator files, and paired
contract tests — **not `routes_series.jl`**. Therefore the route file is written
once, here, in full, and **frozen**. The handlers carry real `apply_event!`
calls, real body validation, and the real `post_state` projection. "Stubbed"
describes their *behavior* (non-functional view-write until the dispatcher branch
exists), not their *code* (complete).

Two route groups are exceptions and work fully on merge of this issue:

- **Messages** (`GET`/`POST /api/series/{id}/messages`) — the `post_message`
  kind is already a known dispatcher branch (`events.jl:367`, a no-op), and the
  route writes the `series_messages` row directly, exactly as the comparison
  messages route does.
- **All `GET` routes and the pin-list read** — plain reads, no events.

### 2.2 The `last_event_at` sort bug (#76)

`comparisons_listing` (`comparisons.jl:669`), `comparisons_for_experiment`, and
`forks_of_comparison` all project:

```sql
COALESCE((SELECT MAX(ua.timestamp) FROM user_actions ua
          WHERE ua.entity_type = 'comparison' AND ua.entity_id = c.id),
         c.updated_at) AS last_event_at
...
ORDER BY last_event_at DESC, c.id DESC
```

`user_actions.timestamp` is space-separated (`2026-05-18 12:34:56`);
`c.updated_at` is `T`-separated ISO with a `Z` suffix
(`2026-05-18T12:34:56.000Z`, from `comparison_now_iso()`). Lexically `T` (0x54)
sorts above space (0x20), so a row whose `last_event_at` falls back to
`updated_at` always sorts above a row carrying a real `MAX(ua.timestamp)` —
regardless of actual recency. The `ORDER BY` is wrong.

The comparison code works around it on the client (`ComparisonSidebar` re-sorts);
`last_event_at` is documented there as a display hint, not a sort key.

**The fix in `series.jl`:** wrap the coalesced expression in SQLite's
`datetime()`, which parses both `YYYY-MM-DD HH:MM:SS` and
`YYYY-MM-DDTHH:MM:SS.SSSZ` into one canonical, lexically-sortable form:

```sql
datetime(COALESCE((SELECT MAX(ua.timestamp) FROM user_actions ua
                   WHERE ua.entity_type = 'series' AND ua.entity_id = s.id),
                  s.updated_at)) AS last_event_at
...
ORDER BY last_event_at DESC, s.id DESC
```

The projected `last_event_at` is the normalized value, so it is also a valid
client-side sort key — the bug is closed end-to-end, not worked around. **Both**
series listing queries — `series_listing` and `forks_of_series` — use this form.
(There is no per-experiment series listing; the folio is corpus-wide, §5.)

## 3. Files & wiring

- Create `packages/HimalayaUI/src/series.jl`.
- Create `packages/HimalayaUI/src/routes_series.jl`.
- `HimalayaUI.jl` — add `include("series.jl")` after `include("comparisons.jl")`
  (line 15) and `include("routes_series.jl")` after `include("routes_comparisons.jl")`
  (line 27).
- `server.jl` — add `register_series_routes!()` after `register_comparisons_routes!()`
  (line 129).

`series.jl` and `comparisons.jl` are both `include`d into the same `HimalayaUI`
module, so `series.jl` **reuses** the in-module generic helpers
(`canonical_json`, `user_id_for_event`, `comparison_now_iso`) directly rather
than duplicating them. Phase 3 (#175 / I3.6) relocates those helpers when it
deletes `comparisons.jl` — explicitly out of scope here. (`user_id_for_event`
and `canonical_json` must survive that deletion regardless: the frozen
`comparison_*` dispatcher branches in `events.jl` keep using them forever.)

## 4. `series.jl` — business logic

Adapted from `comparisons.jl`. New series-specific functions:

- **`series_listing(db)`** — corpus-wide listing for `GET /api/series`. Adapted
  from `comparisons_listing`. `entity_type='series'`; the `last_event_at` bug is
  fixed per §2.2. Reads from `series` / `series_members` (the plate, for
  `member_count` / `member_phases` / `has_stale_members`) and `users`. Includes
  zero-member and orphan-member series defensively, as the comparison original does.

- **`forks_of_series(db, id)`** — adapted from `forks_of_comparison`; rows where
  `forked_from_id = ?`. Same row shape as `series_listing`; same `datetime()` fix.

- **`fetch_series_with_plate(db, id)`** — adapted from
  `fetch_comparison_with_members`. Projects the `series` row + its `series_members`
  (the plate, ordered by `display_order`) + its `series_samples` (the recipe,
  ordered by `position`). Returns `nothing` when the series does not exist. This
  is the `post_state` projection that the `series_plate_committed` route emits.

- **`compute_series_content_hash(db, id)`** — adapted from `compute_content_hash`.
  **Plate only** — hashes `series_members`, never the recipe (`series_samples`),
  per master plan §5.1. `content_hash` is NULL while `state='draft'` and is
  computed on `series_plate_committed`.

- **`current_series_content_hash(db, id)`** — adapted from
  `current_content_hash`. Existence probe + stored hash; returns `nothing` for a
  missing series. Used by the mutating routes for the 404 existence check.

- Projection / normalization helpers, adapted from their comparison equivalents:
  - `_series_listing_rows(rows)` — the lightweight per-row listing shape.
  - `_series_member_payload(db, m)` — normalizes a plate-member entry from a
    request body (fills `snapshot` from `compute_member_snapshot` when omitted).
  - `_series_sample_payload(m)` — normalizes a recipe entry into the
    `series_samples` shape: `{sample_id, position, pinned, excluded}`.

**Not ported:** `is_author` — the editing gate is dropped (§1).

**Reused as-is from `comparisons.jl`** (same module, no duplication):
`canonical_json`, `user_id_for_event`, `comparison_now_iso`,
`compute_member_snapshot` (exposure-based, not comparison-specific).

## 5. `routes_series.jl` — route contract

`register_series_routes!()`, mirroring `register_comparisons_routes!()`.
**Corpus-wide** — there is no `/api/experiments/{eid}/series` route (the series
folio is corpus-wide; the per-experiment listing route has no series counterpart).

Every mutating route wraps in `with_idempotency(db, req)` and calls
`apply_event!(InTransaction(), …)` so the durable event row, any view-table
write, and the idempotency cache row commit or roll back atomically. Body-shape
validation happens **before** `with_idempotency` so a malformed payload returns
an uncached 400, not a cached 500 — the `routes_comparisons.jl` convention.

No route has an `is_author` / 403 gate. Existence (404) and optimistic-concurrency
(409) checks stay — they are not the author gate.

| Route | Method | Event kind | Functional in I2.2 |
|---|---|---|---|
| `/api/series` | `GET` | — | ✅ fully (`series_listing`) |
| `/api/series` | `POST` | `series_created` | degenerate |
| `/api/series/{id}` | `GET` | — | ✅ fully (`fetch_series_with_plate`) |
| `/api/series/{id}/forks` | `GET` | — | ✅ fully (`forks_of_series`) |
| `/api/series/{id}` | `PATCH` | `series_recipe_updated` | degenerate |
| `/api/series/{id}/commit` | `POST` | `series_plate_committed` | degenerate |
| `/api/series/{id}` | `DELETE` | `series_deleted` | degenerate |
| `/api/series/{id}/messages` | `GET` | — | ✅ fully |
| `/api/series/{id}/messages` | `POST` | `post_message` | ✅ fully |
| `/api/users/me/series-pins` | `GET` | — | ✅ fully |
| `/api/series/{id}/pin` | `POST` | `series_pinned` | degenerate |
| `/api/series/{id}/pin` | `DELETE` | `series_unpinned` | degenerate |

"degenerate" — the handler is written in full (real `apply_event!` call,
validation, `post_state`); it writes the `user_actions` log row but no view rows
until the dispatcher branch lands (§2.1).

### 5.1 `GET /api/series`
`series_listing(db)` → `200` JSON array. Empty array on empty tables.

### 5.2 `POST /api/series` — create a draft (`series_created`)
A series is born as a **draft** with a recipe and no plate — a new lifecycle
(`comparison_created` rejects empty members; `series_created` carries zero
`series_members`).

Body validation (before `with_idempotency`):
- `title` — required (mirrors `routes_comparisons.jl`); 400 if missing.
- `samples` — must be an array if present; **may be empty** (a draft can be
  scoped later). 400 if present and not an array.
- view-choice fields type-guarded via the existing `_view_fields_error` shape.

Inside `with_idempotency`:
- Mint the AUTOINCREMENT id with `INSERT INTO series DEFAULT VALUES`, capture
  `lastrowid`. (The placeholder row's `state` defaults to `'committed'` per the
  schema CHECK default; the `series_created` dispatcher in #166 upserts it to
  `'draft'`. Until that lands, a `POST /api/series` leaves a degenerate row at
  `state='committed'` — accepted, see §2.1.)
- Build the payload: `title`, `description`, recipe fields (`ordering_variable`,
  `order_rule`), the full `samples` snapshot (each entry normalized by
  `_series_sample_payload`), `forked_from_id`, `forked_at_hash`, view-choice
  fields. Zero members.
- `apply_event!(InTransaction(), …; kind="series_created", entity_type="series",
  entity_id=new_id, payload)`.
- Project once with `fetch_series_with_plate(db, new_id)` — this is both the SSE
  `post_state` envelope and the HTTP body (identical shapes, so an HTTP-wins /
  SSE-wins race converges). Enqueue the post-commit broadcast.
- `201` with the projection.

### 5.3 `GET /api/series/{id}`
`fetch_series_with_plate(db, id)` → `200`, or `404` when `nothing`.

### 5.4 `GET /api/series/{id}/forks`
`forks_of_series(db, id)` → `200` JSON array.

### 5.5 `PATCH /api/series/{id}` — recipe edit (`series_recipe_updated`)
Edits the recipe (`series_samples`). The payload carries the **full** `samples`
snapshot — never a delta (the dispatcher pure-replaces). May also carry the
recipe fields `ordering_variable` / `order_rule`.

- 404 if `current_series_content_hash(db, id) === nothing`.
- No author gate.
- `apply_event!(…; kind="series_recipe_updated", entity_type="series",
  entity_id=id, payload)`; the payload `samples` entries normalized by
  `_series_sample_payload`.
- Project with `fetch_series_with_plate`; broadcast; `200`.

### 5.6 `POST /api/series/{id}/commit` — commit the plate (`series_plate_committed`)
The series equivalent of the old comparison `submit`. The payload carries the
**full member list** (the plate); members carry **no ids** (the dispatcher mints
them). Sets `state='committed'` and computes `content_hash` from the plate.

- 404 if `current_series_content_hash(db, id) === nothing` (existence before any
  other check — HTTP requires 404 before 409).
- **Optimistic-concurrency check (kept):** if the body carries
  `expected_content_hash` and it does not match
  `current_series_content_hash(db, id)`, return `409` with
  `{error:"conflict", current_hash, current_state}`. This is *not* the author
  gate — it is the same optimistic-concurrency guard the comparison `submit`
  carries, and the Phase 3 series builder reuses `ConflictModal` against it.
- `apply_event!(…; kind="series_plate_committed", entity_type="series",
  entity_id=id, payload)`; members normalized by `_series_member_payload`.
- Project with `fetch_series_with_plate` — the `post_state` envelope and HTTP
  body. Broadcast; `200`.

### 5.7 `DELETE /api/series/{id}` (`series_deleted`)
- 404 if `current_series_content_hash(db, id) === nothing`.
- No author gate.
- `apply_event!(…; kind="series_deleted", entity_type="series", entity_id=id,
  payload=Dict(:id => id))`. The dispatcher (#166) is a one-line
  `DELETE FROM series` — the four child tables cascade.
- Broadcast; `200` with `{id, deleted:true, event_id}`.

### 5.8 Messages — `GET` / `POST /api/series/{id}/messages`
Mirrors `routes_comparisons.jl` exactly; **fully functional in I2.2**.

- `GET` — `SELECT … FROM series_messages m LEFT JOIN users u … WHERE
  m.series_id = ? ORDER BY m.created_at ASC, m.id ASC`. (Sorting messages among
  themselves — all one timestamp format; safe. Per the `migrate_series!`
  docstring caveat, `series_messages.created_at` must never be string-sorted
  against `series` / `series_members` timestamps — this query does not.)
- `POST` — requires the `X-Username` header (401 if absent); non-empty body
  (400 if empty). Inside `with_idempotency`: `get_or_create_user!`, `INSERT INTO
  series_messages` directly, then `apply_event!(…; kind="post_message",
  entity_type="series_message", entity_id=msg_id, payload=msg_json)`. The
  dispatcher is a no-op for `post_message`. Broadcast; `201` with the message JSON.

### 5.9 Pins — `GET /api/users/me/series-pins`, `POST`/`DELETE /api/series/{id}/pin`
Mirrors the comparison pin routes. Pin events are stored with
`entity_type='user'`, `entity_id=user_id` (the `comparison_pinned` precedent) —
the affected series rides in the payload as `series_id`.

- `GET /api/users/me/series-pins` — plain read; 401 without a user;
  `SELECT series_id FROM series_pins WHERE user_id = ? ORDER BY pinned_at DESC,
  series_id DESC` → array of ids. **Fully functional in I2.2.**
- `POST /api/series/{id}/pin` — 401 without a user; inside `with_idempotency`,
  404 if the series does not exist; `apply_event!(…; kind="series_pinned",
  entity_type="user", entity_id=user_id, payload=Dict(:series_id => id))`;
  broadcast; `200` with `{series_id, pinned:true}`. Degenerate until #168.
- `DELETE /api/series/{id}/pin` — symmetric; `kind="series_unpinned"`; idempotent
  at the SQL layer; `200` with `{series_id, pinned:false}`. Degenerate until #168.

## 6. Tests

`packages/HimalayaUI/test/` — a `routes_series` test file, registered in the
suite the way the existing route tests are. Scope follows the issue-map I2.2 card
("GET routes round-trip on empty `series` tables; `last_event_at` sort test;
Julia route tests"); the six-layer paired contract tests belong to #166–#168 and
the `rebuild_views_from_log!` round-trip belongs to #170 — they cannot be written
here because the mutating routes produce no view rows until the dispatcher
branches land.

I2.2 tests:

1. **GET round-trips on empty tables** — `GET /api/series` → `[]`;
   `GET /api/series/{id}` for a missing id → `404`; `GET /api/series/{id}/forks`
   → `[]`.
2. **`last_event_at` sort test** (headline) — seed `series` rows and
   `user_actions` rows (`entity_type='series'`) with **mixed** timestamp formats
   (space-separated and `T`-separated-with-`Z`); assert `series_listing` orders
   them by true recency, i.e. the #76 bug does not recur. Cover the `updated_at`
   fallback path (a series with no events).
3. **Messages full round-trip** — `POST /api/series/{id}/messages` then `GET`
   returns the message; 401 without `X-Username`; 400 on empty body.
4. **`GET /api/users/me/series-pins`** → `[]` for a fresh user; 401 without a user.
5. **Mutating-route smoke tests** — each of `POST /api/series`, `PATCH`,
   `POST .../commit`, `DELETE`, `POST`/`DELETE .../pin`: route is registered and
   reachable; returns 4xx on a malformed body (e.g. `POST` without `title`);
   returns `404` on a missing id where applicable; a `user_actions` row is
   written for a well-formed request.
6. **No `is_author` gate** — a user who is not the series author can reach
   `PATCH` / `commit` / `DELETE` without a `403` (the request returns its normal
   2xx/404, never 403). This proves the gate is gone — under the comparison
   routes the equivalent request returns 403.

## 7. Acceptance criteria (from issue #165)

- [ ] All series routes round-trip on empty `series` tables (GET routes and
      messages fully; mutating routes at the HTTP level — §2.1).
- [ ] Series list ordering is correct — the `last_event_at` bug does not recur —
      covered by the §6.2 test.
- [ ] There is no `is_author` gate (§6.6).
- [ ] Julia route tests pass.
- [ ] The running frontend is unaffected (routes unconsumed by the UI).
- [ ] The route contract is the frozen Phase 2 deliverable #175 re-audits.

## 8. Out of scope

- Series event kinds — the dispatcher branches, SSE handlers, mutators, and
  six-layer contract tests (#166–#168).
- The native `rebuild_views_from_log!` round-trip test (#170).
- The series UI (Phase 3).
- The comparison→series data migration (#172 / master plan §6.1).
- Relocating the shared helpers out of `comparisons.jl` (Phase 3 / I3.6).

## 9. References

- Master plan §5.2 (routes, event kinds, gate):
  `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign.md`
- Issue map I2.2: `docs/superpowers/plans/2026-05-17-himalaya-ui-redesign-issue-map.md`
- Sort bug origin: `packages/HimalayaUI/src/comparisons.jl:669` (issue #76)
- Dropped gate: `packages/HimalayaUI/src/routes_comparisons.jl:255,312`
- Series schema: `migrate_series!` in `packages/HimalayaUI/src/db.jl:723`
- Compare-page originals: `comparisons.jl`, `routes_comparisons.jl`
