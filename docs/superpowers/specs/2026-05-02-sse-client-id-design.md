# Per-Tab SSE Client Identity — Design Spec

**Date:** 2026-05-02
**Status:** Draft
**Related:** [docs/event-log.md](../../event-log.md) §"Client side", [docs/superpowers/specs/2026-05-01-multiplayer-instrumentation-design.md](2026-05-01-multiplayer-instrumentation-design.md) (Plan 7 R5)

## Context

The current SSE self-echo filter ([sseSubscriber.ts:32](../../../packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts:32)) compares `event.actor` to the local `username`. Two browser tabs signed in as the same `X-Username` therefore filter each other's events — by design, on the assumption that "two tabs of the same user = one human collaborating with themselves." Documented in [docs/event-log.md:179-188](../../event-log.md).

The Compare page (in design) invalidates that assumption. Compare's central workflow is multi-pane viewing of several exposures at once, naturally combined with a second tab on Inspect or Index for drill-down on any one of them. With the current filter:

- Tab A (Compare) shows exposures X, Y, Z.
- Tab B (Inspect) — same user — edits a peak on exposure Y.
- Server broadcasts with `actor: "alice"`.
- Tab A's filter drops the frame as a self-echo.
- Compare's panel for Y silently shows stale peaks/indices while the user is actively comparing.

The documented escape hatch (`alice-laptop` vs `alice-desktop`) doesn't fit this workflow — it's not two devices, it's two tabs on one screen, intentionally.

The fix is to **separate audit identity from routing identity**: keep `actor` as the username for the event log and chat mentions, add a per-tab `client_id` for SSE filter purposes only.

## Goals

- Multi-tab same-user workflows (Compare + Inspect, multi-monitor split) see live updates across tabs without sacrificing self-echo suppression within a tab.
- Audit identity (`user_actions.actor`, chat `@mentions`, `index_groups.created_by`) stays human-grained — one row in `users`, one identity for `@alice`.
- Change is single-PR, additive at the schema level, ships independently of any other Plan 7 work.

## Non-goals

- Authentication or identity hardening. `client_id` is a routing tag, not a credential — it's no more spoofable than `X-Username` is today, and the trust model is unchanged.
- Cross-tab state sync beyond SSE-driven cache invalidation. Tabs still maintain independent TanStack Query caches; the change just makes them refetch when relevant events arrive.
- Persistent client identity across browser sessions. `client_id` lives in `sessionStorage` — it survives reload but not closing the tab. This is sufficient for SSE filtering and simpler than the long-lived alternative.
- Conflict resolution (`If-Match`/409 retry — Plan 7 R5b). Still gated on R4 instrumentation; this spec doesn't change that calculus.
- Reconciling the existing `actor`-based filter as a fallback. Once `client_id` ships, it is the sole filter mechanism. The two-mode complexity isn't worth the marginal robustness.

## Pain points in the current design

1. **`actor` is overloaded.** It's the audit identity (queryable in `user_actions`, displayed in chat, FK-able from `index_groups.created_by`) *and* the SSE filter discriminator. Coupling audit to routing means the filter granularity is fixed at the user level, where it should be at the connection level.
2. **`sessionStorage` exists for this exact use case** but isn't yet leveraged. Per-tab durability — survives reload, scoped to one tab — is the natural fit for SSE routing identity.
3. **Compare's UX assumes freshness across panes.** With Compare loading 2-N exposures simultaneously, the cost of a missed broadcast scales with the number of compared items. Stale-on-focus refetch (the current fallback) doesn't help here because the user isn't switching focus — they're staring at a stale Compare grid.

## Design

### Identity model

Two orthogonal identities per request:

| Field | Source | Lifetime | Purpose |
|-------|--------|----------|---------|
| `actor` (`username`) | `X-Username` header → `users.username` | Account-scoped | Audit log, chat mentions, `created_by` FKs |
| `client_id` | `X-Client-Id` header → ephemeral | Tab session | SSE self-echo filter |

`actor` continues to do what it does today; nothing about audit, chat, or FKs changes. `client_id` is **broadcast on the SSE frame** and **persisted on `user_actions`** (nullable, additive column) so the event log records "which tab originated this event" — useful for telemetry and for replay correctness if a future feature broadcasts from `rebuild_views_from_log!`.

### Schema

Single additive change:

```sql
ALTER TABLE user_actions ADD COLUMN client_id TEXT;  -- nullable, no backfill
```

Append to the `stmts` array in `migrate_schema!` ([db.jl:199-215](../../../packages/HimalayaUI/src/db.jl)) — matches the established pattern for the other ~14 ALTER statements there, all wrapped in a single try/catch loop. Also add `client_id TEXT` to the `user_actions` `CREATE TABLE` DDL in `create_schema!` so fresh DBs get the column without relying on the migration. (This dual-add is the project's standard pattern — `payload`, `trace_hash`, `analysis_inputs_hash` all live in both places.) No separate `migrate_rN_*` helper; that ceremony is reserved for table-rebuild migrations.

NULL `client_id` is meaningful and persistent: it represents events with no originating tab — system-emitted events from `_system_request()` (pipeline runs), CLI-emitted events, or pre-migration rows. The SSE filter treats NULL `client_id` as "broadcast to all" (no match against any tab's id), which is the correct behavior for system events.

### Backend

**`actions.jl`** — new helper, mirroring `get_username`:

```julia
function get_client_id(req::HTTP.Request)
    v = HTTP.header(req, "X-Client-Id", "")
    isempty(v) ? nothing : String(v)
end
```

**`events.jl::apply_event!`** ([events.jl:19-67](../../../packages/HimalayaUI/src/events.jl:19)):

- Extract `client_id = get_client_id(req)` alongside `username`.
- Add `client_id` to the `user_actions` INSERT column list and parameter binding.
- Pass `client_id` to `broadcast_event!` as a new positional argument.

**`events.jl::broadcast_event!`** ([events.jl:199-222](../../../packages/HimalayaUI/src/events.jl:199)):

- New `client_id::Union{String, Nothing}` parameter.
- Include `:client_id => client_id` in the JSON Dict written to the SSE frame.

**`events.jl::rebuild_views_from_log!`** ([events.jl:239-254](../../../packages/HimalayaUI/src/events.jl:239)) — unchanged. `client_id` doesn't affect view dispatch; it's pure routing metadata. SELECT list already projects only `id, action, entity_id, payload`.

**`_system_request()`** ([events.jl:76](../../../packages/HimalayaUI/src/events.jl:76)) — unchanged. Synthetic system requests have no `X-Client-Id` header, so `get_client_id` returns nothing, and pipeline-emitted events broadcast to all tabs (correct behavior).

### Frontend

**Mint + persist** — new `clientId.ts` module (or addition to [main.tsx](../../../packages/HimalayaUI/frontend/src/main.tsx)):

```ts
export function getClientId(): string {
  let id = sessionStorage.getItem("himalaya.client_id");
  if (!id) {
    id = crypto.randomUUID();
    sessionStorage.setItem("himalaya.client_id", id);
  }
  return id;
}
```

Called once at module load; the resulting id lives for the tab's session.

**`api.ts`** ([api.ts:42-52](../../../packages/HimalayaUI/frontend/src/api.ts:42)):

- Extend `AuthOpts` with `clientId?: string`.
- Send `X-Client-Id` header on mutating calls (same condition as `X-Username`: `method !== "GET"`).

**`queries.ts`** ([queries.ts:96-98](../../../packages/HimalayaUI/frontend/src/queries.ts:96)):

- Expand `authOpts(username)` to `authOpts(username, clientId)`. All call sites already go through this helper, so the wiring is centralized.
- A module-level `const CLIENT_ID = getClientId()` in `queries.ts` is read once at module load and threaded through every `authOpts(username, CLIENT_ID)` call. No per-render hook — the id is stable for the tab session.

**`sseSubscriber.ts`** ([sseSubscriber.ts:32](../../../packages/HimalayaUI/frontend/src/lib/sseSubscriber.ts:32)):

```ts
// Replace:
if (event.actor && event.actor === ctx.username) return;
// With:
if (event.client_id && event.client_id === ctx.clientId) return;
```

`CurationEvent` interface gains `client_id?: string | null`. `ctx` shape changes from `{ username, qc }` to `{ clientId, qc }`. `clientId` is captured as a stable `const` at App component mount via `getClientId()` — no `useRef` needed (unlike `username`, `clientId` does not change after first load). The existing `usernameRef` in App.tsx, used today only by the SSE handler, can be removed entirely along with the `useRef` import.

`actor` stays on the SSE frame for downstream consumers that may want to display "edited by Alice" hints (none today, but the field is cheap and audit-aligned).

## Migration plan

Per the Plan 7 idempotency contract:

1. **Schema:** `ALTER TABLE user_actions ADD COLUMN client_id TEXT` appended to the existing `stmts` array in `migrate_schema!`. The surrounding try/catch loop handles the duplicate-column case on already-migrated DBs.
2. **No backfill.** Existing rows get NULL `client_id` and that's correct — they predate the concept. NULL also represents system-emitted events (`_system_request()` from `analyze_exposure!`); see "Behavior notes" below.
3. **No rollback path needed.** Additive nullable column; if the feature were reverted, the column sits unused. No data dependency on it from any other migration.

## Behavior notes

- **`_system_request()` events broadcast with NULL `client_id`** and pass the new filter — they refetch in every connected tab including the originator. This is **identical to today's behavior**: the current `actor`-based filter short-circuits on `event.actor && ...`, so events with `actor=null` (which `_system_request()` already produces) refetch in all tabs today. Not a regression.
- **User-triggered routes that internally call `_system_request()`** (e.g. `analyze_exposure!` emitting `analyze_run` sub-events) similarly broadcast with NULL `client_id` for those sub-events. The user-initiating route's primary event still carries the user's `client_id` and self-filters correctly. Sub-event redundant refetches are pre-existing.

## Test plan

### Backend

- **`test_actions.jl`:** new test for `get_client_id` (present, absent, empty-string).
- **`test_events.jl`:**
  - `apply_event!` writes `client_id` column when `X-Client-Id` header present.
  - `apply_event!` writes NULL `client_id` when header absent (existing test's behavior, with new assertion).
  - `broadcast_event!` includes `client_id` in the JSON frame.
  - `rebuild_views_from_log!` round-trip property is unaffected.
- **No new schema test needed** — existing `migrate_schema!` idempotency tests already exercise the additive-column pattern; this migration follows that template.

### Frontend

- **`api.test.ts`:** parallel assertions for `X-Client-Id` matching every existing `X-Username` mutation test (8 cases), plus a new clientId-without-username case.
- **`sse.test.tsx`:** existing file with 6 cases — update all to use new `{ clientId, qc }` ctx shape; add new cases asserting filter matches on `client_id`, ignores `actor` mismatches, accepts events with NULL `client_id` (system events), and accepts events with non-matching `client_id` (other tabs / other users).
- **`clientId.test.ts`:** mint-on-first-call, persist-across-calls, persist-across-`sessionStorage`-reload simulation.
- **No new Playwright E2E.** The behavior is observable via existing E2E tests *only* if we add a multi-tab scenario, which Playwright supports but is heavyweight. Skip for this spec; revisit if Compare's E2E suite naturally needs it.

### Manual verification

- Two tabs as `alice`, both on the same exposure: edits in tab A appear in tab B without manual reload.
- Two tabs as `alice`, tab A on Compare with exposure Y, tab B on Inspect for Y: peak edits in tab B refresh Y's panel in tab A.
- Single tab: editing still does *not* trigger a redundant SSE-driven refetch on top of the local mutation's `invalidateQueries`.

## Implementation plan

See [docs/superpowers/plans/2026-05-02-sse-client-id.md](../plans/2026-05-02-sse-client-id.md) — task-by-task with TDD steps, verified file paths, fixture patterns, and per-task verification commands. Sequential, single-PR; do not parallelize across tasks.

## Open questions

1. **Display `actor` in UI?** Out of scope for this spec; `actor` remains on the SSE frame for any future feature ("Bob is editing" indicators, presence dots) without recommitting to it as filter discriminator.
2. **Coordination with Compare design.** Compare's spec should reference this one as a precondition; this spec should ship before or alongside Compare so Compare's UX assumes correct behavior from day one.
3. **EventSource query-string fallback for server-side pruning?** Not in this spec. EventSource doesn't support custom request headers, so server-side filtering would need `client_id` on the URL: `GET /api/events?client_id=...`. Today's design keeps filtering purely client-side — every subscriber receives every event — which is simple and bandwidth-cheap at lab scale. If/when subscriber count grows or events get heavier (large payloads, frequent broadcasts), revisit. Compare-era follow-up at the earliest.

## Resolved decisions

These came up during plan iteration and are recorded here so the rationale isn't lost:

- **Module-level `getClientId()` over a `useClientId()` hook.** The id is constant for a tab's lifetime, so a hook adds nothing.
- **`clientId` not sent on GET requests.** Same rationale as `X-Username`: GETs are not events and don't appear in `user_actions`.
- **EventSource itself does not carry `client_id`.** EventSource (browser API) does not support custom headers, and filtering is purely client-side — every subscriber receives every event and filters locally. This keeps `broadcast_event!` simple.
- **`actor`-based filter is not retained as a fallback.** Once `client_id` ships, it's the sole filter mechanism. Two-mode complexity isn't worth the marginal robustness.
