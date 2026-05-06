# Compare Page Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Implement the Compare page per [docs/superpowers/specs/2026-05-02-compare-page-design.md](../specs/2026-05-02-compare-page-design.md): a draft/submit-style multi-trace overlay tool with author-only edit, fork-based collaboration, picker-modal selection, per-trace band heights, edit-mode peak click cycling, review-mode annotation toggles, chat thread, and `@comparison:N` mention support.

**Architecture:** Two view tables (`comparisons`, `comparison_members`) written exclusively by `update_view_for_event!`, plus `comparison_messages` written directly by the route handler (same pattern as `sample_messages` — the route does the INSERT, `apply_event!` is called for the event log + SSE broadcast only, and `events.jl` returns `nothing` for `post_message`). Three view-producing event kinds (`comparison_created`, `comparison_submitted`, `comparison_deleted`) routed through `apply_event!(InTransaction(), …)` inside `with_idempotency`. Plus the existing `post_message` kind extended for `entity_type='comparison_message'`. One frontend queue mutator (`saveComparison`) plus a thin `deleteComparison`. Page UI is structured as Sidebar + Plot Panel (plot column + metadata gutter) + Chat, mirroring the Index page's three-card composition.

**Tech stack:** Julia 1.9+, SQLite.jl, Oxygen.jl 1.10, JSON3, SHA (already a dep); React 18, TanStack Query v5, Zustand, Observable Plot, Vitest, Playwright.

**Read first:**
- [docs/superpowers/specs/2026-05-02-compare-page-design.md](../specs/2026-05-02-compare-page-design.md) — full design rationale, schema, every UX decision.
- [docs/event-log.md](../../event-log.md) — dispatcher contract, `apply_event!` semantics, hash invariants, SSE broadcast.
- [docs/mutation-queue.md](../../mutation-queue.md) — queue architecture, `client_op_id` lifecycle, deferred-promise pattern, `with_idempotency`, post-commit broadcast.
- [docs/contract-testing.md](../../contract-testing.md) — six-layer testing rule.
- [CLAUDE.md](../../../CLAUDE.md) — HimalayaUI gotchas (SQLite.jl `Tables.rowtable`/`missing` handling, FK-on-every-connection, AUTOINCREMENT for mention targets, `Manifest.toml` worktree gotcha).
- [docs/superpowers/plans/2026-04-29-chat-mention.md](2026-04-29-chat-mention.md) — pattern for adding a new mention-target type. The compare-mention work in Phase 10 mirrors this.

---

## Migration Architecture

> **This section governs all DB-touching work in this plan.** Each task's "Migration impact" subsection refers back to this.

### Deployment context

- **One central DB**, resolved by `default_db_path()` (env `HIMALAYA_DB_PATH` or `~/.himalaya/himalaya.db`).
- **`open_db` runs `create_schema!` then `migrate_schema!`** on every connection. Both must be idempotent across all DB states this plan introduces.
- **`Manifest.toml` is gitignored** so worktrees re-resolve. No new package deps in this plan.
- **DB shared across user sessions** — migrations may run with multiple readers connected. Each migration runs inside a single `SQLite.transaction`.

### Idempotency contract

Every migration step must be safe to run on:

1. **Fresh DB** — `create_schema!` produces the current schema.
2. **Pre-Compare DB** — Compare tables absent.
3. **Already-migrated DB** — second `open_db` is a no-op.

Concrete patterns:
- `CREATE TABLE IF NOT EXISTS` — idempotent.
- `CREATE INDEX IF NOT EXISTS` — idempotent.
- No backfills required (Compare adds new tables; nothing to transform).

### Ordering constraint

Compare's tables don't depend on each other beyond FK relationships, but they reference existing tables. Single migration function:

```julia
function migrate_compare!(db::SQLite.DB)
    # CREATE TABLE IF NOT EXISTS comparisons (...)
    # CREATE TABLE IF NOT EXISTS comparison_members (...)
    # CREATE TABLE IF NOT EXISTS comparison_messages (...)
    # CREATE INDEX IF NOT EXISTS idx_comparison_members_by_comparison
    # CREATE INDEX IF NOT EXISTS idx_events_by_comparison ... WHERE entity_type='comparison'
    # CREATE INDEX IF NOT EXISTS idx_comparison_messages_comparison
    # CREATE INDEX IF NOT EXISTS idx_comparisons_forked_from
end
```

Called from `migrate_schema!` after Plan 7's `migrate_r4_event_columns!` (so `user_actions.entity_type` is guaranteed to exist before we add the partial index keying on it).

### Rollback story

- **DB backup before any migration:** `cp ~/.himalaya/himalaya.db ~/.himalaya/himalaya.db.pre-compare-backup`.
- **No automated rollback.** Restore from backup if a migration fails mid-deploy.
- **Compare adds new tables only.** Failure leaves the DB without the new tables; existing functionality is unaffected.

### Pre-flight checklist

Before any task that touches `db.jl`:

- [ ] DB backup taken
- [ ] Worktree on the right branch
- [ ] `Pkg.test("HimalayaUI")` green before touching anything
- [ ] Frontend `npm test` and `npm run build` green
- [ ] [docs/event-log.md](../../event-log.md) and [docs/mutation-queue.md](../../mutation-queue.md) re-read recently

---

## E2E selector + accessibility strategy

> Required reading before any frontend phase. Selectors set here are referenced by tests in Phases 4–13.

Per CLAUDE.md frontend gotchas: **never assert against Tailwind class strings.** All test selectors use `data-testid`, `role`, or stable `data-*` attributes that survive styling changes. Each new affordance gets an explicit testid documented here so tests don't drift across phases.

| Affordance | Selector | Notes |
|---|---|---|
| Sidebar list item | `data-testid="comparison-list-item"` `data-comparison-id={id}` | One per row |
| Sidebar scope toggle | `data-testid="comparison-scope-toggle"` `data-scope={"this"\|"all"}` | |
| Sidebar new button | `data-testid="comparison-new"` | |
| Edit / Fork header button | `data-testid="comparison-edit"` / `data-testid="comparison-fork"` | Mutually exclusive based on author |
| Save / Cancel / Discard | `data-testid="comparison-save"` etc | |
| Needs-Review badge | `data-testid="comparison-needs-review"` | Visible only when stale |
| Lineage badge | `data-testid="comparison-lineage"` `data-parent-id={id\|"deleted"}` | |
| Forks popover trigger | `data-testid="comparison-forks-trigger"` | |
| Plot member trace | `data-testid="member-trace"` `data-member-id={id}` | One per band |
| Member metadata row | `data-testid="member-meta-row"` `data-member-id={id}` `data-stale={true\|false}` | |
| Band resize divider | `data-testid="band-divider"` `data-above-id={id}` `data-below-id={id}` | |
| Reorder grip | `data-testid="member-reorder-grip"` `data-member-id={id}` | |
| Peak triangle | `data-testid="peak-tick"` `data-member-id={id}` `data-peak-q={q}` `data-state={"shown"\|"labeled"\|"hidden"}` | Q in row to support hover; peak id is internal |
| Annotation toggles | `data-testid="annotation-toggle-{peaks\|labels}"` | |
| Grouping mode dropdown | `data-testid="grouping-mode"` `data-mode={"bySample"\|"byPhase"\|"distinct"}` | |
| Picker modal root | `data-testid="comparison-picker"` | |
| Picker filter chips | `data-testid="picker-filter-{experiment\|tag}"` | |
| Picker exposure row | `data-testid="picker-row"` `data-exposure-id={id}` `data-locked={true\|false}` | |
| Conflict modal root | `data-testid="conflict-modal"` | |
| Conflict modal actions | `data-testid="conflict-{discard\|overwrite\|fork}"` | |
| Mention chip (comparison) | Reuses existing MentionChip pattern with `data-mention-type="comparison"` `data-mention-id={id}` `data-hash-drift={true\|false}` | |

**`role` usage.** Modals have `role="dialog"` + `aria-labelledby`. Picker uses `role="listbox"`/`role="option"` for the rows. Toggles and checkboxes use native semantics (no extra `role` needed).

**Keyboard:** Esc closes any modal. Tab cycles focus inside the modal (via `useFocusTrap`). Picker rows are arrow-key navigable; Space toggles selection.

---

## File Map (consolidated)

| Action | Path | Phase |
|--------|------|-------|
| Modify | `packages/HimalayaUI/src/db.jl` | 1 (schema), 11 (auth helper) |
| Modify | `packages/HimalayaUI/src/events.jl` | 1 (dispatcher branches) |
| Create | `packages/HimalayaUI/src/comparisons.jl` | 2 (helpers: hash, fetch, lineage queries) |
| Create | `packages/HimalayaUI/src/routes_comparisons.jl` | 2 (REST routes) |
| Modify | `packages/HimalayaUI/src/server.jl` | 2 (route registration) |
| Modify | `packages/HimalayaUI/src/HimalayaUI.jl` | 2 (include new files) |
| Modify | `packages/HimalayaUI/src/json.jl` | 2 (comparison/member serialization) |
| Modify | `packages/HimalayaUI/src/routes_messages.jl` | 10 (comparison_messages support) |
| Modify | `packages/HimalayaUI/src/actions.jl` | 11 (`get_user_id_for_request` helper if missing) |
| Create | `packages/HimalayaUI/test/test_comparisons.jl` | 2 (route tests) |
| Modify | `packages/HimalayaUI/test/test_db.jl` | 1 (schema migration tests) |
| Modify | `packages/HimalayaUI/test/test_events.jl` | 1 (dispatcher round-trip tests) |
| Modify | `packages/HimalayaUI/test/test_idempotency_replay_invariant.jl` | 2 (comparison kinds in matrix) |
| Modify | `packages/HimalayaUI/test/test_route_response_shapes.jl` | 2 (comparison route shapes) |
| Modify | `packages/HimalayaUI/test/runtests.jl` | 2 (include new test file) |
| Create | `packages/HimalayaUI/frontend/src/pages/ComparePage.tsx` | 4 (review mode), 6 (plot), 7 (gutter), 8 (clicks), 9 (toggles) |
| Create | `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx` | 4 (edit mode shell) |
| Create | `packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx` | 4 |
| Create | `packages/HimalayaUI/frontend/src/components/ComparisonPicker.tsx` | 5 |
| Create | `packages/HimalayaUI/frontend/src/components/MultiTracePlot.tsx` | 6 |
| Create | `packages/HimalayaUI/frontend/src/components/MemberTraceLayer.tsx` | 6 |
| Create | `packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx` | 7 |
| Create | `packages/HimalayaUI/frontend/src/components/BandResizeDivider.tsx` | 7 |
| Create | `packages/HimalayaUI/frontend/src/components/AnnotationToggles.tsx` | 9 |
| Create | `packages/HimalayaUI/frontend/src/components/LineageBadge.tsx` | 11 |
| Create | `packages/HimalayaUI/frontend/src/components/ConflictModal.tsx` | 12 |
| Create | `packages/HimalayaUI/frontend/src/components/ExposureListRow.tsx` | 5 (extracted, reused on Inspect) |
| Modify | `packages/HimalayaUI/frontend/src/pages/InspectPage.tsx` | 5 (warm-add menu), uses `ExposureListRow` |
| Modify | `packages/HimalayaUI/frontend/src/components/ChatCard.tsx` | 10 (accept generic entity scope) |
| Modify | `packages/HimalayaUI/frontend/src/components/MentionPicker.tsx` | 10 (Comparisons group) |
| Modify | `packages/HimalayaUI/frontend/src/components/MentionChip.tsx` | 10 (comparison rendering + hash drift indicator) |
| Modify | `packages/HimalayaUI/frontend/src/hooks/useMentionResolution.ts` | 10 (comparison entity type) |
| Modify | `packages/HimalayaUI/frontend/src/api.ts` | 4, 5, 11, 10 |
| Modify | `packages/HimalayaUI/frontend/src/queries.ts` | 4, 5, 10 |
| Modify | `packages/HimalayaUI/frontend/src/state.ts` | 3 (draft slot, mode, toggles) |
| Modify | `packages/HimalayaUI/frontend/src/App.tsx` | 4 (route registration), 12 (conflict modal mount) |
| Create | `packages/HimalayaUI/frontend/src/lib/queue/mutators/saveComparison.ts` | 3 |
| Create | `packages/HimalayaUI/frontend/src/lib/queue/mutators/deleteComparison.ts` | 3 |
| Modify | `packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts` | 3 (comparison kinds) |
| Modify | `packages/HimalayaUI/frontend/src/lib/queue/types.ts` | 3 (extend `OpKind` union with `"comparison_save" \| "comparison_delete"`) |
| Modify | `packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts` | 3 (register both) |
| Create | `packages/HimalayaUI/frontend/src/lib/comparison/draft.ts` | 4 (sessionStorage shape + hooks) |
| Create | `packages/HimalayaUI/frontend/src/lib/comparison/coloring.ts` | 9 (palette / grouping mode) |
| Create | `packages/HimalayaUI/frontend/src/lib/comparison/normalization.ts` | 6 (peak-fit/signal-fit logic) |
| Create | `packages/HimalayaUI/frontend/src/lib/comparison/snapshot.ts` | 3 (client-side `computeMemberSnapshot` — reads TanStack cache, produces `MemberSnapshot` shape) |
| Create | `packages/HimalayaUI/frontend/src/lib/comparison/contentHash.ts` | 3 (canonical client-side hash for citation) |
| Modify | `packages/HimalayaUI/frontend/src/main.tsx` | 4 (router config) |
| Modify | `packages/HimalayaUI/frontend/src/phases.ts` | (no change expected; reuse) |
| Create | `packages/HimalayaUI/frontend/test/queue/saveComparison.test.tsx` | 3 |
| Create | `packages/HimalayaUI/frontend/test/queue/deleteComparison.test.tsx` | 3 |
| Modify | `packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts` | 3 (comparison kinds) |
| Modify | `packages/HimalayaUI/frontend/test/queue/sseEventPayload.contract.test.ts` | 3 |
| Modify | `packages/HimalayaUI/frontend/test/queue/rollbackSymmetry.test.ts` | 3 |
| Modify | `packages/HimalayaUI/frontend/test/queue/authHeaders.test.ts` | 3 |
| Modify | `packages/HimalayaUI/frontend/test/queue/mutatorOnSseWins.test.ts` | 3 |
| Create | `packages/HimalayaUI/frontend/test/MultiTracePlot.test.tsx` | 6 |
| Create | `packages/HimalayaUI/frontend/test/MemberMetaRow.test.tsx` | 7 |
| Create | `packages/HimalayaUI/frontend/test/BandResizeDivider.test.tsx` | 7 |
| Create | `packages/HimalayaUI/frontend/test/ComparisonPicker.test.tsx` | 5 |
| Create | `packages/HimalayaUI/frontend/test/ConflictModal.test.tsx` | 12 |
| Create | `packages/HimalayaUI/frontend/test/draftPersistence.test.ts` | 4 |
| Create | `packages/HimalayaUI/frontend/test/coloring.test.ts` | 9 |
| Create | `packages/HimalayaUI/frontend/test/normalization.test.ts` | 6 |
| Create | `packages/HimalayaUI/frontend/test/contentHash.test.ts` | 3 |
| Create | `packages/HimalayaUI/frontend/e2e/compare.spec.ts` | 13 (smoke) |

---

## Phase 1: Schema + dispatcher

### Task 1.1: Schema migration

**Files:**
- Modify: `packages/HimalayaUI/src/db.jl`
- Modify: `packages/HimalayaUI/test/test_db.jl`

**Migration impact:** Adds three new tables (`comparisons`, `comparison_members` with `snapshot TEXT NOT NULL CHECK (json_valid(snapshot))` and `peak_display TEXT CHECK (peak_display IS NULL OR json_valid(peak_display))`, `comparison_messages`) and four indexes. Pure additions; no backfill.

- [ ] **Step 1: Write the failing migration tests** in `test_db.jl`. Cover:
  - **Fresh DB:** after `open_db`, all three Compare tables exist with correct columns; all four indexes exist.
  - **Already-migrated DB:** second `open_db` is a no-op (no errors, no duplicate tables).
  - **AUTOINCREMENT contract:** `comparisons.id` strict-monotonic — deleting and re-inserting does not reuse the freed id (mention-target stability rule). `comparison_members.id` and `comparison_messages.id` use plain `INTEGER PRIMARY KEY` (neither is `@`-mentioned; `comparison_messages` matches `sample_messages` symmetry).
  - **FK enforcement:** inserting `comparison_members` with `comparison_id` referencing a non-existent comparison fails when `PRAGMA foreign_keys = ON`.
  - **`ON DELETE SET NULL` on `comparison_members.exposure_id`:** delete an exposure that's a member; member row survives with `exposure_id IS NULL`.
  - **`ON DELETE CASCADE` on `comparison_members.comparison_id` and `comparison_messages.comparison_id`:** delete a comparison; both child tables drop their rows.
  - **`ON DELETE SET NULL` on `forked_from_id`:** delete a parent comparison; fork survives with `forked_from_id IS NULL`.

- [ ] **Step 2: Implement `migrate_compare!`** in `db.jl`. Add three `CREATE TABLE IF NOT EXISTS` statements and four `CREATE INDEX IF NOT EXISTS` statements. Wire into `migrate_schema!` *after* the Plan 7 R4 migration (so `user_actions.entity_type` exists before we partial-index on it).

- [ ] **Step 3: Verify `Pkg.test("HimalayaUI")` passes.**

- [ ] **Step 4: Commit** as `feat(compare): schema migration for comparisons / comparison_members / comparison_messages`.

### Task 1.2: Dispatcher branches in `events.jl`

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl`
- Modify: `packages/HimalayaUI/test/test_events.jl`

**Migration impact:** Code-only. Three new branches in `update_view_for_event!`; one new helper for the snapshot diff.

The dispatcher writes to `comparisons` and `comparison_members` exclusively. Per [event-log.md](../../event-log.md) §1, no route or pipeline function may write to these tables directly.

- [ ] **Step 1: Round-trip test in `test_events.jl`** for each of the three view-producing kinds. The property: `apply_event!` writes vs `rebuild_views_from_log!`-fold produce identical state. Test cases:
  - `comparison_created` with 0 members → table rows match payload.
  - `comparison_created` with 3 members → all members inserted; `display_order` matches payload order; `created_by`/`created_at` populated.
  - `comparison_submitted` no-op (same snapshot) → no row changes; `content_hash` may stay or recompute identically.
  - `comparison_submitted` member added → previous + new member both present.
  - `comparison_submitted` member removed → previous minus removed.
  - `comparison_submitted` member updated (same id, new field) → row updated in place.
  - `comparison_submitted` reorder → all `display_order` values match payload.
  - `comparison_deleted` → cascades to `comparison_members` and `comparison_messages`.

- [ ] **Step 2: Implement** the three dispatcher branches in `events.jl`. The `comparison_submitted` branch follows the diff-based pattern documented in the spec (`update_view_for_comparison_submitted!`) — note the **three** member dispositions: DELETE (in DB, not in payload), UPDATE (in both, payload has integer `id`), INSERT (payload `id === nothing`). **Snapshot data comes from the payload, not from the DB** — each member in the payload carries its own `snapshot` JSON, computed by the client at submit time against the cache state the user was looking at. The dispatcher reads `m.snapshot` verbatim and writes it to the `snapshot` column. **Every UPDATE writes the snapshot unconditionally** — not conditionally on whether `exposure_id` changed (the user may have confirmed a new index on the same exposure between submissions). Error on zero-row UPDATE (payload references a member id that doesn't exist for this comparison). Use `Tables.rowtable` for any read-then-modify pattern (per CLAUDE.md SQLite gotcha). Be careful with the JSON `null` ↔ Julia `nothing` ↔ `Tables.rowtable` `missing` triangle: payload-side member ids are typed `Union{Nothing,Int}` (never `Missing`); DB-side reads use `ismissing()` to detect SQL NULL.

- [ ] **Step 3: Implement `compute_content_hash`** helper (will live in `comparisons.jl` once Task 2.1 lands; for now, inline in `events.jl` and refactor later, OR scaffold `comparisons.jl` here). Recommendation: scaffold `comparisons.jl` now, since it's the same effort.

- [ ] **Step 4: Run round-trip tests + idempotency-replay tests.** All green.

- [ ] **Step 5: Commit** as `feat(compare): dispatcher branches + content_hash`.

### Task 1.3: `rebuild_views_from_log!` extension

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl`
- Modify: `packages/HimalayaUI/test/test_events.jl`

`rebuild_views_from_log!(db, exposure_id)` exists today and folds `entity_type='exposure'` events. The Compare equivalent — `rebuild_views_from_log!(db, comparison_id; entity_type="comparison")` — folds `entity_type='comparison'` events.

- [ ] **Step 1: Test** that wiping `comparison_members` for one comparison and replaying its event log reconstructs identical state.

- [ ] **Step 2: Implement** by parameterizing the existing function or adding a sibling. Choose whichever fits the existing code shape best — discuss in PR if unclear.

- [ ] **Step 3: Commit** as `feat(compare): rebuild_views_from_log! supports comparison entity type`.

---

## Phase 2: Backend routes

### Task 2.1: `comparisons.jl` helpers

**Files:**
- Create: `packages/HimalayaUI/src/comparisons.jl`
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl` (include)
- Create: `packages/HimalayaUI/test/test_comparisons.jl`
- Modify: `packages/HimalayaUI/test/runtests.jl` (include)
- Modify: `packages/HimalayaUI/src/actions.jl` (or wherever user-id resolution lives) — add `get_user_id_for_request(req)::Union{Int,Nothing}` if absent

Pure helpers (no HTTP):

- `compute_content_hash(db, comparison_id)` — SHA-256 over canonical serialization (title, description, ordered members tuple per spec, including snapshot JSON).
- `current_content_hash(db, comparison_id)` — fetches stored `content_hash` (returns `nothing` if comparison doesn't exist).
- `compute_member_snapshot(db, exposure_id)` — returns the JSON shape per spec ("Derived analysis state and staleness"): `effective_peaks` + `confirmed_index` (R²-gated) + `analysis_inputs_hash`. **Note:** the existing `effective_peaks(db, …)` helper returns `(q, sharpness, peak_id, peak_kind)` — no `intensity`. The snapshot shape requires intensity, so `compute_member_snapshot` must join the peak ids against `auto_peaks.intensity` (for auto peaks). **Manual peaks carry `intensity: null`** — `peak_curations` has no intensity column, and `get_peaks_for_exposure` returns `NULL AS intensity` for them (see `routes_peaks.jl:86`). This matches the existing API behavior. The normalization layer must handle null intensity gracefully: manual peaks are excluded from peak-fit reference computation (fall back to signal-fit if all peaks in the q_window are manual). The R² gate threshold (0.98) should be extracted as a named constant shared with `PhasePanel` to prevent drift.
- `is_member_stale(db, member)::Bool` — compares `member.snapshot.analysis_inputs_hash` with `current_analysis_inputs_hash(member.exposure_id)`.
- `fetch_comparison_with_members(db, comparison_id)` — returns the full comparison + members shape used in API responses. Handles `exposure_id IS NULL` orphan placeholders. Includes per-member `is_stale: Bool` flag.
- `member_ids_for_comparison(db, comparison_id)` — returns `Set{Int}` of current member ids.
- `recently_used_exposures(db, user_id; limit=20)` — runs the picker history query.
- `comparisons_for_experiment(db, experiment_id)` — derived listing (joins through samples → exposures → members; sorted by `MAX(user_actions.created_at)`).
- `forks_of_comparison(db, comparison_id)` — list comparisons with `forked_from_id = :id`.
- **`is_author(db, comparison_id, user_id)::Bool`** — used by route gates in Task 2.2. Returns `false` if `comparison_id` doesn't exist OR if `created_by IS NULL` OR if `created_by != user_id`. (Moved here from Phase 11; the route gates in Task 2.2 depend on it.)

- [ ] **Step 1:** Write tests covering each helper with a sample DB.
- [ ] **Step 2:** Implement. Use `Tables.rowtable` consistently. Be careful with `ismissing` vs `nothing` per CLAUDE.md.
- [ ] **Step 3:** Commit as `feat(compare): query helpers in comparisons.jl`.

### Task 2.2: `routes_comparisons.jl`

**Files:**
- Create: `packages/HimalayaUI/src/routes_comparisons.jl`
- Modify: `packages/HimalayaUI/src/server.jl` (register routes + include)
- Modify: `packages/HimalayaUI/src/HimalayaUI.jl`
- Modify: `packages/HimalayaUI/src/json.jl` (serializers)

Routes per spec REST API section. Every mutating route wraps in `with_idempotency(db, req)` and uses `apply_event!(InTransaction(), …)`.

Routes:
- `GET /api/experiments/:eid/comparisons` — derived listing.
- `GET /api/comparisons` — global listing.
- `POST /api/comparisons` — create + initial submit. Body must include ≥1 member; otherwise 400. Optional `forked_from_id` + `forked_at_hash`.
- `GET /api/comparisons/:id` — full comparison with members + lineage info.
- `GET /api/comparisons/:id/forks` — list forks.
- `POST /api/comparisons/:id/submit` — author-only (403 if not). 409 on `expected_content_hash` mismatch. Idempotent.
- `DELETE /api/comparisons/:id` — author-only.
- `GET /api/comparisons/:id/messages` — chat thread.
- `POST /api/comparisons/:id/messages` — `post_message` event with `entity_type='comparison_message'`.

- [ ] **Step 1: Write** route-shape tests (`test_route_response_shapes.jl` rows + per-route happy/sad paths in `test_comparisons.jl`).
  - Create-with-zero-members → 400.
  - Create-with-fork-payload → response includes `forked_from_id` + `forked_at_hash`.
  - Submit-as-non-author → 403.
  - Submit-with-wrong-`expected_content_hash` → 409 with `current_hash` + `current_state` body.
  - Submit-with-correct-hash → 200 + new state in response.
  - Delete-as-non-author → 403.
  - Idempotent retry of *successful* submit (status 200) → cached response, no duplicate event. (Status ≥ 400 retries re-evaluate; see 409 retry contract below.)
  - **409 retry contract:** per `with_idempotency`, status ≥ 400 responses are NOT cached. A retry of a submit that originally returned 409 (under same `client_op_id`) re-evaluates the conflict check (correct: the conflict may have resolved between retries). The client breaks out of a stale conflict by minting a fresh `client_op_id` (which the conflict-modal flow does naturally because each modal action issues a new `useQueueMutation.mutate()` call).
  - Cascade test: delete a comparison → children gone (members + messages).

- [ ] **Step 2: Implement** routes. Per [mutation-queue.md](../../mutation-queue.md) §6: use `apply_event!(InTransaction(), …)` and let the post-commit broadcast queue handle SSE.

- [ ] **Step 3: Idempotency-replay matrix** — add comparison kinds to `test_idempotency_replay_invariant.jl` so SSE-fanout-under-retry is verified. Include a **409-retry-re-evaluates row**: submit with wrong `expected_content_hash` → 409; retry with same `client_op_id` → re-evaluates conflict check (per `with_idempotency` contract: status ≥ 400 not cached). Also test that a successful submit (status 200) IS cached on retry.

- [ ] **Step 4: Commit** as `feat(compare): REST routes with idempotency + auth gates`.

### Task 2.3: `post_message` extension for comparisons

**Files:**
- Modify: `packages/HimalayaUI/src/events.jl` (`post_message` dispatch by `entity_type`)
- Modify: `packages/HimalayaUI/src/routes_messages.jl` (or `routes_comparisons.jl` depending on organization)

Per spec: the existing `post_message` kind (used with `entity_type='sample_message'` in `routes_messages.jl:64`) gains `entity_type='comparison_message'` semantics. **The existing handler hard-codes `sample_messages` — this is a new dispatch branch, not a trivial extension.** The INSERT target table must be chosen based on `entity_type`.

- [ ] **Step 1:** Test that `post_message` with `entity_type='comparison_message'` writes to `comparison_messages`, and with `entity_type='sample_message'` still writes to `sample_messages`. Also test that an unknown `entity_type` returns an error (not a silent no-op).
- [ ] **Step 2:** Implement the dispatch. The route handler at `POST /api/comparisons/:id/messages` does the INSERT into `comparison_messages` directly (same pattern as `sample_messages` — see `routes_messages.jl:53`), then calls `apply_event!` for the event log + SSE broadcast only. Wrap in `with_idempotency(db, req) do ... end` per the standard discipline.
- [ ] **Step 2b:** Update `applyRemoteToCache.ts`'s existing `post_message` case. The current implementation hard-codes `payload.sample_id` and `SampleMessage` type. Add an `entity_type` dispatch: `entity_type === 'sample_message'` → existing `queryKeys.messages(sampleId)` path; `entity_type === 'comparison_message'` → new `queryKeys.comparisonMessages(comparisonId)` path. Without this, comparison chat SSE events write to the wrong query key.
- [ ] **Step 3:** Commit as `feat(compare): post_message event routes by entity_type`.

---

## Phase 3: Frontend queue mutators

### Task 3.1: `saveComparison` mutator

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/queue/mutators/saveComparison.ts`
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/mutatorRegistry.ts`
- Modify: `packages/HimalayaUI/frontend/src/lib/queue/applyRemoteToCache.ts`
- Modify: `packages/HimalayaUI/frontend/src/api.ts`
- Create: `packages/HimalayaUI/frontend/test/queue/saveComparison.test.tsx`
- Modify: `packages/HimalayaUI/frontend/test/queue/cache-shape.test.ts`
- Modify: `packages/HimalayaUI/frontend/test/queue/sseEventPayload.contract.test.ts`
- Modify: `packages/HimalayaUI/frontend/test/queue/rollbackSymmetry.test.ts`
- Modify: `packages/HimalayaUI/frontend/test/queue/authHeaders.test.ts`
- Modify: `packages/HimalayaUI/frontend/test/queue/mutatorOnSseWins.test.ts`

Per [contract-testing.md](../../contract-testing.md), every new `OpKind` adds a row to all six contract tests. Note that the **2 OpKinds emit 3 distinct SSE event kinds** (`saveComparison` emits `comparison_created` for create or `comparison_submitted` for update; `deleteComparison` emits `comparison_deleted`). The OpKind-shaped contract tests (`saveComparison.test.tsx`, `deleteComparison.test.tsx`, `authHeaders.test.ts`, `rollbackSymmetry.test.ts`) get **2 rows each** (one per OpKind). The event-shape-shaped contract tests (`cache-shape.test.ts`, `sseEventPayload.contract.test.ts`, `mutatorOnSseWins.test.ts`) get **3 rows each** (one per event kind) — otherwise the `comparison_deleted` `applyRemoteToCache` branch ships untested. Total cases: 2×4 + 3×3 = **17**, not 12.

The mutator (`kind: "comparison_save"`) handles both create and update — `payload.id` absent means create; present means submit. `request` calls `POST /api/comparisons` or `POST /api/comparisons/:id/submit` accordingly. `onMutate` is a no-op (no optimistic write — submission shows a spinner). `onSuccess` writes the canonical comparison + members into the cache and invalidates listing keys.

- [ ] **Step 1: Write** the mutator skeleton + the six contract test rows.
- [ ] **Step 2: Implement** the mutator. Pay attention to:
  - `payload.id` discriminates create vs submit.
  - `expected_content_hash` is included on submits (not creates).
  - 409 response should be detected in `request` and thrown as a typed `ConflictError` so it lands in `onError` and the UI can open the conflict modal. See [mutation-queue.md](../../mutation-queue.md) §10.
  - 403 response is a validation-class error (no retry).

- [ ] **Step 3: `applyRemoteToCache` branches** for `comparison_created` / `comparison_submitted` / `comparison_deleted`. Foreign `comparison_created` and `comparison_submitted` events invalidate `comparison(id)`, `comparisonMembers(id)`, `comparisons` listing. Foreign `comparison_deleted` events use `qc.removeQueries` (NOT `invalidateQueries` — refetching a deleted resource 404s and leaves stale error state) for `comparison(id)` and `comparisonMembers(id)`, plus `qc.setQueryData` to filter the id out of listing keys. **Note:** verify that `replayCoordinator.handleRemoteEvent` forwards comparison events correctly — the coordinator dispatches on `remote.kind` via `applyRemoteToCache`'s switch statement with no `entity_type` filter, so adding new `case` branches is sufficient. Also extend `OpKind` union in `types.ts`.

- [ ] **Step 4: Commit** as `feat(compare): saveComparison queue mutator + contract tests`.

### Task 3.2: `deleteComparison` mutator

Similar shape, fewer edges. `kind: "comparison_delete"`. Same six contract tests.

- [ ] **Steps 1–4:** Mirror Task 3.1 but for `DELETE /api/comparisons/:id`.

### Task 3.3: Client-side `computeMemberSnapshot` + `MemberSnapshot` type

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/comparison/snapshot.ts`
- Modify: `packages/HimalayaUI/frontend/src/api.ts` (`MemberSnapshot` type defined here — shared by HTTP response parser and SSE handler)

`computeMemberSnapshot(exposureId, qc)` reads from the TanStack cache: `effective_peaks` from `queryKeys.peaks(exposureId)`, `confirmed_index` from `queryKeys.indices(exposureId)` (R²-gated — below 0.98 produces `null`), and `analysis_inputs_hash` from the exposure query. Returns the `MemberSnapshot` shape defined in the spec. Called by the save flow at submit time for every member.

The `MemberSnapshot` type is defined in `api.ts` (not in `snapshot.ts`) because it must be the single source of truth for both the HTTP response parser and the SSE `applyRemoteToCache` handler — both paths must produce the same parsed shape to prevent cache divergence during reconciliation.

- [ ] **Step 1: Write** tests asserting that `computeMemberSnapshot` produces the correct shape for a mock cache state, including the R²-gate (below 0.98 → `confirmed_index: null`).
- [ ] **Step 2: Implement.**
- [ ] **Step 3:** Commit as `feat(compare): client-side computeMemberSnapshot`.

### Task 3.4: Client-side `contentHash`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/comparison/contentHash.ts`
- Create: `packages/HimalayaUI/frontend/test/contentHash.test.ts`

Mirrors the server's canonical serialization so the client can verify on render and citation. SHA-256 in browser via `crypto.subtle.digest`.

- [ ] **Step 1: Write** tests asserting that for a fixed input shape, the JS hash matches the Julia hash. Use a small fixture committed to the test directory, hashed by both sides. This fixture also pins the client-side `computeMemberSnapshot` output against the server-side `compute_member_snapshot` output for the same input.
- [ ] **Step 2: Implement.**
- [ ] **Step 3:** Commit as `feat(compare): client-side content_hash for citation parity`.

---

## Phase 4: Frontend list page + create flow

### Task 4.1: Routes + page shells

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/main.tsx` (router)
- Create: `packages/HimalayaUI/frontend/src/pages/ComparePage.tsx`
- Create: `packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx`

Routes:
- `/experiments/:eid/compare` — sidebar list scoped to experiment, no comparison selected (empty state with "+ New" prompt).
- `/experiments/:eid/compare/:id` — review mode of comparison `:id`, sidebar scoped to experiment.
- `/experiments/:eid/compare/:id/edit` — edit mode.
- `/compare/all` — global listing scope (all comparisons across experiments).

`:eid` param feeds the picker default filter and the sidebar scope toggle's "this experiment" option. `/compare/all` skips the experiment context.

- [ ] **Step 1:** Wire routes; render placeholder content.
- [ ] **Step 2:** Add `queryKeys.comparisons(experimentId)`, `queryKeys.comparison(id)`, `queryKeys.comparisonMembers(id)`, `queryKeys.comparisonForks(id)` in `queries.ts`.
- [ ] **Step 3:** Add `getComparisons / getComparison / getComparisonForks` fetchers in `api.ts`.
- [ ] **Step 4:** Commit as `feat(compare): page shells + queries`.

### Task 4.2: Sidebar with scope toggle

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ComparisonSidebar.tsx`

- Comparison list, sorted by latest event.
- Header toggle: "This experiment" (default if `:eid` present) / "All experiments".
- Search box, pin toggle (per-user preference; could defer pin to Phase 13).
- "+ New" button → navigates to `/experiments/:eid/compare/new` which mounts `ComparePageEdit` with empty draft.
- Active comparison highlighted.

- [ ] **Step 1:** Tests with a small comparison list mock; assert sort order, scope toggle behavior, search filtering.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): sidebar with scope toggle`.

### Task 4.3: Draft state + sessionStorage

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/comparison/draft.ts`
- Modify: `packages/HimalayaUI/frontend/src/state.ts`
- Create: `packages/HimalayaUI/frontend/test/draftPersistence.test.ts`

Zustand slot:

```ts
type ActiveDraft = {
  id: number | undefined;        // not optional — TS strict + exactOptionalPropertyTypes
  baseHash: string | undefined;  // server's content_hash at edit-mode entry; frozen until submit
  title: string;
  description: string;
  members: DraftMember[];
};
type ActiveDraftSlot = ActiveDraft | null;
```

The `id: number | undefined` form (rather than `id?: number`) is required by `exactOptionalPropertyTypes: true` — see CLAUDE.md frontend gotchas. Same rule for nullable fields on `DraftMember`:

```ts
type DraftMember = {
  id: number | undefined;           // undefined = new member (INSERT on submit)
  exposure_id: number;
  display_order: number;
  band_height: number;
  y_offset: number;
  normalization: "none" | "max" | "area" | "qwindow";
  color_override: string | undefined;
  label_override: string | undefined;
  q_window_min: number | undefined;
  q_window_max: number | undefined;
  peak_display: { hidden: number[]; labeled: number[] } | undefined;
  snapshot: MemberSnapshot | undefined; // undefined during editing; computed fresh at submit time via computeMemberSnapshot
};
```

`baseHash` is captured at edit-mode entry (when loading an existing comparison's state into the draft slot) and **does not refresh** if a foreign SSE event arrives mid-edit — submission compares this baseline against the server's current hash via `expected_content_hash`, and a mismatch is a real conflict. Auto-updating `baseHash` from SSE events would silently swallow conflicts.

**Create-flow lifecycle:** for a brand-new comparison (never submitted), `id` is `undefined` and `baseHash` is `undefined`. The submit payload uses `POST /api/comparisons` (not `/submit`) and omits `expected_content_hash` entirely — there's no prior server state to conflict with, so the 409 path is unreachable on first submit. After the first submit succeeds and the comparison gets a server-assigned `id`, subsequent edits load with `baseHash` set to the server's `content_hash`.

Mirrored to `sessionStorage` with a schema version. On reload, rehydrate. Closing the tab loses the draft (acceptable for v1 per spec).

- Named actions: `startNewDraft`, `loadDraftFromComparison`, `setDraftTitle`, `setDraftDescription`, `addMember`, `removeMember`, `updateMember`, `reorderMembers`, `resizeBands`, `discardDraft`.
- **`loadDraftFromComparison` recovery path for stale comparisons:** loads saved render parameters (offsets, normalization, colors, labels, q-windows, peak_display) from the server state, but calls `computeMemberSnapshot` against the **current** TanStack cache for each member — not the saved snapshot. This means re-submit = "refresh to current truth." Rationale: edit mode shows live trace data; showing stale peaks that don't match the live trace would be confusing.

- [ ] **Step 1:** Tests for round-trip persistence, schema-version mismatch handling, action correctness.
- [ ] **Step 2:** Implement. Per CLAUDE.md frontend gotchas, no raw `setState`.
- [ ] **Step 3:** Commit as `feat(compare): draft state + sessionStorage persistence`.

### Task 4.4: Edit mode shell

`ComparePageEdit` reads draft from Zustand, renders title/description inputs, plot placeholder, member panel placeholder, "Save" / "Cancel" / "Discard draft" buttons. Save calls `saveComparison` queue mutator. Cancel returns to `/experiments/:eid/compare/:id` (or list page if creating). Discard clears draft.

- [ ] **Step 1:** Tests for save/cancel/discard buttons, blocking save when 0 members.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): edit mode shell + save flow`.

---

## Phase 5: Picker modal

### Task 5.1: `ExposureListRow` extraction

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ExposureListRow.tsx`
- Modify: `packages/HimalayaUI/frontend/src/pages/InspectPage.tsx` (consume)

Single row: exposure name, sample name, sample notes (truncated, full on hover), checkbox or click handler. Used by the picker and the Inspect page's exposure list.

- [ ] **Step 1:** Tests + implement.
- [ ] **Step 2:** Refactor Inspect's exposure list to use it. Verify existing Inspect tests still pass after the extraction — the component interface must match what InspectPage currently renders.
- [ ] **Step 3:** Commit as `refactor: extract ExposureListRow`.

### Task 5.2: Picker modal

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ComparisonPicker.tsx`
- Create: `packages/HimalayaUI/frontend/test/ComparisonPicker.test.tsx`

Search + experiment filter + tag filter + recently-used + already-added locks + sort by exposure name, all per spec.

- Backend route for picker recently-used: `GET /api/users/:id/recently-picked-exposures?limit=20`. Add this in Phase 5 (small backend addition).
- Backend route for tag listing in scope: `GET /api/experiments/:eid/sample-tags` returning distinct `(key, value)` pairs. Add here.

- [ ] **Step 1:** Backend routes + tests. The tag listing route must return distinct `(key, value)` pairs (not just flat values) — add a test that seeds two tags with the same value but different keys and asserts both appear as separate filter options.
- [ ] **Step 2:** Frontend tests covering filters, sort, multi-select, already-added locks, empty state, focus-trap (Tab cycles within modal, Esc closes), restored focus on close.
- [ ] **Step 3:** Implement frontend modal. Wrap the modal container with `useFocusTrap(containerRef, isOpen)` per the existing NavModal/OnboardingFlow pattern.
- [ ] **Step 4:** Wire to edit mode's "Add traces" button.
- [ ] **Step 5:** Commit as `feat(compare): picker modal + backend support routes`.

### Task 5.3: Warm-add accelerator on Inspect

Inspect page gains an "Add to comparison" affordance. Dropdown options:
- "Recent comparison: <title>" (the user's most recent draft, if any in sessionStorage)
- "Pick a comparison..." → opens picker modal scoped to comparisons
- "+ New comparison" → navigates to edit mode of new draft pre-populated with the selected exposure

- [ ] **Step 1:** Tests covering each menu path. Include a **cross-tab boundary test**: adding to "current draft" from Inspect only affects the same tab's Zustand draft — a second tab's draft is unchanged until the user submits and the SSE event fires. (Drafts are `sessionStorage`-scoped per tab; no `BroadcastChannel` in v1.)
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): warm-add from Inspect page`.

---

## Phase 6: Multi-trace plot

### Task 6.1: `MemberTraceLayer`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/MemberTraceLayer.tsx`
- Create: `packages/HimalayaUI/frontend/src/lib/comparison/normalization.ts`
- Create: `packages/HimalayaUI/frontend/test/normalization.test.ts`

Renders one member's marks (line + peak ticks + labels) into a y-band. Receives `member` (which carries `member.snapshot.effective_peaks` and `member.snapshot.confirmed_index`), `yBand: [number, number]`, `peakDisplay`, `onPeakClick?`, `highlightedIndexId?` (for hover-driven phase coloring). Computes peak-fit/signal-fit reference value adaptively against `snapshot.effective_peaks` (peak max in qwindow if peaks exist, else signal max in qwindow). Maps trace y values into `working_band` (inner 70% of `yBand`); clips to total band.

Peak rendering uses `snapshot.effective_peaks` (NOT a live peaks query), rendered in **black** by default. When `highlightedIndexId` is set and matches `snapshot.confirmed_index.id`, peaks belonging to that index render in the phase color from `phases.ts` instead. Non-index peaks stay black regardless.

The trace `(q, I)` data itself is read live from the exposure file via the existing `useTrace(exposureId)` query — only derived analysis state (peaks, index) comes from the snapshot.

Normalization library:

- `computeReference(trace, peaks, qWindow)` → number (the reference value)
- `applyNormalization(trace, reference, normalization, yBand, workingBandFraction=0.7)` → array of `{q, y}` points in plot coordinates.

- [ ] **Step 1:** Test normalization: peak-bearing trace normalizes on peak max; peakless trace falls back to signal max; tail clipping at total band envelope; working band fraction respected.
- [ ] **Step 2:** Test `MemberTraceLayer` rendering: line points correct, peak ticks at q positions, labels for `peak_display.labeled` only.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Commit as `feat(compare): MemberTraceLayer + adaptive normalization`.

### Task 6.2: `MultiTracePlot`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/MultiTracePlot.tsx`
- Create: `packages/HimalayaUI/frontend/test/MultiTracePlot.test.tsx`

Single Observable `<Plot>` with shared q-scale + zoom domain. Computes y-bands from `band_height` ratios (per spec formula). Composes N `MemberTraceLayer` marks. Brush + double-click reset zoom; mouse-wheel pan/zoom. Same control surface as existing `TraceViewer` — reuse hooks where possible. **Note:** brush-to-zoom requires the Observable Plot `.scale(name).invert(px)` cast pattern documented in CLAUDE.md frontend gotchas (`(el as unknown as { scale: ... })`) — extract a shared `invertQ(plot, px)` helper for both `MultiTracePlot` and `TraceViewer`.

- Aspect ratio: `0.3` hardcoded constant exported as `COMPARE_PLOT_ASPECT` for future surfacing.
- `panel_height` derived from container; `MultiTracePlot` is `flex-1` and queries its own bounding rect via `useElementSize` or similar.

- [ ] **Step 1:** Test y-band layout from `band_height` ratios; reorder shifts bands correctly; zoom behavior; double-click reset.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): MultiTracePlot composition`.

---

## Phase 7: Metadata gutter + reorder + resize

### Task 7.1: `MemberMetaRow`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx`
- Create: `packages/HimalayaUI/frontend/test/MemberMetaRow.test.tsx`

Single line per trace in review mode (label / phase chip / `a` / R² / `K` for cubics). Hover or click expands into a detail card overlay (peak count + other secondary metadata).

In edit mode, the row exposes:
- Drag handle (left edge) — for reorder.
- Inline-or-expand controls: label override, color swatch, normalization dropdown, q-window numeric inputs.
- "Reset color" button when override is set.

Aligned to plot baseline: row's vertical center matches the y-baseline of the corresponding band (shared y-stacking math).

- [ ] **Step 1:** Tests for review-mode read-only rendering, edit-mode controls, expand/collapse.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): MemberMetaRow + read-only/edit modes`.

### Task 7.2: `BandResizeDivider`

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/BandResizeDivider.tsx`
- Create: `packages/HimalayaUI/frontend/test/BandResizeDivider.test.tsx`

Thin grabbable divider between adjacent metadata rows in edit mode. Drag adjusts `band_height` of the two adjacent members in normalized space (push semantics, total Σ stays constant). Clamp at minimum (`0.1`). Snap to neighbor when within dead zone.

- [ ] **Step 1:** Tests for drag math (Δ resize → both neighbors updated correctly), clamp behavior, snap behavior.
- [ ] **Step 2:** Implement (likely with `@dnd-kit` or a small custom drag handler).
- [ ] **Step 3:** "Reset heights" button in edit-mode header — sets all `band_height` to `1.0`.
- [ ] **Step 4:** Commit as `feat(compare): band resize divider + reset heights`.

### Task 7.3: Reorder via metadata-row drag handle

Reorder uses the existing `@dnd-kit` patterns from Inspect's thumbnail drag. On drop, draft's member ordering updates; `display_order` re-sequenced contiguously. Plot + gutter re-render aligned via shared layout math.

- [ ] **Step 1:** Tests for reorder semantics.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): metadata-row reorder`.

---

## Phase 8: Edit-mode peak click + label rendering

### Task 8.1: Peak click cycle

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MemberTraceLayer.tsx`
- Modify: `packages/HimalayaUI/frontend/src/lib/comparison/draft.ts` (cyclePeakDisplay action)

Triangle markers clickable in edit mode only. Click cycles `shown → labeled → hidden → shown` for the clicked peak. `alt+click` → directly to `hidden`. State persisted in `comparison_members.peak_display` JSON.

- [ ] **Step 1:** Tests for cycle semantics, alt+click direct-to-hidden, state persistence into draft.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): peak click cycle`.

### Task 8.2: Label placement with leader-line dodge

Labels render vertically above triangle markers. 1D dodge layout pass collects all labeled peaks per trace (in q order), spreads them horizontally to avoid collisions, emits leader lines to triangles.

- [ ] **Step 1:** Tests for layout: no overlap, leader lines connect text to original triangle position.
- [ ] **Step 2:** Implement (possibly via Observable Plot's `dodge` mark or a custom layout function feeding raw `text` and `link` marks).
- [ ] **Step 3:** Commit as `feat(compare): peak label placement with leader-line dodge`.

### Task 8.3: Hover tooltip

Triangles get a hover tooltip showing q-value (and optionally peak id when a developer-only `?showPeakIds` URL flag is set, to keep the production UI clean).

- [ ] **Step 1:** Tests.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): peak hover tooltip`.

---

## Phase 9: Review-mode annotation toggles + coloring

### Task 9.1: Coloring library

**Files:**
- Create: `packages/HimalayaUI/frontend/src/lib/comparison/coloring.ts`
- Create: `packages/HimalayaUI/frontend/test/coloring.test.ts`

API:
- `colorFor(member, mode, palette)` — returns the resolved color for a member given the current grouping mode (`bySample` | `byPhase` | `distinct`) and a categorical palette.
- Resolution order: `member.color_override` → grouping-mode default → fallback gray for orphans.
- Palette: ~10–12 distinguishable, accessibility-checked hues. Pick from existing project palette if there is one; otherwise vendor a small constant.

- [ ] **Step 1:** Tests for each mode + override fallback + orphan fallback.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): trace coloring with grouping mode`.

### Task 9.2: Grouping mode toggle (per-tab Zustand)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/state.ts`
- Add UI control to `ComparePage.tsx` header.

Action: `setGroupingMode("bySample" | "byPhase" | "distinct")`. State per-tab (not persisted on the comparison; not persisted in `sessionStorage` since it's a viewing preference).

- [ ] **Steps 1–3.** Test, implement, commit.

### Task 9.3: Annotation toggles

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/AnnotationToggles.tsx`
- Modify: `packages/HimalayaUI/frontend/src/state.ts`

Two checkboxes in review-mode header: peak ticks / per-trace labels. Per-tab Zustand. Toggling re-renders plot marks. Hidden in edit mode.

No predicted-phase-ratio toggle — per spec, we don't render predicted-q ticks at all in v1. The figure is the result of curation, not helpful suggestions.

- [ ] **Step 1:** Tests.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): annotation toggles`.

### Task 9.4: Per-member color override picker

Color swatch grid in `MemberMetaRow`'s expanded state. Click sets `color_override`. "Reset color" clears it.

- [ ] **Steps 1–3.**

### Task 9.5: Hover-driven phase coloring of peaks

**State location: Zustand single-setter** (`setHighlightedCompareMemberId(id: number | undefined)`) — pass `undefined` to clear, matching the existing `hoveredIndexId` pattern on the Index page. Not component-local — the highlight must coordinate across `MemberMetaRow` (hover source) and `MemberTraceLayer` (render target).

**Render mechanism:** `MemberTraceLayer` re-renders its peak marks when `highlightedMemberId` changes — switching from black fill to phase-color fill on existing marks. No full plot recompute.

Hovering a metadata row sets `highlightedCompareMemberId` in Zustand to the member's id. `MemberTraceLayer` reads this and colors that member's `snapshot.confirmed_index.peak_ids` peaks in the phase color; non-index peaks stay black. Click-to-pin makes the highlight sticky; click again unpins. Cleared on page navigation, entering edit mode, or member removal.

**Keyboard equivalent:** rows are focusable (`tabindex="0"`); focus triggers transient highlight same as hover. Enter to pin, Esc to clear.

Members with `confirmed_index === null` have no hover affordance — nothing to highlight. Hover association always reads from `snapshot.confirmed_index.peak_ids` (snapshot data, not live).

- [ ] **Step 1:** Tests for hover state, click-to-pin, keyboard (Tab + Enter + Esc), no-confirmed-index inert path.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): hover-driven phase coloring of peaks`.

### Task 9.6: "Needs Review" badge

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePage.tsx`
- Modify: `packages/HimalayaUI/frontend/src/api.ts` (response includes `is_stale` per member)

Backend `GET /api/comparisons/:id` response includes `is_stale: bool` per member, computed by joining each member's `snapshot.analysis_inputs_hash` against the current exposure's `analysis_inputs_hash`. The comparison is stale if any member is. Per-member staleness shows a warning icon on each stale metadata row (`data-stale="true"`) so the user knows which of N members drifted.

Review-mode header shows a "Needs Review" badge (amber/warning tone, same visual language as `StaleIndicesBanner`) with a tooltip when stale. Clicking the badge (only available to the author) navigates to edit mode. Non-authors see the badge as informational only.

- [ ] **Step 1:** Backend test: create comparison, change underlying exposure, fetch comparison, assert `is_stale: true` on affected member. Also test **low-R² snapshot path**: confirm an index with R²=0.95 on the source exposure, submit comparison, fetch, assert `confirmed_index === null` in snapshot (below-gate confirmations are not snapshotted).
- [ ] **Step 2:** Frontend test: badge visible when any member stale, hidden otherwise; clickability author-gated; per-member `data-stale` attribute correct.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Commit as `feat(compare): Needs Review badge for stale snapshots`.

---

## Phase 10: Chat thread + `@comparison:N` mentions

### Task 10.1: ChatCard generalization for comparison context

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ChatCard.tsx`
- Modify: `packages/HimalayaUI/frontend/src/api.ts` + `queries.ts` for `comparison_messages`.

`ChatCard` today is sample-scoped. Add a context prop: `entityType: 'sample' | 'comparison'` + `entityId: number`. The hook layer chooses which API to call.

- [ ] **Step 1:** Tests verifying both contexts route to correct endpoints.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Mount on review-mode `ComparePage`. Hidden in edit mode.
- [ ] **Step 4:** Commit as `refactor(chat): generic entity context + compare integration`.

### Task 10.2: `@comparison:N` mention support

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MentionPicker.tsx` — Comparisons group + scope rules per [chat-mention-design.md](2026-04-29-chat-mention.md).
- Modify: `packages/HimalayaUI/frontend/src/components/MentionChip.tsx` — render comparison chip with title + member count tooltip + drift indicator when `@comparison:42@<oldhash>` differs from current.
- Modify: `packages/HimalayaUI/frontend/src/hooks/useMentionResolution.ts` — comparison entity type.
- Modify: `packages/HimalayaUI/frontend/src/api.ts` — `getComparison(id)` already exists from Phase 4.
- Modify: `packages/HimalayaUI/frontend/src/lib/renderMentions.tsx` — eager hash insertion on `@comparison:N`.

Eager `@comparison:42@<hash8>` insertion in compose path. Rendering shows live state with a "(changed)" indicator if hash drifts.

- [ ] **Steps 1–4.** Test, implement (mirror existing mention types), commit as `feat(compare): @comparison:N mentions`.

### Task 10.3: Backend mention resolution route

`GET /api/comparisons/:id` already returns the data. The mention resolution layer just needs the `comparison` type added to `useMentionResolution`'s lookup table.

- [ ] **Step 1:** Done as part of Task 10.2.

---

## Phase 11: Authorship + forking

### Task 11.1: Author check route gates

**Files:**
- Modify: `packages/HimalayaUI/src/routes_comparisons.jl`
- Modify: `packages/HimalayaUI/test/test_comparisons.jl`

`is_author(db, comparison_id, user_id)` was added in Task 2.1. This task wires it into the submit and delete routes (these already returned 403 for non-authors per Task 2.2's tests; this is bookkeeping confirming the gate is consistent across all author-gated affordances).

`created_by IS NULL` (orphaned author) → no user matches `is_author`; route returns 403 with "no author" detail.

- [ ] **Step 1:** Verify the existing tests (added in Task 2.2) still cover all author-gated paths.
- [ ] **Step 2:** Add explicit tests for the orphaned-author case (delete the user row, attempt to submit as anyone, expect 403).
- [ ] **Step 3:** Commit as `test(compare): orphaned-author 403 case`.

### Task 11.2: Frontend Edit-vs-Fork affordance

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/ComparePage.tsx`
- Modify: `packages/HimalayaUI/frontend/src/api.ts` (current user query)

Header shows Edit button when `comparison.created_by === currentUser.id`; otherwise Fork. Fork button creates a new draft pre-populated with the parent's snapshot, navigates to edit mode.

- [ ] **Step 1:** Tests for visibility logic + fork-creates-correct-draft.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): author-only edit + fork action`.

### Task 11.3: Lineage badge + forks popover

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/LineageBadge.tsx`

Review-mode header shows "Forked from <parent title> by <author> [view parent →]" when `forked_from_id` is set; "Forked from a deleted comparison" when `forked_from_id IS NULL` and `forked_at_hash` is set.

Forks popover: "Forks (N) →" link that opens a list of child comparisons via `GET /api/comparisons/:id/forks`.

- [ ] **Step 1:** Tests.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): lineage badge + forks popover`.

---

## Phase 12: Conflict modal

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ConflictModal.tsx`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx`
- Create: `packages/HimalayaUI/frontend/test/ConflictModal.test.tsx`

Mounted at `App.tsx`. Listens for `ConflictError` from `saveComparison` mutator. Modal shows side-by-side: server's current state (from the **409 response body**, never `qc.getQueryData`) vs the local draft. Wrap the modal container with `useFocusTrap(containerRef, isOpen)`. Actions:

- *Discard my changes* — clears draft, navigates to review of server state.
- *Overwrite with mine* — re-submits with `expected_content_hash` set to server's `current_hash`. **Second-409 race:** if another tab submits between the 409 and the re-submit, the re-submit itself returns a second 409 with updated `current_state`. The modal must handle re-opening with the new server state (refresh the diff view, let the user pick again). Rare but must not break.
- *Fork to a new comparison* — calls saveComparison with no `id`, with `forked_from_id` set to the original. New comparison opens.

Each modal action mints a fresh `client_op_id` (via `useQueueMutation`'s per-`mutate()` mint). Note: per `with_idempotency`, 409 responses (status ≥ 400) are not cached, so a retry with the same `client_op_id` re-evaluates the conflict check — no stale cached 409 to break out of.

- [ ] **Step 1:** Tests including: basic 409 → modal opens; overwrite succeeds; overwrite hits second 409 → modal re-opens with updated state; fork creates new comparison.
- [ ] **Step 2:** Implement.
- [ ] **Step 3:** Commit as `feat(compare): conflict modal with second-409 handling`.

---

## Phase 13: Polish + e2e

### Task 13.1: Skeleton loading

Per CLAUDE.md boneyard rules: every load-gated card on Compare gets a `<Skeleton>` wrapper gated on `query.isLoading` (NOT `isPending` — disabled queries and background refetches stay skeleton-free). List page, plot card, member panel, chat. Boneyard captures committed.

- [ ] Implement.

### Task 13.2: Sidebar pin (deferred from Task 4.2)

Per-user pinned comparisons surface at the top of the sidebar. Backend: small `comparison_pins(user_id, comparison_id)` table or per-user JSON in a settings table. Discuss in PR which fits.

- [ ] Backend table + routes + frontend pin/unpin.

### Task 13.3: e2e smoke

**Files:**
- Create: `packages/HimalayaUI/frontend/e2e/compare.spec.ts`

Smoke test: create comparison with 2 members → submit → reload → fork → edit → submit. Mocked `/api` per existing Playwright pattern (`page.route`).

- [ ] Implement.

### Task 13.4: Final polish

- Empty state for sidebar.
- Keyboard shortcuts (Esc to cancel edit, Cmd+Enter to submit).
- Tooltip cleanup.
- Visual review against [himalayaui-design.md](../../himalayaui-design.md) §1 (plot-dominant, recede chrome).
- **Verify `data-testid` coverage**: every affordance in the E2E selector table (Phase 0 section) should have its `data-testid` added in the same phase that builds the component — not bolted on here.

---

## Validation

Before opening the implementation PR:

- [ ] `Pkg.test("HimalayaUI")` green
- [ ] `Pkg.test()` (core Himalaya) green
- [ ] `npm test` (Vitest) green
- [ ] `npm run build` green (TS strict)
- [ ] `npm run e2e` green
- [ ] Manual exploratory test of: create → submit → fork → edit → conflict (open same comparison in two tabs as same user, edit both, submit both) → conflict modal flow → discard / overwrite / fork-to-new
- [ ] Visual sanity check across light + dark mode

---

## Open follow-ups (post-v1)

Captured here so they aren't lost:

- **Generalize `sample_messages` + `comparison_messages` into one `messages` table** with `entity_type` discriminator. Out of scope for this plan; the parallel-tables approach was chosen to keep this work isolated.
- **Server-side draft persistence** — survives across devices and tabs. Probably a `comparison_drafts(user_id, comparison_id?, snapshot, updated_at)` table.
- **PR-style merge from forks back to parent.** Requires a diff UI and a curation workflow.
- **Aspect-ratio control surfacing** — only if user feedback demands it.
- **`working_band_fraction` schema field** — only if per-comparison tuning is asked for.
- **Replicate-shading** — when the broader codebase grows a "what counts as a replicate" concept.
- **Phase chips in picker** — if pick-by-trace-property turns out to be missed in practice.
- **Color-blind safe palette repo-wide** — possibly forced by this work; coordinate with [himalayaui-design.md](../../himalayaui-design.md) §3 open question.
