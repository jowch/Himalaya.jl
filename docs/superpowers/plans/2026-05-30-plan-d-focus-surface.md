# Plan D — Focus Surface Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the redesigned Focus surface — an assignment cart (0..N phases + the 3-state indexed/form_factor/null model), reflection combs with an indexing-space residual toggle, phase-coloured detector rings, the q-link triple, a plot-only hypothetical-candidate preview, and a snap-to-peaks custom-index modal — moving the frontend onto the native assignment endpoints and retiring the legacy `/groups` machinery + dual-write.

**Architecture:** A new `useAssignment` data layer (TanStack Query against `GET /api/exposures/{id}/assignment`) plus four queue mutators (`assignment_add`, `assignment_remove`, `assignment_set_state`, and the custom-index commit) replace the `GroupEntry.active` single-active model. A `deriveActiveIndices(assignment, indices)` helper replaces all `groups.find(g=>g.active)` reads. The combs panel and detector rings render on `PlotSurface.overlay` (Plan C) using `peakMark` glyphs and `phaseColor`. The hypothetical preview is pure ephemeral Zustand (`previewPhaseId`), never a mutator. The custom-index modal composes `ui/` primitives. The Bonnet badge consumes `IndexEntry.bonnet` (Plan B). When the cart is live, the legacy `index_groups` tables, `/groups` routes, dual-write, and `ensure_custom_group!` are deleted (the cleanup half of Plan A's transitional shim).

**Tech Stack:** TypeScript, React, TanStack Query, the mutation queue (`useQueueMutation`/`applyRemoteToCache`/`mutatorRegistry`), Zustand, Observable Plot via `PlotSurface`, the `ui/` design system, Vitest + the six-layer contract-test rule, the mockup `docs/redesign-mockups/2026-05-29-focus-plot.html` as the pixel/interaction spec.

---

## Dependencies & scope

- **Depends on Plan A** (assignment endpoints + 3 event kinds), **Plan B** (`IndexEntry.bonnet`), **Plan C** (`PlotSurface`, `peakMark`).
- **Backend item B4** (lattice-driven speculative build path for the custom-index commit) is delivered as Task D-9 here (small backend addition) since it's tightly coupled to the modal.
- **Cleanup:** retiring legacy `/groups` + dual-write is the final task; until then both run (Plan A guarantee).
- The mockup is the visual/interaction contract for combs, residual view, preview ghosting, custom-index snap, and the cart block layout (× remove in the header, `+ custom index…` footer button, contextual Bonnet note appearing only when substantive).

## The 3-state model on the wire (from Plan A)

`GET /api/exposures/{id}/assignment` → `{ exposure_id, state: "indexed"|"form_factor"|"null", members: number[] }`. **State is explicit, never inferred from `members.length`.** `indexed` with 0 members = "call in progress"; `form_factor`/`null` always have 0 members. The cart UI renders: phase blocks (indexed), or a single form-factor/null declaration chip.

## Must-preserve queue invariants (from the survey)

- New mutators mint `client_op_id` **inside** `mutationFn` (inherited via `useQueueMutation`), negative optimistic ids, rollback symmetry, per-iteration replay try/catch.
- Assignment SSE events need their own `applyRemoteToCache` arms **and** a **distinct post_state shape** (`{assignment:{state,members}}`) so the `CurationPostState` cast (applyRemoteToCache.ts:18-41, reads `ps.indices`) never clobbers cache. Register all 4 kinds in **both** `resolveMutator` and `resolveMutatorForEvent`.
- `affectsExposurePeaks: () => false` for assignment mutators (they don't touch peaks); the preview never sets a deferred or emits an event (`affectsExposurePeaks` stays false; hover must not gate on `useExposureHasPendingPeakOps`).

## File structure

| File | Change | Responsibility |
|---|---|---|
| `src/api.ts` | Modify | `Assignment` type; `getAssignment`/`setAssignmentState`/`addAssignmentPhase`/`removeAssignmentPhase` fetchers |
| `src/queries.ts` | Modify | `useAssignment`; `useAddAssignmentPhase`/`useRemoveAssignmentPhase`/`useSetAssignmentState` hooks; `queryKeys.assignment` |
| `src/lib/queue/mutators/assignment.ts` | **Create** | 3 mutators (`assignment_add`/`assignment_remove`/`assignment_set_state`) |
| `src/lib/queue/mutators/customIndex.ts` | **Create** | custom-index commit mutator |
| `src/lib/queue/mutatorRegistry.ts` | Modify | register the 4 kinds in both dispatch tables |
| `src/lib/queue/applyRemoteToCache.ts` | Modify | `assignment_*` case arms + assignment post_state arm |
| `src/lib/queue/types.ts` | Modify | extend `OpKind` with the 4 kinds |
| `src/lib/assignment.ts` | **Create** | `deriveActiveIndices(assignment, indices)`, `assignmentState(assignment)` helpers |
| `src/state.ts` | Modify | `previewPhaseId` (ephemeral); custom-index modal gate; remove single-active assumptions |
| `src/components/AssignmentCart.tsx` | **Create** | the cart (phase blocks + 3-state + Bonnet note + custom-index footer) |
| `src/components/PhasePanel.tsx` | Rewrite | candidate list + cart + preview wiring (replaces `activeGroup` reads) |
| `src/components/CombPanel.tsx` | **Create** | reflection combs + residual toggle (replaces `FocusReflectionsTable`) |
| `src/components/DetectorRingOverlay.tsx` | Modify | phase-colour rings, concentric under coexistence, hollow ghost ring, comb-tooth q-link node |
| `src/components/CustomIndexModal.tsx` | **Create** | symmetry + lattice + live preview + snap-to-peaks |
| `src/components/PlotCard.tsx` | Modify | `losingPeakIds` → derive from cart; preview ghost row |
| `packages/HimalayaUI/src/routes_analysis.jl` | Modify | B4: accept a client-fitted basis on the speculative POST; **delete** `/groups` routes + dual-write + `ensure_custom_group!` (cleanup) |
| `packages/HimalayaUI/src/db.jl` | Modify | (cleanup) drop `index_groups`/`index_group_members` after migration confidence |
| tests | Create/modify | six-layer per mutator + component tests |

---

## Task sequence

> **Sequencing (survey):** cart state model + SSE post_state shape FIRST; then the 4 OpKinds + applyRemoteToCache arms + registry BEFORE cart UI; cart semantics BEFORE detector rings (render per-phase from cart) and BEFORE preview.

### Task D-1: `Assignment` type + `useAssignment` data layer

**Files:** `src/api.ts`, `src/queries.ts`.

```typescript
// api.ts
export type AssignmentState = "indexed" | "form_factor" | "null";
export interface Assignment {
  exposure_id: number;
  state: AssignmentState;
  members: number[]; // index ids
}
export const getAssignment = (id: number) =>
  request<Assignment>("GET", `/api/exposures/${id}/assignment`);
export const setAssignmentState = (id: number, state: AssignmentState, opts?: AuthOpts) =>
  request<AssignmentMutationResponse>("POST", `/api/exposures/${id}/assignment/state`, { state }, opts);
export const addAssignmentPhase = (id: number, index_id: number, opts?: AuthOpts) =>
  request<AssignmentMutationResponse>("POST", `/api/exposures/${id}/assignment/members`, { index_id }, opts);
export const removeAssignmentPhase = (id: number, index_id: number, opts?: AuthOpts) =>
  request<AssignmentMutationResponse>("DELETE", `/api/exposures/${id}/assignment/members/${index_id}`, undefined, opts);
```

> **Backend addendum (small):** Plan A delivered `GET /assignment` + `POST /assignment/state`. Add the native `POST/DELETE /api/exposures/{id}/assignment/members` routes (emit `assignment_add`/`assignment_remove`, mirror the Plan A `/assignment/state` route shape with `with_idempotency` + `apply_event!(InTransaction())` + `_assignment_body` response). This replaces the dual-write target.

- [ ] TDD: backend route tests (mirror Plan A's `assignment/state` test) → implement native member routes. Frontend: `useAssignment` query hook + `queryKeys.assignment(id)`.
- [ ] Commit: `feat(focus): assignment data layer + native member routes`.

### Task D-2: `deriveActiveIndices` + replace all `find(g=>g.active)` reads

**Files:** Create `src/lib/assignment.ts`; modify `PhasePanel.tsx:13`, `PlotCard.tsx:231`, `FocusReflectionsTable.tsx:45` (the latter superseded by CombPanel later, but unblock now).

```typescript
// lib/assignment.ts
export function deriveActiveIndices(a: Assignment | undefined, indices: IndexEntry[]): IndexEntry[] {
  if (!a || a.state !== "indexed") return [];
  const ids = new Set(a.members);
  return indices.filter((ix) => ids.has(ix.id));
}
```

- [ ] TDD: unit test `deriveActiveIndices` (empty for form_factor/null; filters by members). Replace the three reads to consume `useAssignment` + `deriveActiveIndices`. Keep behavior identical (still single coherent active set, now assignment-sourced).
- [ ] Commit: `refactor(focus): deriveActiveIndices replaces groups.find(g=>g.active)`.

### Task D-3: The 4 mutators + applyRemoteToCache arms + registry (six-layer)

**Files:** `src/lib/queue/mutators/assignment.ts`, `mutatorRegistry.ts`, `applyRemoteToCache.ts`, `types.ts`.

Model each mutator on `addIndexToGroupMutator` (captured exemplar), keyed on `queryKeys.assignment(exposureId)`:

```typescript
export const addAssignmentPhaseMutator: Mutator<AddPhaseInput, AssignmentScope, AssignmentMutationResponse> = {
  kind: "assignment_add",
  onMutate: (p, qc) => {
    const key = queryKeys.assignment(p.exposureId);
    const prev = qc.getQueryData<Assignment>(key);
    if (prev) qc.setQueryData<Assignment>(key, {
      ...prev, state: "indexed",
      members: prev.members.includes(p.indexId) ? prev.members : [...prev.members, p.indexId],
    });
    return { restore: () => { if (prev !== undefined) qc.setQueryData(key, prev); } };
  },
  request: (p) => api.addAssignmentPhase(p.exposureId, p.indexId, buildAuth(p)),
  onSuccess: (p, response, qc) => {
    const { payload: row } = stripQueueMetadata(response);
    qc.setQueryData<Assignment>(queryKeys.assignment(p.exposureId), row as Assignment);
  },
  synthesizeFromSse: (remote, base) => ({ ...base, /* {state, members} from post_state */ } as AssignmentMutationResponse),
  affectsExposurePeaks: () => false,
};
```

`assignment_remove` mirrors (filter member; `treats404AsSuccess: true`). `assignment_set_state` sets `state` + clears `members` when non-indexed.

`applyRemoteToCache` arms:
```typescript
case "assignment_add":
case "assignment_remove":
case "assignment_set_state": {
  const ps = remote.post_state as { assignment?: Assignment } | undefined;
  if (ps?.assignment) qc.setQueryData(queryKeys.assignment(id), ps.assignment);
  else qc.invalidateQueries({ queryKey: queryKeys.assignment(id) });
  break;
}
```

**Distinct post_state:** the backend assignment routes must emit `post_state = {assignment: {state, members}}` (NOT the curation `{indices}` shape) so `applyPostStateOnly`'s `Array.isArray(ps.indices)` guard skips it. Add this to the Plan A/D-1 routes.

- [ ] TDD six-layer per kind (route emit, SSE frame, applyRemoteToCache merge, cache-row shape, onMutate↔rollback inverse, onSuccess convergence) — mirror `mutatorOnSseWins.test.ts` + `hooks.test.tsx`.
- [ ] Register all 4 in `resolveMutator` + `resolveMutatorForEvent`; extend `OpKind`.
- [ ] Commit: `feat(queue): assignment_add/remove/set_state mutators + SSE arms`.

### Task D-4: `AssignmentCart` component

**Files:** Create `src/components/AssignmentCart.tsx`; compose `PhaseChip`/`ScoreBar`/`Button`/`IconButton` (× remove in block header), `SegmentedControl` for the 3-state (Indexed / Form factor / No scattering), the contextual Bonnet note (only when substantive: suggestion when a cubic is in and peaks unexplained; consistency `ratio ≈ 1.279` when two cubics coexist — consume `IndexEntry.bonnet`), and a `+ custom index…` footer `Button` (must read as a button on hover — fixes the mockup-noted hover gap). Visual contract: `docs/redesign-mockups/2026-05-29-focus-plot.html`.

- [ ] TDD: component tests — renders N phase blocks for `indexed`; renders the form-factor/null declaration for those states; × removes via `useRemoveAssignmentPhase`; SegmentedControl switches state via `useSetAssignmentState`; Bonnet note appears only when substantive; footer button has hover affordance.
- [ ] Commit: `feat(focus): AssignmentCart (3-state, Bonnet note, custom-index footer)`.

### Task D-5: `CombPanel` (combs + indexing-space residual toggle)

**Files:** Create `src/components/CombPanel.tsx` (replaces `FocusReflectionsTable`). Render per-phase teeth on a shared log-q ruler from `IndexEntry.predicted_q` (all orders), full-height gridlines tying teeth to observed peaks, hollow-caret (`peakMark` predictedAbsent) for predicted-but-absent, a leftover row of unexplained observed peaks; a `SegmentedControl` toggles to the residual lollipop view (`index_peaks.residual` from `IndexEntry.peaks[].residual`, peaks bending off the Δq/q=0 line). Drive via `PlotSurface.overlay` or a dedicated SVG ruler. q-link: comb teeth respond to `hoveredQ` (the third node) — grow + terracotta ring, not recolour.

- [ ] TDD: combs render a tooth per predicted order; absent orders render hollow carets; leftover peaks appear in the leftover row; residual toggle switches views; hovering `hoveredQ` lights the matching tooth.
- [ ] Commit: `feat(focus): CombPanel with indexing-space residual toggle`.

### Task D-6: Detector rings in phase colour + comb-tooth q-link node

**Files:** `src/components/DetectorRingOverlay.tsx` (ring colour logic at 97-114). Colour rings by assigned phase (`phaseColor`); render concentric ring sets under coexistence; hollow ghost ring for predicted-but-absent; keep the `hoveredQ` sink (74-75) and add the comb tooth as the third q-link node so peak↔ring↔tooth light together. Threads phase colour through SVG attributes (no `style` literal → design-guard safe, matching the existing pattern).

- [ ] TDD: rings colour by phase; two concentric sets under coexistence; ghost ring for absent order; q-link lights ring+tooth together for a hovered q.
- [ ] Commit: `feat(focus): phase-colour detector rings + comb-tooth q-link node`.

### Task D-7: Hypothetical-candidate preview (plot-only, ephemeral)

**Files:** `src/state.ts` (`previewPhaseId`/`previewIndexId` ephemeral, omitted from `partialize`), `PhasePanel.tsx` (candidate hover sets/clears preview), `PlotCard.tsx` (rewrite `losingPeakIds` `alreadyActive` check from "in activeGroupIndices" to "in any cart member"; render a dashed ghost comb row + trace highlight from the hovered candidate's `predicted_q`). **Never** a mutator, never an SSE event, `affectsExposurePeaks` stays false; `onHoverLeave`/blur MUST clear (missing clear = stale ghost masking the real cart). Reuse the `losingPeakIds` machinery.

- [ ] TDD: hovering a candidate sets `previewPhaseId` and renders a ghost comb (plot-only); leaving clears it; preview never enqueues a mutation; a weak candidate's ghost shows an all-caret (explains few peaks) row.
- [ ] Commit: `feat(focus): plot-only hypothetical-candidate preview`.

### Task D-8: PhasePanel rewrite (cart + candidates + preview wired)

**Files:** Rewrite `src/components/PhasePanel.tsx` to compose `AssignmentCart` + the candidate list (with `inCall` derived from the assignment, Bonnet badge from `IndexEntry.bonnet`, hover → preview) + the `CustomIndexModal` gate. Remove all `activeGroup`/`useAddIndexToGroup`/`useRemoveIndexFromGroup` usage in favor of the assignment hooks.

- [ ] TDD: candidate `inCall` reflects assignment membership; clicking toggles via assignment mutators; Bonnet badge shows on consistent candidates; custom-index footer opens the modal.
- [ ] Commit: `refactor(focus): PhasePanel on the assignment cart`.

### Task D-9: Custom-index modal + B4 lattice-build backend

**Files:** Create `src/components/CustomIndexModal.tsx` (compose `ModalShell` + `SegmentedControl` symmetry picker [Pn3m/Im3m/Ia3d/Lamellar/Hexagonal] + lattice numeric input/slider + a live-preview comb computed client-side via real physics [`2π√N/a` cubic, `2πn/d` lamellar, `4π√M/(√3·a)` hex] + a running "lands on N of M observed peaks" fit + snap-to-peaks [magnetic drag to values where a predicted order lands on an observed peak; click an observed peak to snap first reflection]). Commit via a new `customIndexMutator`. **Backend B4:** extend the speculative POST (`routes_analysis.jl:297+`) to accept a client-fitted `{phase, basis}` directly (not only anchor_peak+ratio), persist as a candidate index, then add it to the assignment.

- [ ] TDD: backend — speculative POST accepts `{phase, basis}`, persists an index, returns it. Frontend — modal computes reflections from physics; snap pulls lattice to a peak (verify 200→198, 196→197, 215 stays free per the mockup); commit adds a candidate + assignment member.
- [ ] Commit: `feat(focus): custom-index modal + lattice-driven speculative build (B4)`.

### Task D-10: Retire legacy `/groups` + dual-write (cleanup)

**Files:** `routes_analysis.jl` (delete `GET/POST/DELETE /groups*`, `ensure_custom_group!`, the dual-write in the member routes), `db.jl` (drop `index_groups`/`index_group_members` once the assignment is the sole source — guard behind a confirmed migration), remove `addIndexToGroupMutator`/`removeIndexFromGroupMutator` + their registry entries + `applyRemoteToCache` index_confirmed/unconfirmed arms, delete `FocusReflectionsTable.tsx`.

- [ ] TDD: full backend suite green without the legacy routes; full frontend suite + e2e green without the group mutators; verify no consumer references `GroupEntry`/`useGroups`.
- [ ] Commit: `refactor: retire legacy index_groups machinery + dual-write (assignment is sole source)`.

### Task D-11: Full green + visual acceptance + a11y floor

- [ ] Backend + frontend suites green (capture once).
- [ ] `npm run build` (tsc + lint:design) + `npm run e2e`.
- [ ] Visual acceptance vs the Focus mockup (combs, residual, rings, preview, custom-index snap, cart × button alignment, footer hover).
- [ ] a11y floor from the mockup audit P1: cart blocks / candidate rows are `<button>`s (keyboard + role); q-link has a focus-driven path (focus a peak → same highlight), not hover-only.
- [ ] Commit: `test(focus): visual + a11y acceptance for the Focus surface`.

---

## Self-Review

**1. Spec coverage** (F1–F9, design §3/§8.1/§8.2):
- F1 cart 0..N + 3-state + form-factor → D-1..D-4, D-8. ✓
- F2 ranked candidates + a/d label → D-8. ✓
- F3 Bonnet badge/note (consumes Plan B) → D-4, D-8. ✓
- F4/F5 combs + residual → D-5. ✓
- F6/F9 phase rings + q-link triple → D-6. ✓
- F7 plot-only preview → D-7. ✓
- F8 custom-index snap + B4 → D-9. ✓
- Legacy retirement → D-10. ✓

**2. Placeholder scan:** the load-bearing new contracts (Assignment type, mutator shape, applyRemoteToCache arms, deriveActiveIndices, distinct post_state) have full code; the breadth of components is task-level with explicit acceptance + the mockup as the pixel spec + the captured `ui/` primitive signatures to compose. No "add error handling"/"similar to".

**3. Type/name consistency:** `Assignment`/`AssignmentState`/`assignment_add|remove|set_state`/`queryKeys.assignment`/`deriveActiveIndices`/`previewPhaseId`/`AssignmentMutationResponse` consistent across api, queries, mutators, registry, cache, components. `post_state = {assignment:{state,members}}` consistent backend↔frontend.

**Risks:** (a) issue #37 Bug 1 — a fresh custom index id not yet in cache: on custom-index commit, invalidate `queryKeys.assignment` rather than splice. (b) Don't ship D-10 until migration confidence is high (Plan A's `rebuild_views_from_log!` round-trip + a backfill check on the real DB).

---

## Review findings — required fixes (queue-reviewer, himalaya-reviewer, saxs-physics-reviewer, 2026-05-30)

1. **[HIGH] Extend the `SseEvent.post_state` type union (D-3 types.ts).** `types.ts:124` types `post_state?: CurationPostState | Comparison | Series` — the plan's `remote.post_state as {assignment?: Assignment}` cast won't compile. Add `export interface AssignmentPostState { assignment: Assignment }` and widen the union to include it (the File-structure table says "extend `OpKind`" only — also extend post_state). The backend (Plan A fix #1) must emit `post_state = {assignment:{state,members}}` or D-3's cache arm never sees `ps.assignment`.
2. **[HIGH] Assignment `onSuccess` writes ONLY `queryKeys.assignment`, never `queryKeys.exposure`.** `synthesizeResponseFromSse` lifts `analysis_inputs_hash` from post_state (`replayCoordinator.ts:162`) → `undefined` for an assignment frame. Unlike `peakAdd` (which writes the exposure hash), the assignment mutators must not, or they wipe the exposure's hash with `undefined`. Pin in the D-3 mutator spec.
3. **[HIGH] D-10 must no-op-return the BACKEND dispatcher branches, not delete them.** Dropping `index_group_members` while historical `index_confirmed`/`index_unconfirmed` events remain in `user_actions` makes `rebuild_views_from_log!` throw (`no such table`) on replay. Convert the `events.jl:344-358` branches to `kind == "…" && return nothing` guards (the retired-kind pattern at events.jl:375-376) instead of removing them. Also **bump `schema_version`** in the D-10 commit so a tab with a queued `index_confirmed` op that reloads post-cutover drops-with-toast rather than silently no-op-resolving.
4. **[MED] D-9 custom-index `basis` convention + Ia3d round-trip test.** The modal must submit `basis` = q₁ slope = `2π/a × first(phaseratios(P))`, **not** `a` and **not** `2π/a` — else the speculative route's `fit()` re-derivation disagrees by a factor of `first(phaseratios(P))`. Add a contract test asserting `predicted_q_for_phase(phase, basis)` reproduces the modal's client-side comb for the same `a`, using **Ia3d** specifically (`√6` first reflection maximizes the convention-mismatch signal; Pn3m's `√2` can mask it). Risk (a) — invalidate `queryKeys.assignment` on a fresh custom-index id (issue-#37 pattern), never splice a negative placeholder.
5. **[MED] Pin two tests + split the six-layer post_state checkbox.** (a) `applyPostStateOnly` is a **no-op on an assignment frame** (the single highest-value test — currently implicit). (b) Two pending `assignment_add` ops on the same exposure: reverse-rollback + insertion-replay converge to both members present. Split D-3's six-layer "route emit + SSE frame" into two checkboxes asserting the distinct `{assignment}` post_state at **both** the Julia route-emit test and the `sseEventPayload.contract.test.ts` layer.
6. **[LOW] D-9 `QNumInput` focus-gate:** the custom-index numeric input is externally updated by snap-drag, so follow the `QNumInput` focus-gate pattern (it writes to the same draft, so not a defect — just a one-line note).
7. **Confirmed sound:** client_op_id mint site (inherited from `useQueueMutation`), no optimistic-id concern (members are existing positive index ids; rollback snapshots the whole `Assignment` so it's symmetric across members+state), distinct-post_state-shape dodge + invalidate fallback, dual-table registration of all 4 kinds, replay participation (disjoint `assignment`/`peaks` caches, per-iteration try/catch), preview isolation (ephemeral Zustand, `affectsExposurePeaks:false`, omitted from partialize, blur-clear), D-10 replay reads the op log via `resolveMutator` (not a hardcoded kind set).

## Execution Handoff
1. **Subagent-Driven (recommended)** — fresh subagent per task; the queue tasks (D-3) and cleanup (D-10) get careful two-stage review.
2. **Inline Execution** — batch with checkpoints after D-3, D-8, D-10.
