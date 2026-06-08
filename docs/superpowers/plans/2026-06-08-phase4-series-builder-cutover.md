# Greenfield Series-Builder Cutover (`/series/:id`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the legacy `/series/:id` builder page body with a greenfield `src/print/pages/SeriesBuilderPage.tsx` — a single always-"Compose" surface assembled from the built `BuilderRail` / `SeriesPlate` / `MemberList` composites, on a **lazy-draft** model, reusing the carried series/draft logic and deleting the legacy presentation.

**Architecture (Plan 6b — the last Series cutover):**
- **Lazy draft:** opening `/series/:id` shows the committed series in the Compose chrome, read-only. View controls (offset, scale, grouping, annotation toggles, collapse) are **local `useState`** and never start a draft. The **first recipe edit** (title/desc/add/remove/reorder a sample/ordering variable/order rule) calls `startSeriesDraftFromSeries(series)` (idempotent) then the matching draft mutator.
- **Confirm = Save → Commit chain (the correctness-critical path):** the backend commit takes `members` verbatim with NO re-resolution, while Save (PATCH recipe) is what re-resolves the plate. So "Confirm series" must: `save.mutate(buildSeriesSaveBody(draft))` → on `save.isSuccess`, read the **fresh `members` from the save response** → `commit.mutate(buildSeriesCommitBody(freshMembers))` → on `commit.isSuccess`, `discardSeriesDraft()` + stay on `/series/:id` (now read state). Chained via pending refs (queue mutations expose `.mutate` + `isSuccess`/`data`, not `mutateAsync`). Never commit a stale plate.
- No autosave (mockup shows only "Confirm series", no Save button). No conflict UI (LWW — Plan 6a relaxed the 409). The draft is sessionStorage-mirrored for crash recovery (already built).

**Tech Stack:** React 18 + TypeScript (strict, `exactOptionalPropertyTypes`), TanStack Query (queue-routed mutations), Zustand (seriesDraft + annotation toggles), Vitest + RTL, Playwright.

---

## The wiring table (the heart of the page)

| User action | Carried mutator/hook | Draft? | Save/Commit/Local |
|---|---|---|---|
| Edit title / description | `setSeriesDraftTitle` / `setSeriesDraftDescription` | start-if-none + mutate | local draft only (flushed on Confirm) |
| Add sample | `addSeriesSample(sampleId)` (id from `useCorpusSamples`, corpus-wide) | start + mutate | local draft |
| Remove sample | `removeSeriesSample(rowId)` (by local row id) | mutate | local draft |
| Reorder sample | `reorderSeriesSample(from, to)` (array indices) | mutate | local draft |
| Change ordering variable / order rule | `setSeriesOrderingVariable(v\|null)` / `setSeriesOrderRule(rule)` | start + mutate | local draft |
| Change offset / scale / grouping / collapse | local `useState` | **no** | local-only |
| Toggle peak ticks / labels | `<AnnotationToggles/>` (self-wires Zustand) | **no** | local-only (Zustand) |
| **Confirm series** | `useSaveSeries` → then `useCommitSeriesPlate` (chain) | reads + discards on success | **Save→Commit** then `discardSeriesDraft()` (no navigation — stay, now read) |
| Cancel / discard | `discardSeriesDraft()` | drops draft | local-only (no request, no nav) |

`OrderRule` = `"ascending" | "descending" | "manual"`. `SeriesRecipeRow` = `{ id; sample_id; position; pinned; excluded }` (`id` is a local render key, never wired). `buildSeriesCommitBody(members)` is one-arg/hash-free (post-6a).

---

## Provenance ledger

**NEW (create):** `src/print/pages/SeriesBuilderPage.tsx`, `src/print/pages/builderAdapters.ts`, tests `test/print-pages/builderAdapters.test.ts` + `test/print-pages/SeriesBuilderPage.test.tsx`.

**CARRIED (reuse, don't touch):** `useSeries`/`useSeriesTraces`/`useSaveSeries`/`useCommitSeriesPlate`/`useCorpusSamples` (`queries.ts`); the `seriesDraft` state machine + `startSeriesDraftFromSeries`/`discardSeriesDraft`/`setSeriesDraft*`/`addSeriesSample`/`removeSeriesSample`/`reorderSeriesSample` (`state.ts`); `buildSeriesSaveBody`/`buildSeriesCommitBody` (`lib/series/`); `toWaterfallRows`/`dominantPhase`/`resolveState` (`print/waterfall/waterfallModel.ts`); `GroupingModeToggle`/`AnnotationToggles`/`FigureExportControls`/`SeriesPhaseStrip` (`src/components/` — carried logic-bearing controls, reused as-is) + `buildMultiTraceExportSpec` (export spec). The greenfield `BuilderRail`/`SeriesPlate`/`MemberList`/`MemberRow`/`AutoGroup` composites.

**OLD (delete on cutover — grep-enumerate, don't trust this list):** legacy `src/pages/SeriesBuilderPage.tsx`, `src/components/SeriesRecipeEditor.tsx`, `SeriesBuilderRail.tsx`, `SeriesMemberRow.tsx`, `SeriesReadingPanel.tsx`, `MultiTracePlot.tsx`, `MemberMetaGutter.tsx`, `OffsetDock.tsx` + their tests — IF grep shows no surviving non-legacy consumer. (Keep `GroupingModeToggle`/`AnnotationToggles`/`FigureExportControls`/`SeriesPhaseStrip` — the greenfield page imports them.)

---

## Scope decisions (Rule-3 — settled; don't re-litigate)

| # | Decision | Why |
|---|---|---|
| 1 | **Lazy draft** (start on first recipe edit; view controls local) | The chosen mode model (2026-06-08); reconciles the always-Compose mockup with "view controls are free, recipe edits are drafted". |
| 2 | **Confirm = Save→Commit chain**, no autosave, no separate Save button | Mockup shows only "Confirm series"; the backend commits `members` verbatim, so the recipe must be saved (plate re-resolved) before the commit reads `members`. |
| 3 | **Representation (waterfall/heatmap) + Track-reflections OMITTED** | `HeatmapChart` ⛔ out-of-scope, `TrackingLine` ⏸ deferred; `BuilderRail` already omits both. Only the scale (log/lin) + offset + grouping view controls ship. |
| 4 | **Keep annotation toggles (peak ticks/labels); drop grouping-mode persistence to the draft** | Resolved decision (2026-06-07): annotation toggles are in the mockup + schema; grouping is a local view control, not a recipe edit. |
| 5 | **`AutoGroup` "Adjust" = start editing; "Confirm series" = the commit chain.** In read state (no draft) show the grouping summary + "Adjust"; with a live draft, "Confirm series" runs the chain and "Cancel" discards. | The compose card's two actions map onto the lazy-draft read/edit sub-states without a separate page mode. |
| 6 | **Inline-editable title/description** via the draft setters (dotted-underline affordance per mockup `.editable`). | Mockup; matches the carried `setSeriesDraftTitle/Description`. |

---

## Task 1: Builder view-model adapters (`builderAdapters.ts`)

Pure mappers from carried API data → the composites' props, + the Confirm-chain gate helpers. Pure + tested.

**Files:** Create `src/print/pages/builderAdapters.ts`; Test `test/print-pages/builderAdapters.test.ts`.

- [ ] **Step 1: Write the failing test** — cover:
  - `membersToMemberData(members)` → `MemberDatum[]` (`{key, phases, variableValue, dataLine}`) for `MemberList`. Derive `phases` from each member's snapshot via the carried `dominantPhase`; `key` = stable member id string; `variableValue` = `label_override ?? ordering label`; `dataLine` = the mono lattice·q₁ line (match `SeriesMemberRow`'s `MemberDatum` shape — read `MemberList.tsx`/`SeriesMemberRow.tsx` first and mirror it exactly).
  - `recipeRowView(row, sampleNameById)` → `{ name, position }` for an editable recipe row.
  - `addableSamples(corpusSamples, draftRecipe)` → corpus samples whose id is NOT already in `draft.recipe` (the add-sample dropdown options).
  - `groupingSummary(series)` → the AutoGroup body copy (e.g. "Himalaya read N samples as one series, ordered by {ordering_variable}").
  - `legendPhasesOf(rows)` → distinct phases for the plate legend.
  Use real shapes from `src/api.ts` (`Series`, `SeriesMember`, `Sample`) — read them; do not invent fields.

- [ ] **Step 2: Run → FAIL.** `npx vitest run test/print-pages/builderAdapters.test.ts`.

- [ ] **Step 3: Implement `builderAdapters.ts`** — pure functions only; import types from `../../api`, `dominantPhase` from `../waterfall/waterfallModel`. Match the `MemberDatum` interface exported by `MemberList.tsx` exactly (conditional-spread optional fields for `exactOptionalPropertyTypes`). No React, no side effects.

- [ ] **Step 4: Run → PASS.**

- [ ] **Step 5: Gate + commit** — `tsc -p tsconfig.build.json` (0), `lint:design` (0); `git add` the two files; commit (footer `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`).

---

## Task 2: The greenfield builder page (`SeriesBuilderPage.tsx`)

One implementer owns the whole page (it's one coherent component — splitting view/edit/confirm across implementers would thrash a single file). Build it from the wiring table + the composite prop contracts.

**Files:** Create `src/print/pages/SeriesBuilderPage.tsx`; Test `test/print-pages/SeriesBuilderPage.test.tsx`. Reference (read, don't edit): the legacy `src/pages/SeriesBuilderPage.tsx` + `src/components/SeriesRecipeEditor.tsx` (logic-wiring ONLY, never presentation), the page-sim story for the builder if present, `BuilderRail`/`SeriesPlate`/`MemberList`/`AutoGroup` source, `src/print/pages/SeriesFolioPage.tsx` (sibling pattern), `test/AGENTS.md`.

- [ ] **Step 1: Write the failing component test** — mirror the folio/scoping test harness (mock `../../queries`: `useSeries` returns a committed `Series` with `members`; `useSeriesTraces` returns a `Record<number,Trace>`; `useCorpusSamples` returns sample rows; `useSaveSeries`/`useCommitSeriesPlate` return controllable `{ mutate, isSuccess, data, error, isPending }`). Assert the CORE behaviors (add the rest after the page exists):
  - renders the committed series read-state (SeriesPlate title + MemberList rows) with NO draft (`seriesDraft` null), and the "Confirm series" commit is NOT active in read state.
  - editing the title (or adding a sample) **starts a draft** (`startSeriesDraftFromSeries` called) and reflects the edit.
  - **Confirm chain:** with a draft, clicking "Confirm series" calls `save.mutate` with `buildSeriesSaveBody(draft)`; then on `save.isSuccess` (flip the mock + rerender) it calls `commit.mutate` with members from the save response; then on `commit.isSuccess` it `discardSeriesDraft()`. Assert the ORDER (save before commit) and that commit's members come from the save response, not the stale draft.
  - Cancel discards the draft without a request.
  - view controls (offset/scale/grouping) change local state and do NOT start a draft.

- [ ] **Step 2: Run → FAIL** (module not found).

- [ ] **Step 3: Implement the page.** Structure (compose, don't reinvent — every appearance lives in the composites; the page is placement + state):

  **Mount + data:** `const id = Number(useParams().id)`; `const seriesQ = useSeries(id)`; `const tracesQ = useSeriesTraces(id)`; honest `isError`/`isLoading` states (return-early PageFrame + EmptyState on error, Skeleton on load — mirror the folio page). `const series = seriesQ.data`.

  **Draft (lazy):** `const draft = useAppState(s => s.seriesDraft)` + the mutators (selectors). A helper `ensureDraft()` = `if (!draft) startSeriesDraftFromSeries(series!)` called at the top of every recipe-edit handler before the matching mutator. The "effective" title/desc/recipe/order read from `draft ?? series` (draft wins when present).

  **Local view state:** `offset` (default 1.2), `scale` ("log"), `groupingMode` ("bySample"), `collapsed` (false). Annotation toggles via `<AnnotationToggles/>` (self-wires Zustand) — read `showPeakTicks`/`showPeakLabels` from `useAppState` for the export spec.

  **Render model:** `const members = series.members`; `const rows = toWaterfallRows(members, tracesQ.data ?? {})`; `const memberData = membersToMemberData(members)`; pass `offsetScale={offset}` to `SeriesPlate`.

  **Layout:** `PageFrame width="builder"` → the figure plate (`SeriesPlate`, inline-editable title/desc via the draft setters wired into the plate `title`/`subtitle` slots) on the main column + `BuilderRail` aside. `BuilderRail` props: `grouping={<AutoGroup variant="compose" title="Auto-grouped" actions={draft ? [{label:"Confirm series", onClick: onConfirm}, {label:"Cancel", onClick: onCancel, muted:true}] : [{label:"Adjust", onClick: ensureDraft, muted:true}]}>{groupingSummary(series)}</AutoGroup>}`, `orderedBy={effectiveOrderingVariable ?? "—"}`, the order-rule control, `offset`/`onOffsetChange`, `scale`/`onScaleChange`, `traces={<recipe/member rows>}`, `onAddSample`, `onCopyPng` (wire `FigureExportControls`/copy), `onCollapse`.

  **The traces slot:** in read state, render `MemberList` (presentational, hover-linked to the plate via `hoveredKey`/`hoveredQ`). In draft state, render the editable recipe rows (grip/name + ▲/▼ reorder IconButtons → `reorderSeriesSample` by index + dismiss → `removeSeriesSample` by row id) + the "+ Add sample" control (a `Field`/native select of `addableSamples(corpusSamples, draft.recipe)` → `addSeriesSample`). (Reorder is button-based per the carried mutator's index contract; drag is a deferred enhancement — `reorderSeriesSample` takes indices.)

  **The Confirm chain (correctness-critical — implement exactly):**
  ```tsx
  const save = useSaveSeries();
  const commit = useCommitSeriesPlate();
  const stage = useRef<"idle" | "saving" | "committing">("idle");

  const onConfirm = (): void => {
    if (!draft || stage.current !== "idle") return;
    stage.current = "saving";
    save.mutate({ id, ...buildSeriesSaveBody(draft) });
  };
  // Save landed → commit the FRESH plate from the save response (NOT the stale draft).
  useEffect(() => {
    if (stage.current !== "saving" || !save.isSuccess || !save.data) return;
    stage.current = "committing";
    commit.mutate({ id, ...buildSeriesCommitBody(save.data.members) });
  }, [save.isSuccess, save.data, id]);
  // Commit landed → drop the draft; stay on the page (now read state).
  useEffect(() => {
    if (stage.current !== "committing" || !commit.isSuccess) return;
    stage.current = "idle";
    discardSeriesDraft();
  }, [commit.isSuccess]);
  // Either error → reset the chain so the user can retry.
  useEffect(() => {
    if (save.error || commit.error) stage.current = "idle";
  }, [save.error, commit.error]);

  const onCancel = (): void => { stage.current = "idle"; discardSeriesDraft(); };
  const confirmBusy = save.isPending || commit.isPending || stage.current !== "idle";
  ```
  > VERIFY before trusting: that `useSaveSeries` (queue-routed) exposes `.data` as the updated `Series` (with re-resolved `members`). If the queue does NOT surface the response `.data`, fall back to: on `save.isSuccess`, the `useSeries(id)` cache is invalidated → read `seriesQ.data.members` once it reflects the save (track via a members-identity change), then commit. Pick whichever the real `useSaveSeries` supports and note it in the report. The invariant is non-negotiable: **commit the post-save members, never the pre-save plate.**

  **Honest states:** corpus/series load error (return-early), loading skeleton, and a commit/save-error notice (token-classed `role="alert"`). Gate "Confirm series" on `!confirmBusy` and a live draft.

  **Design guard:** `src/print/pages/**` is NOT exempt — placement-only className; appearance lives in the composites. Semantic token classes OK; bracket-literals/raw-hex/side-stripes FAIL.

- [ ] **Step 4: Run → PASS;** add the deferred assertions (read-state has no active Confirm; reorder/remove mutate the draft; export/copy wired; error notice on `commit.error`).

- [ ] **Step 5: Gate + commit** — tsc 0, lint:design 0, `npx vitest run test/print-pages/SeriesBuilderPage.test.tsx`. Commit (co-author footer).

---

## Task 3: Cutover — repoint route + grep-delete the legacy builder

**Files:** Modify `src/components/AppRoutes.tsx`; delete the OLD set (grep-enumerated).

- [ ] **Step 1: Repoint** — in `AppRoutes.tsx`, change the `/series/:id` import from `../pages/SeriesBuilderPage` to `../print/pages/SeriesBuilderPage`. Leave the Route line.

- [ ] **Step 2: Grep-enumerate the OLD set** (don't trust the ledger):
  ```bash
  grep -nE 'from "\.\./components/(SeriesRecipeEditor|SeriesBuilderRail|SeriesMemberRow|SeriesReadingPanel|MultiTracePlot|MemberMetaGutter|OffsetDock)' src/pages/SeriesBuilderPage.tsx
  for f in SeriesRecipeEditor SeriesBuilderRail SeriesMemberRow SeriesReadingPanel MultiTracePlot MemberMetaGutter OffsetDock; do
    echo "== $f =="; grep -rln "components/$f\"" src --include=*.tsx | grep -v "src/pages/SeriesBuilderPage.tsx"
  done
  ```
  Delete only confirmed orphans. CRITICAL: `MultiTracePlot` may be imported by other surfaces (it was the shared Compare/Series render core) — if grep shows a surviving consumer, LEAVE it and report. Same caution for `SeriesMemberRow` (a greenfield `SeriesMemberRow` ALSO exists under `src/print/` — only delete the `src/components/` one).

- [ ] **Step 3: Delete** the legacy page + confirmed-orphan components + their tests via `git rm`. KEEP `GroupingModeToggle`/`AnnotationToggles`/`FigureExportControls`/`SeriesPhaseStrip` (greenfield consumes them).

- [ ] **Step 4: Guards**
  ```bash
  grep -rnE 'from "(\.\./)+(pages|components)/(Series(Builder|RecipeEditor|Reading)|MultiTracePlot|MemberMetaGutter|OffsetDock)' src/print/pages/SeriesBuilderPage.tsx && echo LEAK || echo clean
  grep -rn "SeriesRecipeEditor\|SeriesBuilderRail\|MemberMetaGutter\|OffsetDock" src test && echo DANGLING || echo clean
  ```
  Both `clean` (modulo legitimate greenfield imports — read the page's import block by eye).

- [ ] **Step 5: FULL-suite gate** — `tsc` 0, `lint:design` 0, `npm test` (WHOLE vitest — green; migrate any test referencing a deleted component or legacy-only testid to the greenfield DOM), `npm run build`. Do not leave the suite red.

- [ ] **Step 6: Commit** (co-author footer; `git add` the route + migrated tests; `git rm` already staged deletions).

---

## Task 4: Migrate the builder e2e spec

**Files:** Modify `e2e/series-builder.spec.ts` (exists). Reference `e2e/series-folio.spec.ts` + `e2e/series-scoping.spec.ts` (migrated siblings) + `e2e/AGENTS.md`.

- [ ] **Step 1: Read the current spec + siblings** — map legacy selectors (`recipe-title`/`recipe-add-sample`/`recipe-save`/`recipe-commit`/`recipe-cancel`) to the greenfield DOM.
- [ ] **Step 2: Rewrite** the route mocks (`**/api/series/:id`, `**/api/series/:id/traces`, `**/api/samples` corpus, `**/api/series/:id` PATCH save, `**/api/series/:id/commit`, drain SSE + identity seed) and the flow for the lazy-draft + Confirm-chain DOM: edit (starts draft), add a sample, **Confirm series** → assert a PATCH `/api/series/:id` fires THEN a POST `/api/series/:id/commit` fires (order), and the page returns to read state. There is no separate Save button and no conflict modal.
- [ ] **Step 3: Run** `npx playwright test e2e/series-builder.spec.ts` → green. Update `e2e/smoke.spec.ts` if it visits `/series/:id`.
- [ ] **Step 4: Commit** (co-author footer).

---

## Final gate (orchestrator)
- [ ] `tsc -p tsconfig.build.json` (0) · `lint:design` (0) · full `npm test` (green) · `npm run build` · `npx playwright test e2e/series-builder.spec.ts e2e/smoke.spec.ts`.
- [ ] Final `frontend-reviewer` over the whole cutover — pay special attention to the **Confirm Save→Commit chain** (no stale-plate commit; pending-ref reset on error; commit disabled during save) and the lazy-draft trigger (view controls never start a draft).
- [ ] Do NOT run `finishing-a-development-branch`. After this lands, the three Series pages are all greenfield — the legacy `src/pages/Series*` are gone and the branch is one cutover (dead-code track 7 + shell reskin 8) from mergeable.

## Deferred (out of scope for 6b)
- Drag-to-reorder recipe rows (the carried `reorderSeriesSample` is index-based; button reorder ships, drag is an enhancement).
- Representation (waterfall/heatmap) toggle + Track-reflections (HeatmapChart ⛔ / TrackingLine ⏸).
- The vertical `#phasestrip` companion (needs an L1 renderer on WaterfallChart row geometry).
- Boneyard `builder.bones.json` capture + live visual-fidelity pass (need a live dev server).
