# Phase-4 Cutover Strategy — greenfield "The Print" HimalayaUI rebuild

**Date:** 2026-06-03
**Branch:** `worktree-greenfield-ui-rebuild` (unmerged, stays unmerged until the whole rebuild is complete)
**Evidence base:** [`docs/greenfield-phase4-cutover-ledger.md`](../../greenfield-phase4-cutover-ledger.md) — the 12-area survey (carry/drop/new/refactor/gap) this strategy is built on.

## Purpose

Layers 0–3 of the greenfield rebuild are complete: every "The Print" composite is built and assemblable in Storybook, plus the figure-export control. What does **not** exist is the page layer — the six route bodies that wire those composites to live data. Phase 4 is the **cutover**: replace each legacy page body with its greenfield equivalent, route by route, until the live app runs entirely on The Print.

This document fixes the *strategy* — scope, the new-vs-old discipline, the cutover model, the repeatable per-route recipe, the build sequence, and the carried risks. It is the input to the implementation plan(s) that follow.

## Scope & non-goals

**In scope:**
- Swap all six page bodies legacy → greenfield, one route at a time.
- One new backend read route: `GET /api/series/:id/traces` (gates the folio).
- Relax the series-commit 409 optimistic-concurrency gate + drop the conflict modal.
- Per-route deletion of the OLD code each swap replaces.
- Final shell reskin once all six bodies are greenfield.

**Non-goals (each its own downstream sub-project, brainstormed separately):**
- **Edit-tracking → undo/redo → versioning.** This replaces the *conflict model*, not any page swap. The cutover ships on last-write-wins (LWW) — exactly what is live today, so no regression. Undo/redo gets its own spec after the app is on The Print.
- Keyboard-reorder a11y for the drag lists (pointer-only today).
- The optional `GET /api/samples/bulk-exposures` endpoint (the `useCorpusExposures` fan-out is bounded by `useQueries`; not a blocker).

## Provenance model — the "deliberate new vs old" rule

The central discipline. At every commit, every file is unambiguously in one of **three named buckets** — never a two-way old/new blur:

- **NEW** — `src/print/**`. The greenfield pages and composites.
- **OLD** — `src/pages/**`, `src/components/**` (legacy page bodies + presentational components). **Deleted per-route as each is swapped out** — never left wired-but-dead.
- **CARRIED** — provenance-neutral infra: `queries.ts` / `api.ts`, `lib/queue/**`, the SSE mount, the custom hooks. This is *not* "old code we tolerate" — it is code that **already is what the greenfield implementation would write**.

**The repurpose test:** a mechanism may be CARRIED only if the answer to *"would we write it this way in the new surface?"* is **yes**. The data plane, mutation queue, SSE wiring, and hooks pass — they are plumbing, identical either way. Anything that fails the test is OLD and gets rebuilt in `print/`, not imported across the line.

**Two hard boundary rules** (keep the line from blurring during the mixed-state window):
1. **`src/print/pages/**` imports NEW + CARRIED only — never OLD.** If a swapped page reaches for a legacy presentational component, that is a *gap to build in `print/`*, not a cross-boundary import. This is mechanically checkable: a grep for legacy-component imports inside `src/print/pages/` must return nothing.
2. **No partial-greenfielding of a legacy page.** We never edit a legacy `src/pages/*.tsx` to "start using print bits." Each swap is atomic: the new body lands in `src/print/pages/`, the route repoints, the old body is deleted.

The result is a genuinely new rebuilt surface that *repurposes* existing mechanisms — not a hybrid of both systems.

## Cutover model — incremental-by-route under the carried shell

The ledger frames "no greenfield `App.tsx`/`AppRoutes`/`CorpusShell`" as a gap implying a parallel greenfield app flipped big-bang. That gap dissolves: the SSE lifecycle, the mutation queue, and the QueryClient are all CARRIED — the *mechanisms* survive untouched; they merely live in the legacy shell files today.

So:
1. Keep the existing `App.tsx` + `AppRoutes` + SSE mount + queue + `CorpusShell` as-is.
2. Swap one route's page body legacy → greenfield at a time. `AppRoutes` routes by path, so a greenfield page slots in wherever the legacy one was; navigation between mixed pages just works (it is route-based).
3. Greenfield bodies mount **body-only** inside the existing `CorpusShell` `<Outlet>` (no chrome of their own during migration).
4. After each swap the whole app still runs — visually check that one page against its mockup, recapture its skeleton, commit. A clean checkpoint.
5. **Reskin the shell last** — `CorpusShell` / `CorpusTopbar` / `NavModal` / etc. are presentational; swap them in one deliberate step once all six bodies are greenfield.

This satisfies the ledger's hard constraint (atomic per route) while keeping a runnable, demoable app at every step and never producing a giant final-merge diff. The only cost is a transient mixed-look app mid-migration — seen only by the developer, on an unmerged branch.

**Cutover order — Strategy 1 (harness-first → fail-fast → door → cluster):**

`Loupe → Focus → Samples → Series (folio → scoping → builder)`

Rationale: prove the swap *mechanics* on the smallest page first (Loupe), then deliberately hit the highest-uncertainty integration early (Focus: the TracePlot engine swap + assignment model) so problems surface before five more pages are invested — and because the engine work de-risks the Series waterfall. Samples (the door) next. The Series trio last as a cluster, because that is where both backend changes land.

## Per-route recipe (the repeatable unit of work)

Each route is one self-contained, just-in-time unit. Build its gaps *with* the page that consumes them (not as an upfront batch), ending at a runnable, gated checkpoint:

1. **Build the gaps** — the route's missing `print/` composites. TDD, design-guard-clean (placement-only; appearance in `ui/` primitives), Storybook stories. Each composite its own commit.
2. **Assemble the page** — `src/print/pages/<Name>Page.tsx` composes the route's composites into the full body (the shell-level composition the assembly stories deferred). Imports NEW + CARRIED only.
3. **Wire CARRIED data** — connect to the existing query hooks / queue / SSE / route params. The page owns the cross-child state the panels expect (selection, active band, …). Add the page-level pure adapters the ledger named (`Exposure → GalleryExposure` URL adapter, `WaterfallRow` adapter) as small helpers.
4. **Repoint the route** — swap the route's element legacy → greenfield in `AppRoutes`. App still runs (mixed state).
5. **Delete the OLD** — remove the legacy page body + its now-orphaned legacy-only components (grep confirms no remaining importers), respecting the ledger's *extract-first* traps. CARRIED is untouched.
6. **Recapture boneyard** — capture the new page's `bones.json` so the loading skeleton matches the new layout.
7. **Gate & commit** — **visual check vs mockup** (developer eyeball of the rendered page against the relevant `docs/redesign-mockups/*.html`; there is no automated fidelity harness in this branch — an automated screenshot gate is a possible later sub-project, not a cutover prerequisite), unit tests, e2e (mocked), `npm run lint:design`, `npx tsc --noEmit -p tsconfig.build.json`, `npm run build`. Commit the swap.

## Build sequence

1. **Loupe** (`/samples/loupe/:sampleId`) — prove the swap mechanics on the smallest, read-mostly page. Also establishes `src/print/pages/` + the body-in-`CorpusShell`-outlet mount convention. Gaps: small ports (`LoupeSidePanel` tags editor, reject chips), greenfield `ConflictModalShell`/`ScopingConfirmModal` ports if reached.
2. **Focus** (`/sample/:sampleId`) — fail-fast on the riskiest integration: the TracePlot **Observable Plot → d3** engine swap + the durable assignment model. Gaps: `StaleIndicesBanner` (hash-stale debounce), `SpeculativeBuilder` (anchor + ratio dialog), `FocusNotesMargin` (q-ref notes + focus-gated textarea + drawer), `losingPeakIds` dim, q-link/`onHoverQ` hover circuit.
3. **Samples** (`/samples`) — the door. Gaps: corpus → focus door (`StatusCell`), gallery dbl-click / shift-click → loupe, X/Esc cull shortcuts, `CorpusSample.screened`/`.phase` field wiring.
4. **Backend:** implement `GET /api/series/:id/traces` (a single `exposure_id → Trace` map per series; Julia, TDD) + `useSeriesTraces` hook + the `CardFigure`/`WaterfallRow` adapter. **Then Series-folio** (`/series`). Gaps: series-pins wiring (`useSeriesPins`/`usePinSeries`/`useUnpinSeries` over the already-implemented routes; folio pin affordance), `useSeriesForksOf` hook.
5. **Series-scoping** (`/series/new`) — gaps: greenfield scoping rows/sparkline/order/confirm composites; `FlagButton`. **Decide the parked product questions here** (see below).
6. **Backend:** relax the series-commit 409 gate — delete the `expected_hash` extraction + the 409 response branch + the `current_series_content_hash` read in `routes_series.jl:225-244`; **keep** `compute_series_content_hash`/`current_series_content_hash` (still write `content_hash` for the `series_plate_committed` post_state + future fork/stale checks). Do this **atomically** with dropping `SeriesCommitConflictModal`/`ConflictModalShell` and the frontend ceasing to send/handle `expected_content_hash`. **Then Series-builder** (`/series/:id`). Gap: `SeriesRecipeEditor` (the largest single behavioral gap — title/desc/add-sample/order/Save/Commit wired to `useSaveSeries` + `useCommitSeriesPlate`; `BuilderRail` is only a presentational traces slot today).
7. **Backend dead-code track** — relocate the shared `comparisons.jl` helpers (`canonical_json`, `comparison_now_iso`, `compute_member_snapshot`, `is_member_stale`, `recently_used_exposures`, `picker_samples`, `_topk_phases`, …) to a neutral helpers file **before** any deletion (Julia include-order is load-bearing); then delete `routes_messages.jl` + its registration/include/test and the `series/:id/messages` route registrations. **Keep** `comparison_*` tables + the comparison dispatcher branches in `events.jl` — `rebuild_views_from_log!` re-folds historical comparison events, so dropping them breaks replay on any DB with history. On the frontend, drop the retired comparison/message hooks + their tests, bump `SCHEMA_VERSION` so stale persisted ops are ignored, and remove the comparison SSE arms after the mutators are gone.
8. **Shell reskin** — swap `CorpusShell` / `CorpusTopbar` / `NavModal` / `OnboardingFlow` / `InfrastructureBanner` / toast container onto print primitives; wire `CorpusTopbar onChange → navigate` (derive active from pathname). Delete the remaining OLD shell. Full-suite gate. The app is now entirely greenfield.

Backend changes attach to their gating route (traces → folio; 409 relax → builder); the dead-code track runs once the consuming frontend hooks are gone.

## Parked product decisions (decided at the Series stage, not now)

- **Series-commit gate:** the strategy relaxes it (step 6) to fulfill the no-conflict-UI decision. Confirm at the builder stage that LWW-on-commit is acceptable, or keep the gate.
- **Experiment-scoped picker + "recently used":** keep the `usePickerSamples`/`useSampleTags`/`useRecentlyPickedExposures` hooks + experiment-scoped picker routes (future builder picker) or drop as currently-unused? Decide when wiring the builder picker. (Corpus-wide `GET /api/sample-tags` + `GET /api/picker-samples` are unambiguously carry — used by scoping.)
- **`GroupingModeToggle` / `AnnotationToggles`:** keep-vs-drop for the builder grouping mode. Decide at the builder stage.

## Risks carried from the ledger (honored in execution)

- **TracePlot engine swap** (Observable Plot → d3) — deep entanglement: peak hit-testing, scroll-zoom, label placement, `losingPeakIds` inversion. Surfaces early by design (Focus is route 2).
- **`CardFigure` N+1** — resolved by the traces route gating folio; otherwise O(N×M) fetches on folio load.
- **`comparisons.jl` include-order** — Julia loads it before `series.jl`; relocate shared helpers before deleting or get `UndefVarError` at module load.
- **Keep `comparison_*` tables + dispatcher branches** — required for `rebuild_views_from_log!` replay exhaustiveness; every new event kind also needs a dispatcher branch + round-trip test.
- **409 relax must be atomic** with the modal drop + the frontend ceasing to handle 409, or transition-window `ConflictError`s go unhandled.
- **SSE-after-commit ordering** — broadcast fires after the SQLite transaction commits (`_flush_post_commit_broadcasts!`).
- **`useCorpusExposures` O(1)-observer** — must not regress to per-row `useExposures` (caused the `ERR_INSUFFICIENT_RESOURCES` bug on the 139-sample corpus).
- **Focus notes cache-coherence** — read `useSamples(experimentId)`, not `useCorpusSamples`, or notes silently revert after save.
- **`defaultExposureId` ordering** (rep → first-accepted → first) — load-bearing for which frame opens; not caught by render tests.
- **`StageTabs`/topbar must `navigate`, not `setActiveSample`** — or the `,`/`.` step shortcuts break on `/sample/:id`.
- **Boneyard recapture before each ship** — `bones.json` must be (re)captured + committed before a greenfield page ships, or the skeleton falls back to plain text.
- **`SCHEMA_VERSION` / persist version bumps** — bump `sessionStorage` `SCHEMA_VERSION` when OpKinds change; bump Zustand persist version + migrate when adding persisted fields.

## Success criteria

- All six routes render greenfield bodies.
- `src/pages/**` and the legacy presentational components are deleted; the shell is reskinned.
- **Zero** legacy imports inside `src/print/**` (grep-clean).
- Full suite green: Julia (HimalayaUI), Vitest, Playwright (mocked) E2E; each surface visually checked against its mockup at swap time.
- Dead message/comparison **routes** removed; comparison **tables + dispatcher branches** retained for replay.
- The app runs on last-write-wins with **no** conflict UI.

## Decomposition note

This strategy spans six page swaps + two backend changes + a dead-code track + a shell reskin. It is large but cohesive — each route is an independently-testable unit sharing one recipe. The implementation plan may be split per stage (e.g. one plan for Loupe+Focus+Samples, one for the Series cluster + its backend changes, one for the dead-code track + shell reskin) so each plan produces a runnable, gated checkpoint. The undo/redo/versioning sub-project is explicitly **out** and gets its own spec → plan cycle later.
