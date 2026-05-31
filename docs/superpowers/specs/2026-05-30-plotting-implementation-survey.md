# Plotting redesign — implementation survey & architecture decisions

> **Status:** survey complete 2026-05-30. Two multi-agent codebase surveys (backend build-readiness, ~620K tokens / 7 agents; frontend integration-risk, ~776K tokens / 8 agents) plus three settled architecture decisions. This is the **input to `writing-plans`** — the plan is downstream of this doc. Companion to `2026-05-29-plotting-redesign-design.md` (the design spec) and the two mockups.

## Executive summary

The redesign is **mostly a frontend build on data the backend already computes**, plus a small number of net-new backend items and a clean rebuild of the plot layer. Two surveys established this:

- **Backend ~70% ready.** Ranked candidates with score/basis/`lattice_d`/phase, full `predicted_q` (all orders incl. absent), per-reflection `residual`+`q_observed`, `ngc`, the speculative builder, and waterfall member layout are **already computed and already serialized** to the client (`GET /api/exposures/{id}/indices`, verified `routes_analysis.jl:111-136`). The "physics-heavy" features (combs, residuals, q-linked rings) need **zero backend work** — only re-presentation.
- **Frontend ~40% surgery / ~60% additive.** The deep change is the assignment cart replacing the single-active model; the rest (peakMark, PlotSurface, rings, comb panel, preview, custom-index modal, series views) is additive *provided the cart's state model + SSE shape are designed first*.

## Settled architecture decisions (2026-05-30)

1. **One durable assignment per exposure — retire the multi-group machinery.** The cart **is** the assignment: 0..N phase-indices as members (each carries its own lattice/peaks, so coexistence = 2 members), plus an explicit `state: 'indexed' | 'form_factor' | 'null'`. Candidates remain a separate ranked list. *We are free to simplify the backend* — retire the auto-group / custom-group / `idx_one_custom_group_per_exposure` complexity rather than work around it (user: "simpler is better, don't keep guessing multiple groupings"). Less surgery **and** less standing complexity.
2. **Three distinct states, not a nullable phase.** `indexed` (≥1 phase) · `form_factor` (structured scattering, no lattice — a positive observation worth showing) · `null` (no interesting scattering — flat, featureless). Form-factor and null are **different** and render/read differently. Explicit enum on the wire, **never inferred** from `members.length` (inference is fragile under SSE replay).
3. **Cart is durable backend state.** The committed cart persists in SQLite, round-trips the event log, loads from the API on `/sample/:id`. **No URL serialization** (dodges the deliberate one-way URL→store contract / handshake race). The hypothetical hover **preview** stays ephemeral Zustand, orthogonal to the cart, never a mutator or SSE event.
4. **The plot layer is a greenfield rebuild, not a preserve-refactor.** The existing custom SVG overlay state machine is **not sacred** (user: "we're not beholden to keeping the plot infrastructure… if it's not necessary we don't have to keep it"). PlotSurface owns peak drag/snap/hit-test in whatever architecture is cleanest. We preserve the **behaviors** (q-link feel + tolerances — don't regress #180; drag-to-add/move/exclude; snap-to-peak), not the **code**.

## Backend gap matrix

| Feature | Status | Backend need |
|---|---|---|
| F2 ranked candidates (score + lattice a/d + phase) | **ready** | none — already sorted, serialized |
| F4/F5 combs + indexing-space residual | **ready** | none — `predicted_q` (all orders) + `index_peaks.residual`/`q_observed` already shipped |
| F6/F9 detector rings + q-link | **ready** | none — q values already flow through `peaks`/`predicted_q` |
| S1 waterfall/heatmap; S2 migration; S7 export | **ready** | none — member layout (`y_offset`/`band_height`) + trace arrays served |
| S4/S5 series analytics (phases-present, member rows) | **partial** | client-derivable from member snapshots; *optional* `q1` snapshot field or aggregation endpoint to avoid per-exposure fan-out |
| **B1 assignment container** (the cart) | **rework** | Replace single-active group model with one durable assignment per exposure (decision #1). Touches `db.jl` (schema + retire `idx_one_custom_group_per_exposure:130`), `events.jl` (index_confirmed/unconfirmed `344-358`), `routes_analysis.jl` (`ensure_custom_group!` `14-40`, member routes `146-229`). **Needs schema migration + `rebuild_views_from_log!` round-trip tests.** |
| **B2 form-factor / null states** | **missing** | `indices.phase` is `TEXT NOT NULL` (`db.jl:75`). Add the 3-state enum to the assignment (NOT a nullable phase); ripples through `effective_peaks`, `persist_analysis`, `series.jl` snapshot composition, every phase-string consumer |
| **B3 Gauss-Bonnet predictor** | **missing** | Self-contained, **no schema change**. Given an assigned bicontinuous cubic's lattice, predict the coupled cubic (a_Im3m≈1.279·a_Pn3m), flag/boost matching candidates. Today only observational κ exists (`_ngc_for_phase:94-109`), no ratio coupling, no score term |
| **B4 lattice-driven speculative build** | **partial** | Custom-index builder picks symmetry+lattice, computes reflections client-side; existing route is anchor-peak driven (`routes_analysis.jl:297+`) — persistence must accept a client-fitted basis |
| **B5 assignment event kinds** | **missing** | `assignment_add` / `assignment_remove` / `set_form_factor` / `set_null` following the `apply_event!` dispatcher contract (sole writer to view tables), with idempotency + SSE |

## Frontend: must-preserve invariants (cross-cutting, survive the plot rebuild)

These are **data/build layer**, unaffected by the plot freedom:

- **`client_op_id` minted inside `mutationFn`**, never at hook construction (`useQueueMutation.ts:126`) — else retries share an id and break backend idempotency (cart stuck at N−1).
- **Foreign-event replay** = rollback pending reverse, re-apply insertion order, per-iteration try/catch (`replayCoordinator.ts:76-116`). New cart ops auto-participate.
- **Negative optimistic ids = the confirmed/pending flag** — `peakMark()` renders `id<0` outline-only; no consumer filters `id>0` or derefs `/api/peaks/${negative}`.
- **`peakMark()` colour-parameterized** (no CSS-var reads) — export renders literal `#hex` on canvas (`presets.ts:28-40`); on-screen passes resolved colour, export passes hex.
- **`check-design` 0-baseline gate** — new phase-colour files pass only if colour threads through JS mark options (like TraceViewer does), not `style`/Tailwind; else allowlist (`check-design.mjs:37-54`) or compose a `ui/` primitive.
- **SSE post_state shape is context-dependent** (`applyRemoteToCache.ts:18-41` casts `CurationPostState`, guards `Array.isArray(ps.indices)`). Assignment events MUST define a **distinct** post_state shape or the cast reads `undefined` and clobbers cache.
- **URL→store sync is one-way** (`useSyncActiveSampleFromRoute`) — no store→URL reverse. `/sample/:id` is the permalink; the cart loads from the API (decision #3).

## Frontend: HIGH-severity conflicts → resolutions

1. **Cart vs `GroupEntry.active`** (3 readers: `PhasePanel.tsx:14`, `PlotCard.tsx:231`, `FocusReflectionsTable.tsx:45` do `find(g=>g.active)`) → explicit `state` discriminator separate from phases array; one `deriveActiveIndices(assignment)` helper replaces all three reads.
2. **SSE wipes new fields** (`members` splice spreads old shape, `applyRemoteToCache.ts:165/178`) → new assignment case arms (or wholesale `invalidate(groups)`) + distinct post_state shape; register all 4 kinds in **both** `resolveMutator` and `resolveMutatorForEvent` (`mutatorRegistry.ts:58-206`).
3. **peakMark vs 5 drifted defs** (TraceViewer SVG polygon `PEAK_OFFSET_PX=7`; MemberTraceLayer `symbol='triangle'` offset 5; traceExportMarks `'triangle2'`; PlotCard `TriangleSvg`) → one colour-parameterized builder returning both `Plot.Markish` and SVG geometry; converge all 5 + retire `--color-peak-manual` magenta in **one sweep** (incl. `traceAdapter.test.ts:59` manual-legend assertion → diamond). `offsetPx` parameterized (7/5) to avoid the visual jump.
4. **New colour files vs lint:design** → thread `phaseColor()` into JS mark options only, or allowlist, or compose `PhaseStrip`/`PhaseChip`/`ScoreBar`.
5. **Series fan-out perf cliff** (~62 round-trips for a 20-member series; S5 adds `useMemberIndices`) → keep intentional fan-out, bound concurrency like `useCorpusExposures` (scheduled `useQueries`, `queries.ts:198-214`); don't snapshot-embed (de-syncs from live curation).

## Reuse inventory (most of the redesign rides existing infra)

- `useQueueMutation` framework — 4 new cart mutators just implement `kind`/`onMutate`/`request`/`onSuccess`/`synthesizeFromSse`.
- `nextOptimisticId` + `replacePlaceholder` — optimistic cart blocks.
- `resolveMutator`/`resolveMutatorForEvent` dispatch tables — register the 4 kinds.
- `phaseColor()` (`phases.ts:46`) + `colorFor()` (`lib/comparison/coloring.ts`) — single colour source for peakMark/rings/teeth/strip/blocks.
- `ui/` primitives: `PhaseChip`→cart blocks, `PhaseStrip`→series strip, `SegmentedControl`+`ModalShell`→custom-index modal, `ScoreBar`→quality.
- Plot helpers (now **as references for the rebuild**, not preserve targets): `nearestClickablePeak` (filters `id<0`), `invertQ`/`applyQ`, `makeXScale`, `computeYBands`, `buildMemberMarks`.
- Figure-export `traceAdapter`/`multiTraceAdapter` + `LIGHT_PALETTE` + `ExportSpec` thunk — S7 clean preset.
- `hoveredQ` channel (`state.ts:73-78`) + `DetectorRingOverlay` sink — comb-tooth as third q-link node.
- `losingPeakIds` (`PlotCard.tsx:246-259`) — F7 preview dim/orphan logic (rewrite `alreadyActive`→"in any cart block").
- Six-layer contract-test templates (`test/AGENTS.md`, `mutation-queue.md §11`) — every new cart kind.

## Build order (backend + frontend interleaved, dependency-respecting)

1. **Cart state model + SSE post_state shape** (design contract first — everything F1/F6/F7 depends on it).
2. **B1+B2+B5 backend**: durable single-assignment schema + 3-state enum + 4 event kinds, with `rebuild_views_from_log!` round-trip tests. Gates all Focus cart / custom-index UI and S6.
3. **B3 Gauss-Bonnet** (self-contained, no schema — can run parallel to 2).
4. **X1 peakMark()** — converge all 5 sites + retire magenta in one sweep. Colour-parameterized.
5. **X2 PlotSurface** — greenfield rebuild owning the Plot instance + scales + gestures + rAF resize + peak interaction (drag/snap/hit-test) + `overlay`/`hitTest` API. Refactor TraceViewer/MillerPlot/MultiTracePlot onto it. Preserve behaviors, not the old overlay.
6. **Focus F4/F5 combs+residual** and **F6/F9 phase-colour rings + comb-tooth q-link** (pure frontend on shipped `predicted_q`/`residual`, once peakMark+PlotSurface exist).
7. **Focus F1 cart UI + mutators** and **F7 preview** (needs the cart contract + 4 kinds).
8. **Focus F8 custom-index modal** (needs B4 lattice-build path).
9. **Series S1/S2** on PlotSurface; **S3/S4/S5 + S6 form-factor/null members** (needs B2).
10. **Series S7 clean export** last, reusing peakMark at the export boundary.

## Open questions for the plan

- **B4 commit path:** does a fresh custom index route through the existing `createSpeculativeMutator` (then translate to assignment-add) or a new `assignmentAddSpeculativePhase` kind? Watch issue #37 Bug 1 (`ensure_custom_group!` mints a fresh id not yet in cache — invalidate, don't splice).
- **Q-link tolerance unification:** trace uses per-peak relative (`Q_LINK_REL_TOL=0.01`); rings use span-relative (2% clamped). The comb-tooth third node + a unified formula: span-relative (consistency) or per-peak (feel near small q)? Design call — re-tune deliberately, don't regress #180.
- **Bonnet recompute (F3):** confirm never persisted; on cross-tab assignment change, does `applyRemoteToCache` invalidate `queryKeys.indices` so rankings recompute, or recompute every render?
- **S5 perf:** confirm per-exposure `useMemberIndices` fan-out (fresh, slower) over snapshot-embed (stale, faster) — or add a batch endpoint?
- **Spec reconcile:** §8.5 still recommends keeping manual-peak magenta; §5.1 retires it (gray diamond). Flip §8.5 so `peakMark()` consumers don't fork.

## Source

Raw survey syntheses (both `{surveys, synthesis}` structured outputs) were produced by the `plotting-impl-scope` (run `wf_9d222d68-ab7`) and `frontend-integration-survey` (run `wf_62aa66b8-414`) workflows. All file:line claims spot-verified against source before this doc was written.
