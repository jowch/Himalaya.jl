# Plan E — Series Surface Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the redesigned Series surface — a waterfall/heatmap thinking-view, light anchors as migration-tracking handles keyed by `(phase, reflection order)` with hollow ghost rings where an order is absent, a phase-strip companion, a derived "phases present" reading, member rows showing lattice (a/d) + first-peak q₁, form-factor/null members handled end-to-end, and a decoupled clean scientific export preset.

**Architecture:** Most of Series is **client-derivable** from data the member snapshot already carries (`confirmed_index{phase, lattice_d, peak_ids}` + `effective_peaks`). The existing `MultiTracePlot`/`MemberTraceLayer`/`MemberHeatmapLayer`/`CrossTraceTrackingLayer` (which already keys tracking by `(phase, Miller-order)`) are refactored onto `PlotSurface` (Plan C) and extended: the tracking layer gains the hover-anchor handle affordance + hollow ghost ring at predicted-q where an order is absent; a new `SeriesPhaseStrip` and `seriesReading` derivation read per-member snapshots; member rows extend `MemberMetaRow` with lattice + q₁; the export adapter gains a clean preset. The per-exposure fan-out perf cliff is bounded via the existing `useStableQueryMap`/scheduled-`useQueries` pattern — snapshots are NOT embedded (they'd de-sync from live curation).

**Tech Stack:** TypeScript, React, TanStack Query (`useStableQueryMap`), Observable Plot via `PlotSurface` + `peakMark`, the `ui/` system (`PhaseStrip`), the figure-export adapter pipeline, Vitest, the mockup `docs/redesign-mockups/2026-05-29-series-plot.html` as the pixel/interaction spec.

---

## Dependencies & scope

- **Depends on Plan C** (`PlotSurface`, `peakMark`) and **Plan A's B2 three-state** (form-factor/null members must round-trip through the member snapshot — the series member's `confirmed_index` is null and the assignment state distinguishes form_factor vs null).
- Independent of Plan D's cart UI, but reuses the 3-state vocabulary.
- The mockup is the visual/interaction contract: light plate-ringed anchor beads, terracotta migration connector on hover, hollow ghost ring for absent orders, two-phase gradient phase-strip cells, the derived reading lines, the "no Bragg peaks · q₁ —" form-factor row, the clean export (white, 2px traces, Arial axes/title/footnote).

## What the snapshot already gives us (verified, no backend needed)

`SeriesMember.snapshot` → `{ effective_peaks: [{id,q,intensity,sharpness,source}], confirmed_index: {id, phase, lattice_d, r_squared, ngc, peak_ids[]} | null, analysis_inputs_hash }`. So per member, client-side: phase (or null), lattice `a/d` (= `lattice_d`), first-peak q₁ (= `min(q)` of effective_peaks), and `(phase, order)` keying (via `peak_ids` Miller order). The derived reading + member rows + migration tracking all flow from this. **Optional backend (S5):** add `q1` to the snapshot only if `min(effective_peaks.q)` proves too costly client-side — default is client-derive.

## File structure

| File | Change | Responsibility |
|---|---|---|
| `src/components/MultiTracePlot.tsx` | Refactor onto PlotSurface (Plan C) | hero; waterfall/heatmap; compose layers |
| `src/components/CrossTraceTrackingLayer.tsx` | Modify | hover-anchor handle affordance + hollow ghost ring at predicted-q for absent orders; keep `(phase, k)` keying |
| `src/lib/series/seriesReading.ts` | **Create** | derive "phases present": per-phase variable span + lattice trend, coexistence/form-factor lines |
| `src/components/SeriesPhaseStrip.tsx` | **Create** | one cell per sample along the variable; coexistence gradient; hollow dashed form-factor cell |
| `src/components/SeriesMemberRow.tsx` | **Create** (or extend `MemberMetaRow`) | lattice a/d (both under coexistence) + first-peak q₁, phase-coloured; "no Bragg peaks · q₁ —" for form-factor/null |
| `src/lib/series/anchors.ts` | **Create** | `(phase, order)` anchor map + absent-order detection from snapshots |
| `src/lib/figure-export/presets.ts` | Modify | `CLEAN_SCIENTIFIC` preset (white, 2px, Arial) |
| `src/lib/figure-export/adapters/multiTraceAdapter.ts` | Modify | clean-preset path; include form-factor members; use `peakMark` |
| `src/pages/SeriesBuilderPage.tsx` | Modify | mount phase-strip + reading + member rows; keep bounded fan-out |
| tests | Create/modify | derivation + layer + export tests |

---

## Task sequence

### Task E-1: Refactor MultiTracePlot onto PlotSurface

**Files:** `MultiTracePlot.tsx`. (Overlaps Plan C Task 5 — if done there, this is a no-op verification.) Preserve `Representation` waterfall/heatmap, `computeYBands`, the layer composition, `offsetToBandFraction`/`computeMaxPlotWidth` exports. Peaks via `peakMark`.

- [ ] TDD: keep existing MultiTracePlot tests green on PlotSurface; waterfall + heatmap both render.
- [ ] Commit: `refactor(series): MultiTracePlot on PlotSurface` (skip if landed in Plan C).

### Task E-2: `(phase, order)` anchor map + absent-order ghost rings

**Files:** Create `src/lib/series/anchors.ts`; modify `CrossTraceTrackingLayer.tsx`. The layer already keys `(phase, k)` (Miller order) and draws polylines. Add: (a) `buildAnchorMap(members)` → `Map<"phase:order", Array<{memberPos, q | null}>>` where `q=null` means the phase is present but that order is unobserved (absent); (b) render a hollow ghost ring (`peakMark` predictedAbsent) at the predicted-q on that member's baseline and route the connector through it; (c) where the phase is absent entirely, the track ends (no vertex). Predicted-q for an absent order comes from the member's `confirmed_index` basis × phase ratios (same physics as Focus combs).

- [ ] TDD: `buildAnchorMap` keys by phase+order; a member with the phase but missing order k yields a `q=null` (absent) vertex; the layer draws a ghost ring there and routes through it; a member without the phase ends the track. Verify against the mockup's 1:0.25 Pn3m √3-absent case.
- [ ] Commit: `feat(series): (phase,order) migration tracking with absent-order ghost rings`.

### Task E-3: Hover-anchor handle affordance

**Files:** `CrossTraceTrackingLayer.tsx` + `MultiTracePlot` overlay. Anchors are light plate-ringed beads; hovering one tracks that `(phase, order)` across the whole stack (terracotta connector). Wire via `PlotSurface.overlay` hit-testing on anchor beads (reuse `hitTest`); the hover state is ephemeral (Zustand or local). a11y: a focus path to the same highlight (mockup audit P1 — tracking is currently hover-only).

- [ ] TDD: hovering an anchor sets the tracked `(phase,order)`; the connector threads every carrying member; leaving clears; keyboard focus on an anchor triggers the same track.
- [ ] Commit: `feat(series): hover/focus anchor handles drive migration tracking`.

### Task E-4: Derived "phases present" reading

**Files:** Create `src/lib/series/seriesReading.ts` + a readout component. From per-member snapshots compute, per phase: the variable span it covers + its lattice trend (e.g. `a 205 → 195 Å`), plus `coexistence at …` and `form factor only at …` lines. Generalizes to any series — no hand-written "X → Y" narration (the user's concern). For coexistence rows, colour each phase name by `phaseColor` so they self-decode.

```typescript
export interface SeriesReading {
  phases: Array<{ phase: string; spanLabel: string; latticeTrend: string }>;
  coexistenceAt: string[];     // variable values with >1 phase
  formFactorOnlyAt: string[];  // variable values that are form-factor
}
export function seriesReading(members: SeriesMember[], variableOf: (m: SeriesMember) => number | string): SeriesReading
```

- [ ] TDD: a synthetic series (Pn3m a 205→195 across 4 members, one coexistence member, one form-factor member) yields the right span/trend/coexistence/form-factor lines; a no-phase member doesn't appear as an indexed phase.
- [ ] Commit: `feat(series): derived phases-present reading`.

### Task E-5: `SeriesPhaseStrip` companion

**Files:** Create `src/components/SeriesPhaseStrip.tsx` composing `ui/PhaseStrip` (which takes `PhaseSegment[]` with `phase`/`coexistWith`). One cell per sample along the variable; coexistence → two-phase gradient (`coexistWith`); form-factor/null → hollow dashed cell (extend PhaseStrip or render a sibling cell). Derive segments from member snapshots + assignment state.

- [ ] TDD: N samples → N cells in variable order; coexistence cell shows the gradient; form-factor cell renders hollow/dashed; null cell distinct from form-factor.
- [ ] Commit: `feat(series): phase-strip companion`.

### Task E-6: Series member rows (lattice a/d + q₁)

**Files:** Create `src/components/SeriesMemberRow.tsx` (or extend `MemberMetaRow` review-mode). Per member: lattice (`a`/`d`, both shown under coexistence), first-peak q₁ (`min(effective_peaks.q)`), phase name(s) `phaseColor`-coloured so coexistence rows self-decode. Form-factor/null members: "no Bragg peaks · q₁ —".

- [ ] TDD: indexed member shows `a` + q₁; coexistence shows both lattices; form-factor shows the no-Bragg line; phase names are phase-coloured.
- [ ] Commit: `feat(series): member rows with lattice + first-peak q1`.

### Task E-7: Form-factor / null members end-to-end

**Files:** `SeriesBuilderPage.tsx`, `MemberTraceLayer`/`MemberHeatmapLayer` (render a neutral broad-shouldered trace, no anchors), `SeriesPhaseStrip` (hollow cell), `SeriesMemberRow` (no-Bragg line). A member whose assignment state is `form_factor` shows a real trace with no peak anchors; `null` shows a flat/featureless trace. Distinguish the two (decision #2). Depends on Plan A B2 representing the states in the snapshot.

- [ ] TDD: a form-factor member renders a trace + hollow strip cell + no-Bragg row + no anchors; a null member is distinct (flat); neither produces migration vertices.
- [ ] Commit: `feat(series): form-factor and null members end-to-end`.

### Task E-8: Decoupled clean export preset

**Files:** `presets.ts` (`CLEAN_SCIENTIFIC`: white ground, 2px traces, Arial `q (Å⁻¹)`/`Intensity` axes, title, footnote), `multiTraceAdapter.ts` (clean-preset path; include form-factor members; route peaks through `peakMark` with literal `#hex` colours — the export renderer can't resolve CSS vars). One-click Copy PNG / Save SVG, no tuning UI; same samples/order/offset as the live figure.

- [ ] TDD: the adapter emits an `ExportSpec` with white bg + Arial + 2px traces; form-factor members are included; peaks use `peakMark` with resolved hex; the spec round-trips through `buildExportSvg`.
- [ ] Commit: `feat(series): decoupled clean scientific export preset`.

### Task E-9: Wire SeriesBuilderPage + bounded fan-out + full green

**Files:** `SeriesBuilderPage.tsx` — mount the phase-strip, reading, and member rows; **keep the bounded fan-out** (`useMemberTraces`/`useMemberIndices`/`useMemberExposures` via `useStableQueryMap`); do NOT snapshot-embed indices. Preserve the `useStableQueryMap` signature-hash identity (wheel/brush smoothness depends on the per-id data-ref nonce). Guard `ERR_INSUFFICIENT_RESOURCES` if the contact sheet also batches.

- [ ] TDD: builder mounts all new components; the fan-out stays bounded (assert `useQueries` single-tree, not N observers); MultiTracePlot wheel/brush unaffected.
- [ ] `npm test` + `npm run build` + `npm run e2e` (series/compare specs) green.
- [ ] Visual acceptance vs the Series mockup (waterfall, anchor hover track, absent-order ghost, phase-strip, reading, member rows, form-factor handling, clean export).
- [ ] Commit: `feat(series): wire builder page; bounded fan-out; series acceptance`.

---

## Self-Review

**1. Spec coverage** (S1–S7, design §8.3):
- S1 waterfall/heatmap thinking-view → E-1. ✓
- S2 `(phase,order)` migration + absent-order ghost + anchor handles → E-2, E-3. ✓
- S3 phase-strip companion → E-5. ✓
- S4 derived phases-present reading → E-4. ✓
- S5 member rows (lattice + q₁) → E-6. ✓
- S6 form-factor/null members → E-7 (depends on Plan A B2). ✓
- S7 decoupled clean export → E-8. ✓
- Perf cliff bounded, no snapshot-embed → E-9. ✓

**2. Placeholder scan:** the derivations (`seriesReading`, `buildAnchorMap`) have typed signatures + concrete acceptance against the mockup's specific cases (1:0.25 Pn3m √3-absent; a 205→195 trend); component tasks are task-level with the captured props (`PhaseStrip`/`MemberMetaRow`/`multiTraceAdapter`) + the mockup as spec. No vague steps.

**3. Type/name consistency:** `seriesReading`/`SeriesReading`, `buildAnchorMap` (`"phase:order"` key), `SeriesPhaseStrip`/`PhaseSegment`, `SeriesMemberRow`, `CLEAN_SCIENTIFIC` consistent across libs, components, and the export adapter. `(phase, order)` keying matches the existing `CrossTraceTrackingLayer` `(phase, k)` convention.

**Risks:** (a) anchor hit-testing must reuse `PlotSurface.hitTest` so q-link tolerance stays consistent with Focus. (b) form-factor vs null distinction depends on Plan A B2 landing the explicit state in the member snapshot — block E-7 on that. (c) clean-export colour must be literal hex (no CSS var) or canvas renders blank — enforced by `peakMark`'s colour-parameterization (Plan C).

---

## Execution Handoff
1. **Subagent-Driven (recommended)** — fresh subagent per task; E-2 (tracking) and E-8 (export) get careful review.
2. **Inline Execution** — batch with checkpoints after E-4 and E-8.
