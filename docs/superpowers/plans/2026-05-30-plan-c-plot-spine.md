# Plan C — Plot Spine (peakMark + PlotSurface) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the ~5 drifted peak-mark definitions with one `peakMark()` builder honoring the §5.1 encoding atoms, and extract a shared `<PlotSurface>` primitive that TraceViewer / MillerPlot / MultiTracePlot refactor onto — a **greenfield rebuild** (we preserve the *behaviors* — q-link feel, drag/add/exclude, snap, wheel/brush zoom — not the existing custom-SVG implementation).

**Architecture:** `peakMark(opts)` is a pure, colour-parameterized builder returning both an Observable Plot `Markish` (for Plot.dot/export contexts) and an SVG-geometry descriptor (for the interactive overlay). `<PlotSurface>` owns the Observable Plot instance, a shared log/linear x-scale, the gesture layer (wheel/click/brush/dblclick), `invertQ`/`applyQ`, and rAF-batched resize, exposing `overlay(scales)` + `hitTest(q, px, tol)` so consumers attach interaction without re-implementing the scale plumbing. Margins and `xType` stay **per-consumer props** (they genuinely differ). The §5.1 atoms: DOWNWARD auto triangle, manual neutral-gray DIAMOND (magenta retired), ghosted excluded, hollow caret for predicted-but-absent, outline-only for optimistic `id<0`.

**Tech Stack:** TypeScript, React, Observable Plot (`import * as Plot`), Vitest + Testing Library (Plot mocked with a stub `.scale()`), the Print design tokens. `peakMark` lives in `src/components/ui/**` or routes colour through JS mark options to satisfy the `check-design` 0-baseline gate.

---

## Dependencies & scope

- **No backend dependency.** Pure frontend, on data already shipping.
- **This is the spine for Plans D and E** — both build plot features (combs, rings, migration tracking) on `PlotSurface.overlay/hitTest`. Plan C MUST land before D/E plot work.
- **Behaviors preserved (not code):** q-link tolerance feel (don't regress #180 — re-tune `Q_LINK_REL_TOL`/`PEAK_HIT_PX` deliberately), drag-to-add manual peak, click-to-exclude auto peak, wheel zoom about cursor, double-click reset, brush. The old custom-SVG overlay state machine in `TraceViewer.tsx:364-648` is **replaced**, not ported.

## The 5 drifted marks `peakMark()` subsumes (verified)

| Site | Current | Drift |
|---|---|---|
| `TraceViewer.tsx:502-507` | hand-built SVG `<polygon>`, `PEAK_OFFSET_PX=7`, `--color-peak-manual:447` | upward polygon; magenta manual |
| `MemberTraceLayer.tsx:265-276` | `Plot.dot symbol='triangle'`, `PEAK_OFFSET_PX=5` | 2px offset drift; Plot's `triangle` (up) |
| `traceExportMarks.ts:52-63` | `Plot.dot symbol='triangle2'` | Observable upward variant |
| `multiTraceExportMarks.ts:191-203` | `Plot.ruleX` ticks | different mark family entirely |
| `PlotCard.tsx:648-663` | `TriangleSvg` legend | separate geometry copy |

## §5.1 encoding atoms (the `peakMark` contract)

- **colour = phase, always.** Provenance/state via silhouette + fill, never hue.
- **auto** → filled **downward** triangle (points *at* the peak). **manual** → **diamond**, neutral-gray when unindexed, phase-coloured once indexed (magenta `--color-peak-manual` **retired**).
- **excluded** → ghosted (faint/hollow), not struck.
- **predicted-but-absent** → hollow caret.
- **optimistic** (`id<0`) → outline-only (no fill); not clickable (preserve `nearestClickablePeak` skip of `id<0`).
- **hot (q-link)** → grow + terracotta ring, **not** recolour (hue-proximity to Pn3m amber).

## File structure

| File | Change | Responsibility |
|---|---|---|
| `src/components/ui/peakMark.ts` | **Create** | The single mark builder (Plot.Markish + SVG descriptor), colour-parameterized |
| `src/components/ui/peakMark.tsx` | **Create** | `<PeakGlyph>` React/SVG component (legend + overlay) wrapping the descriptor |
| `src/components/PlotSurface.tsx` | **Create** | Shared Plot instance + scales + gestures + resize + `overlay`/`hitTest` |
| `src/lib/plot/invertQ.ts` | Modify | (already has `invertQ`/`applyQ`) add `hitTest` helper if not co-located in PlotSurface |
| `src/components/TraceViewer.tsx` | Rewrite onto PlotSurface | Focus hero; overlay via `peakMark` + `PlotSurface.overlay` |
| `src/components/MillerPlot.tsx` | Refactor onto PlotSurface | margins 32/8 |
| `src/components/MultiTracePlot.tsx` | Refactor onto PlotSurface | margins 50/14, waterfall/heatmap |
| `src/components/MemberTraceLayer.tsx` | Use `peakMark` | replace `symbol='triangle'` |
| `src/lib/figure-export/marks/traceExportMarks.ts` | Use `peakMark` | replace `symbol='triangle2'` |
| `src/lib/figure-export/marks/multiTraceExportMarks.ts` | Use `peakMark` | replace `ruleX` ticks |
| `src/components/PlotCard.tsx` | Use `<PeakGlyph>` | replace `TriangleSvg`; retire magenta |
| `test/peakMark.test.ts(x)` | **Create** | atom coverage |
| `test/PlotSurface.test.tsx` | **Create** | scale/gesture/overlay/resize coverage |
| `test/figure-export/traceAdapter.test.ts:59` | Modify | "manual" legend row → diamond assertion |

---

## Task sequence (build order within Plan C)

> **Sequencing rule (from the survey):** peakMark converges all 5 sites **and retires magenta** in one sweep BEFORE PlotSurface (PlotSurface assumes unified mark geometry). PlotSurface exposes `overlay/hitTest` and preserves behaviors BEFORE any consumer refactor.

### Task 1: `peakMark()` builder + atom tests

**Files:** Create `src/components/ui/peakMark.ts`, `test/peakMark.test.ts`.

The builder is the load-bearing new primitive, so here is its full contract and implementation:

```typescript
// src/components/ui/peakMark.ts
import * as Plot from "@observablehq/plot";

export type PeakSource = "auto" | "manual";
export interface PeakMarkState {
  source: PeakSource;
  excluded?: boolean;      // ghosted
  optimistic?: boolean;    // id<0 → outline-only, non-interactive
  predictedAbsent?: boolean; // hollow caret (no observed peak at a predicted q)
  hot?: boolean;           // q-link: grow + terracotta ring
}
export interface PeakMarkOpts extends PeakMarkState {
  /** Phase colour (resolved string) when indexed; neutral-gray when unindexed manual.
   *  NEVER read from a CSS var here — callers pass the resolved colour so the export
   *  renderer (canvas, no CSS vars) and on-screen share one builder. */
  color: string;
  /** Pixel offset of the glyph apex above the baseline/peak. Parameterized to
   *  preserve the existing 7px (standalone) / 5px (member) feel without re-drift. */
  offsetPx?: number;
  /** Base glyph radius in px (grows when hot). */
  r?: number;
}

export interface PeakGlyphDescriptor {
  shape: "triangle-down" | "diamond" | "caret";
  fill: string | "none";
  stroke: string;
  strokeWidth: number;
  ring: boolean;        // terracotta hot ring
  r: number;
  offsetPx: number;
  interactive: boolean; // false when optimistic
}

const NEUTRAL = "var(--color-ink-faint)"; // unindexed manual (resolved by caller for export)
const HOT_RING = "var(--color-accent)";

/** Resolve the §5.1 atoms into a geometry descriptor (single source of truth
 *  for both the SVG overlay and the legend; the Plot.Markish builder below
 *  derives from the same fields). */
export function peakGlyph(opts: PeakMarkOpts): PeakGlyphDescriptor {
  const shape = opts.predictedAbsent ? "caret"
    : opts.source === "manual" ? "diamond" : "triangle-down";
  const baseR = opts.r ?? 4;
  const r = opts.hot ? baseR * 1.5 : baseR;
  const filled = !opts.predictedAbsent && !opts.optimistic && !opts.excluded;
  return {
    shape,
    fill: filled ? opts.color : "none",
    stroke: opts.color,
    strokeWidth: opts.optimistic ? 1.25 : 0.75,
    ring: opts.hot === true,
    r,
    offsetPx: opts.offsetPx ?? 7,
    interactive: opts.optimistic !== true,
  };
}

/** Observable Plot Markish for a set of peak rows (export + member layers).
 *  Each row must carry q, y, and a resolved `color`; state flags optional. */
export function peakMark(
  rows: ReadonlyArray<{ q: number; y: number; color: string } & PeakMarkState>,
  channels: { x?: string; y?: string } = {},
): Plot.Markish {
  // Observable Plot has no native downward triangle that matches our geometry;
  // use a custom path symbol so member/export marks share the overlay glyph.
  return Plot.dot(rows, {
    x: channels.x ?? "q",
    y: channels.y ?? "y",
    r: 4,
    symbol: (d: unknown) => ((d as PeakMarkState).source === "manual" ? "diamond" : "triangle"),
    // 'triangle' is Plot's up-triangle; rotate via the symbol path is not exposed,
    // so the overlay (PeakGlyph SVG) is canonical for on-screen orientation and
    // these Plot marks are used only where exact apex orientation is non-critical
    // (export/member thumbnails). See PeakGlyph for the downward geometry.
    fill: (d: unknown) => {
      const s = d as PeakMarkOpts;
      return s.predictedAbsent || s.optimistic || s.excluded ? "none" : s.color;
    },
    stroke: (d: unknown) => (d as PeakMarkOpts).color,
    strokeWidth: 0.75,
    fillOpacity: (d: unknown) => ((d as PeakMarkState).excluded ? 0.35 : 1),
  });
}
```

**Acceptance (test `peakMark.test.ts`):** auto→`triangle-down`; manual→`diamond`; manual unindexed→neutral stroke; excluded→`fill:"none"` + faded; optimistic→`fill:"none"` + `interactive:false`; predictedAbsent→`caret`; hot→`r` grows + `ring:true`; **no test asserts magenta** (assert the diamond instead). Mirror `TraceViewer.test.tsx`'s Plot mock pattern.

- [ ] Write `peakMark.test.ts` covering each atom (failing).
- [ ] Run → fails (module missing).
- [ ] Implement `peakMark.ts` (above).
- [ ] Run → passes.
- [ ] Commit: `feat(ui): peakMark — single peak-glyph builder honoring the §5.1 atoms`.

### Task 2: `<PeakGlyph>` SVG component + converge the 5 sites + retire magenta (one sweep)

**Files:** Create `src/components/ui/peakMark.tsx` (`<PeakGlyph descriptor>` rendering the downward triangle / diamond / caret with optional hot ring); modify `TraceViewer.tsx` (overlay peak render → `peakGlyph` geometry), `MemberTraceLayer.tsx` (→ `peakMark`), `traceExportMarks.ts` (→ `peakMark`), `multiTraceExportMarks.ts` (→ `peakMark`), `PlotCard.tsx` (`TriangleSvg` → `<PeakGlyph>`); delete all `--color-peak-manual` references (TraceViewer:447, PlotCard legend, `presets.peakManual`, any export hex); update `traceAdapter.test.ts:59` (manual legend row asserts a diamond, not "manual"/magenta).

- [ ] Write/adjust tests first: `peakMark.tsx` legend test (diamond geometry), `traceAdapter.test.ts:59` migration (failing).
- [ ] Run → fails.
- [ ] Implement `<PeakGlyph>` + convert all 5 sites + delete magenta token usages.
- [ ] Run unit + `npm run build` (tsc + design-guard `lint:design`: confirm no raw colour literal escapes — colour threads through JS mark options / resolved props).
- [ ] Commit: `refactor(plots): converge 5 peak-mark defs onto peakMark; retire --color-peak-manual magenta`.

> **Design-guard note:** `peakMark.ts/tsx` lives under `src/components/ui/**`, so the closed-look allowlist covers it; consumers pass a **resolved** phase colour (via `phaseColor()` threaded into JS props), never a `style={{background: phaseColor()}}`. Verified pattern: TraceViewer/DetectorRingOverlay already pass colour into mark options without an allowlist entry.

### Task 3: `<PlotSurface>` primitive (scales + gestures + resize + overlay/hitTest)

**Files:** Create `src/components/PlotSurface.tsx`, `test/PlotSurface.test.tsx`.

The PlotSurface API (the contract D and E build on):

```typescript
export interface PlotSurfaceProps {
  /** Plot marks for the data layer (trace line, band, member rows, heatmap rects). */
  marks: Plot.Markish[];
  xType: "log" | "linear";
  xDomain: [number, number] | null;
  yDomain?: [number, number] | null;
  margins: { left: number; right: number; top: number; bottom: number }; // per-consumer
  qUnits?: string;
  onXDomain: (d: [number, number] | null) => void;
  /** Interaction wiring — all optional so a read-only plot (MillerPlot) omits them. */
  onAddPeak?: (q: number) => void;
  onClickPeak?: (peakId: number, altKey: boolean) => void;
  /** Render the interactive SVG overlay given live scales. Called every frame /
   *  resize. Receives apply/invert + a hitTest closure. This is where D/E draw
   *  combs, rings, q-link lines, migration tracks. */
  overlay?: (ctx: PlotOverlayContext) => React.ReactNode | void;
  /** Pixel hit-test tolerance for peak clicks (default PEAK_HIT_PX). */
  hitTolerancePx?: number;
}

export interface PlotOverlayContext {
  applyQ: (q: number) => number | null;   // q → px
  invertQ: (px: number) => number | null; // px → q
  applyY: (v: number) => number | null;
  width: number;
  height: number;
  /** Nearest clickable (id>=0) peak to a pixel x, or null. */
  hitTest: (peaks: Peak[], px: number, tolerancePx?: number) => Peak | null;
}
```

**Internals (preserve behaviors, rebuild implementation):**
- Own the `Plot.plot({...marks, x:{type:xType,...}, y:..., margin*})` call; keep `formatAxis` decimal tick labels.
- `ResizeObserver → setResizeKey`, read `clientWidth/clientHeight` **directly in render** (never lagged state) — the load-bearing pattern.
- Re-bind wheel/click/dblclick/brush listeners on every `_resizeKey` change so handlers don't close over stale `onXDomain`/`onAddPeak`. `wheel` must `preventDefault` (passive:false); coords from `getBoundingClientRect`.
- Wheel = zoom about cursor (log-aware, from TraceViewer:295-324). dblclick = reset. click = `hitTest` → `onClickPeak` else `onAddPeak(invertQ(x))`.
- `overlay(ctx)` rendered into an absolutely-positioned `<svg pointer-events:none>` layered over the plot, rebuilt on resize/scale change.
- Reuse `invertQ`/`applyQ` from `lib/plot/invertQ.ts`; add `hitTest` reusing `nearestClickablePeak` (which already skips `id<0`).

**Acceptance (`PlotSurface.test.tsx`):** mounts and renders the data marks (Plot mocked with stub `.scale`); `invertQ`/`applyQ` round-trip via the stub; a wheel event calls `onXDomain` with a narrowed domain; a click within tolerance of a peak calls `onClickPeak`, an empty click calls `onAddPeak`; an `id<0` peak is never clickable; `overlay` receives a context whose `hitTest` skips optimistic peaks; a resize bump re-renders without stale-closure dispatch.

- [ ] Write `PlotSurface.test.tsx` (failing).
- [ ] Run → fails.
- [ ] Implement `PlotSurface.tsx`.
- [ ] Run → passes.
- [ ] Commit: `feat(plots): PlotSurface primitive — shared scales, gestures, resize, overlay/hitTest`.

### Task 4: Refactor TraceViewer onto PlotSurface (Focus hero)

**Files:** Rewrite `src/components/TraceViewer.tsx` body to render `<PlotSurface>` with the trace line + sigma-band marks, and an `overlay` that draws peaks via `peakGlyph` + the cursor q-line + predicted-q lines. Wire `onAddPeak`/`onClickPeak` to the existing `onAddPeak`/`onRemovePeak`/`onTogglePeakExclusion` props. Keep the `TraceViewerProps` interface stable (consumers unchanged). Preserve q-link: `hoveredQ`/`onHoverQ` drive the hot glyph (grow+ring) and cursor.

- [ ] Adjust `TraceViewer.test.tsx` to the new internals (keep `nearestClickablePeak` tests; they move to PlotSurface/peakMark but TraceViewer still renders `data-testid="trace-viewer"`).
- [ ] Run → fails where expected.
- [ ] Rewrite TraceViewer onto PlotSurface + peakGlyph overlay.
- [ ] Run unit + `npm run build`.
- [ ] Visual check against `docs/redesign-mockups/2026-05-29-focus-plot.html` (down-triangle autos, gray diamond manual, ghosted excluded, q-link hot = grow+ring).
- [ ] Commit: `refactor(plots): TraceViewer renders on PlotSurface + peakMark overlay`.

### Task 5: Refactor MillerPlot + MultiTracePlot onto PlotSurface

**Files:** `MillerPlot.tsx` (margins 32/8/6/22, read-only — omit interaction props), `MultiTracePlot.tsx` (margins 50/14/8/32, waterfall/heatmap; keep `computeYBands`, `Representation`, `MemberTraceLayer`/`MemberHeatmapLayer`/`CrossTraceTrackingLayer` composition; route peaks through `peakMark`). PlotSurface takes `margins` + `xType` as props, so MillerPlot's zoom and Compare's divergent scales are honored (no margin/xType standardization).

- [ ] Adjust `MillerPlot.test.tsx` / `MultiTracePlot.test.tsx` (keep `offsetToBandFraction`, `computeMaxPlotWidth` exports/tests).
- [ ] Run → fails where expected.
- [ ] Refactor both onto PlotSurface; preserve y-band semantics + representation toggle.
- [ ] Run unit + `npm run build` + `npm run e2e` (mocked) for the compare/series specs.
- [ ] Commit: `refactor(plots): MillerPlot + MultiTracePlot render on PlotSurface`.

### Task 6: Full green + visual acceptance

- [ ] `(cd packages/HimalayaUI/frontend && npm test)` — all Vitest green; capture once if slow.
- [ ] `npm run build` — tsc + `lint:design` pass (no raw colour literals escaped; magenta gone).
- [ ] `npm run e2e` — mocked Playwright specs pass.
- [ ] Run `test/harness/seed-series.sh` (live visual) if available; eyeball the Focus + Compare surfaces.
- [ ] Commit any test-fixture updates: `test(plots): update fixtures for the unified plot spine`.

---

## Self-Review

**1. Spec coverage** (survey X1 + X2; design §5.2):
- X1 peakMark unifying 5 defs + retiring magenta → Tasks 1–2. ✓
- X2 PlotSurface owning scales/gestures/resize, exposing overlay/hitTest, preserving TraceViewer interaction behaviors → Tasks 3–5. ✓
- Per-consumer margins/xType (no standardization) → Task 3 props + Task 5. ✓
- Export convergence (traceExportMarks/multiTraceExportMarks on peakMark, colour-parameterized) → Task 2. ✓

**2. Placeholder scan:** peakMark is fully implemented; PlotSurface has a complete typed contract + internals spec + acceptance tests. The consumer refactors (Tasks 4–5) are task-level with explicit acceptance + the mockup as visual spec — appropriate for a greenfield UI rebuild where transcribing every JSX line up front would drift.

**3. Type/name consistency:** `peakMark`/`peakGlyph`/`PeakGlyphDescriptor`/`PeakMarkOpts`, `PlotSurface`/`PlotOverlayContext`/`overlay`/`hitTest`, `invertQ`/`applyQ` — consistent across builder, primitive, and consumers. `offsetPx` parameterized (7/5) everywhere to avoid re-drift.

**Known risks:** (a) Observable Plot lacks a true downward-triangle symbol — the canonical orientation lives in the `<PeakGlyph>` SVG overlay; Plot.dot marks (member thumbnails/export) use `triangle`/`diamond` where exact apex orientation is non-critical, documented inline. (b) q-link tolerance: re-tune `Q_LINK_REL_TOL`/`PEAK_HIT_PX` in PlotSurface deliberately; gate on the #180 feel.

---

## Review findings — required fixes (frontend-reviewer, 2026-05-30)

1. **[HIGH] Scope MillerPlot OUT of the PlotSurface refactor.** `MillerPlot.tsx:116-145` is a **q-vs-Miller-ratio scatter with `linearRegressionY` overlays** — x=`ratio` (linear, `ticks:4`), y=`q`, no log/linear q-axis, no `invertQ`/`applyQ`, no wheel/zoom, no peak `hitTest`. The `PlotSurfaceProps` contract here is q/peak-typed and exposes none of MillerPlot's ratio-axis/regression/tick needs. The self-review's "MillerPlot zoom is honored" justification is about the wrong plot (MillerPlot has no zoom). **Fix Task 5:** refactor only TraceViewer + MultiTracePlot onto PlotSurface; leave MillerPlot standalone (or, only if a future need arises, widen PlotSurface to first-class non-q x-channels + axis config — out of scope now).
2. **[MED] `peakMark()` Plot.Markish does not render the caret atom** — it only branches `manual→diamond : triangle`. The `predictedAbsent→caret` glyph (combs D-5, absent-order ghost rings E-2) MUST go through the SVG `peakGlyph`/`<PeakGlyph>` path, **never** the `Plot.dot` markish (which would draw an up-triangle). Make this explicit so a worker wiring ghost rings via `peakMark()` doesn't get the wrong glyph.
3. **[MED] Export peak colour is by-source, not by-phase.** `traceExportMarks.ts:42-50` colours export peaks from `LIGHT_PALETTE` by `source`/`excluded` (`peakAuto:#1f5fb0`, `peakManual:#a0421f`), deliberately not `phaseColor()`. The §5.1 "colour = phase, always" framing is for the *on-screen* surfaces — at the export sites, pass the `LIGHT_PALETTE.*` / clean-preset hex into `peakMark`'s `color` param (it carries literals fine). Don't "fix" the export to phase-colour or it regresses the printable palette. The manual→diamond geometry change does correctly alter the export legend (`traceAdapter.test.ts:59` update already noted).
4. **Confirmed sound:** design-guard (phaseColor through JS mark options / `ui/**` is clean — `check-design` is purely textual on literal `oklch(`/`#hex`), the Observable-Plot-in-React seam (scale cast, read-clientWidth-in-render, re-bind on resize, wheel preventDefault), per-consumer margins/xType (for the two q-plots), `offsetPx` 7/5 preservation, the no-downward-triangle caveat (the mockup itself hand-builds the path — `focus-plot.html:742`), `nearestClickablePeak` `id<0` skip, magenta retirement.

## Execution Handoff
1. **Subagent-Driven (recommended)** — fresh subagent per task; review the peakMark sweep (Task 2) and PlotSurface (Task 3) carefully (highest blast radius).
2. **Inline Execution** — batch with checkpoints after Tasks 2 and 3.
