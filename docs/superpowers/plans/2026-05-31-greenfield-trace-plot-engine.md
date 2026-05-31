# Greenfield Trace-Plot Engine — Implementation Plan (Plan #1: engine + hero/mini proof)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build an owned, declarative, d3-backed single-trace plot engine under `src/print/plot/` — the projection, interaction math, marks, axes, and the `<TracePlot>` component — and prove it in Storybook against the captured real fixtures at two feature levels (the Focus hero: axes + wheel-zoom + click-to-add/select peaks; and the mini: axes and interaction off).

**Architecture:** We own the q→pixel projection via `d3-scale` (no Observable Plot, no DOM scale-readback). Once we own the projection, the entire plot — trace line, ±σ band, peak glyphs, axes, cursor overlay — is one declarative SVG tree driven by `(samples, peaks, projection)`. The trusted interaction *behaviours* (peak hit-test, log-aware zoom-about-cursor) are ported verbatim from `PlotSurface.tsx` as pure functions. `PlotFrame` owns the SVG root + `ResizeObserver` width + non-passive wheel listener; `TracePlot` owns the projection and translates pixel gestures back to q.

**Tech Stack:** React 18 + TypeScript strict (`exactOptionalPropertyTypes` — use conditional-spread for optional props, never explicit `undefined`); `d3-scale` / `d3-shape` (promoted to direct deps); Vitest + React Testing Library; Storybook (CSF3, `@storybook/react-vite`). Tailwind v4 closed-look design system — `src/print/plot/**` is exempted from the inline-appearance guard because it authors SVG appearance.

**Scope note — this is Plan #1 of two.** Plan #1 delivers the engine working and proven on the **single-trace** layout (`traces.length === 1`; >1 overlays in shared scales, which the engine already supports but which no consumer needs yet). **Deferred to Plan #2:** the band/waterfall stacked layout, migrating real consumers (`TraceViewer`/Focus hero, `MultiTracePlot` + member layers), the figure-export migration (four Plot-importing files: `types.ts`, `renderer.ts`, `marks/traceExportMarks.ts`, `marks/multiTraceExportMarks.ts`), deleting the orphaned `MillerPlot` + its test, retiring `PlotSurface`/`peakMark.ts` ×2, and dropping `@observablehq/plot` from `package.json`.

**Verified facts this plan relies on** (checked against live source 2026-05-31):
- `d3-scale@4.0.2`, `d3-shape@3.2.0`, `d3-array@3.2.4`, `d3-format@3.1.2` resolve in `node_modules` **transitively only** (via `@observablehq/plot@0.6.17`). Task 1 promotes the two we use to direct deps + adds their `@types`.
- `hitTestPeaks` lives at `src/components/PlotSurface.tsx:104–123` (skips `id < 0`, skips non-finite px, `<=` so a later equidistant peak wins). `PEAK_HIT_PX = 10` (`PlotSurface.tsx:26`). Wheel zoom math: `PlotSurface.tsx:258–284` (log branch clamps cursor q to `≥ 1e-6`; min-span guard `1e-4`; `factor = Math.exp(deltaY * 0.001)`). dblclick reset: `:290–292` (`onReset?.() ?? onXDomain(null)`).
- `PeakGlyph` (`src/print/ui/PeakGlyph.tsx`) emits a `<g data-role="peak-glyph">` (NOT a standalone `<svg>`); props `{ descriptor, x, y, haloStroke?, opacity?, dataPeakId? }`; color is pre-resolved in the descriptor. The descriptor is built by `peakGlyph(opts)` in `src/print/ui/peakMark.ts` (`{ source, color, predictedAbsent?, excluded?, hot? }`; shape decided there: predicted-absent→caret, manual→diamond, else triangle-down).
- `phaseColor(phase: string): string` (`src/phases.ts:45`) returns an **OKLCH string** (fallback `"oklch(0.65 0.02 270)"`).
- `formatAxis(d: number): string` (`src/lib/plot/formatAxis.ts`) — pure SAXS-axis number formatter, no d3/Plot dep.
- `Trace = { q: number[]; I: number[]; sigma: number[] }` (`src/api.ts:186–190`). `MemberSnapshotPeak = { id; q; intensity: number|null; sharpness; source: "auto"|"manual" }` (`:431–437`).
- Fixtures: `realTraces: Record<number, Trace>` (ids 37/65/66/67/93, ~956 pts each) and `realMembers` / `transitionSeries` in `src/print/fixtures/`.
- Design guard: `scripts/check-design.mjs` — the directory exemption is `isExcluded(relPath)` at **line 134–137**; Storybook glob is `../src/print/**/*.stories.@(ts|tsx)` (`.storybook/main.ts`), so `src/print/plot/**/*.stories.tsx` auto-discovers. Vitest tests live under `test/print-*` (e.g. `test/print-ui/`), not colocated — this plan uses `test/print-plot/`.
- `npm` scripts: `test` = `vitest run`; `build` = `npm run lint:design && tsc --noEmit -p tsconfig.build.json && vite build`; `lint:design` = `node scripts/check-design.mjs`; `storybook` = `storybook dev -p 6006 --no-open`.

**Commit convention:** every commit message in this plan must end with the repo's required trailer:
```
Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
```
(Omitted from the example commands below for brevity — add it on each commit.)

---

## File Structure

All paths relative to `packages/HimalayaUI/frontend/`.

| File | Responsibility |
|---|---|
| `src/print/plot/projection.ts` (create) | Pure d3-scale wrapper: `makeAxis`, `makeProjection`, `positiveExtent`. No DOM. |
| `src/print/plot/interaction.ts` (create) | Ported pure fns + constants: `hitTestPeaks`, `zoomXDomain`, `PEAK_HIT_PX`, `BRUSH_DEADZONE_PX`, `MIN_SPAN_FRAC`. |
| `src/print/plot/PlotFrame.tsx` (create) | SVG root + `ResizeObserver` width + non-passive wheel + click/dblclick; projection-agnostic (emits pixel events). |
| `src/print/plot/Axis.tsx` (create) | One axis (bottom/left) as JSX over `axis.ticks()`; mono tick labels via `formatAxis`. |
| `src/print/plot/marks/TraceLine.tsx` (create) | Trace line + ±σ band `<path>`s via `d3-shape`. |
| `src/print/plot/marks/PlotPeaks.tsx` (create) | Peak marks — composes `print/ui/PeakGlyph` inside the plot SVG. Owns `PlotPeak` type. |
| `src/print/plot/TracePlot.tsx` (create) | The engine: owns projection + gesture→q translation; assembles Axis + marks + overlay. |
| `src/print/plot/index.ts` (create) | Barrel. |
| `src/print/plot/TracePlot.stories.tsx` (create) | Hero + mini stories against fixtures. |
| `scripts/check-design.mjs` (modify, line 134–137) | Add `print/plot/` to `isExcluded`. |
| `package.json` (modify) | Promote `d3-scale`/`d3-shape` to direct deps + `@types/*` devDeps. |
| `test/print-plot/*.test.ts(x)` (create) | Unit + RTL tests for each unit above. |

---

## Task 1: Foundation — d3 deps + design-guard exemption

**Files:**
- Modify: `package.json`
- Modify: `scripts/check-design.mjs:134-137`

- [ ] **Step 1: Promote d3 modules to direct dependencies**

The engine imports `d3-scale` and `d3-shape` directly. They currently resolve only transitively via `@observablehq/plot`; depending on a hoisted transitive is fragile (and breaks when Plot is removed in Plan #2). Install the two runtime modules + the four `@types` packages (d3 ships no bundled types):

Run (from `packages/HimalayaUI/frontend/`):
```bash
npm install d3-scale@^4.0.2 d3-shape@^3.2.0
npm install -D @types/d3-scale @types/d3-shape @types/d3-array @types/d3-format
```
Expected: `package.json` gains `"d3-scale"` and `"d3-shape"` under `dependencies` and the four `@types/d3-*` under `devDependencies`; `npm` exits 0.

- [ ] **Step 2: Exempt `src/print/plot/**` from the inline-appearance guard**

The plot engine authors SVG appearance (stroke widths, phase colours, axis fonts), exactly like `print/ui/`. Add it to the directory exclusion. Open `scripts/check-design.mjs` and change `isExcluded` (currently lines 134–137):

Before:
```js
// src/components/ui/** and src/print/ui/** are excluded entirely (where appearance is authored).
function isExcluded(relPath) {
  return relPath.startsWith("components/ui/") || relPath.startsWith("print/ui/");
}
```
After:
```js
// src/components/ui/**, src/print/ui/**, and src/print/plot/** are excluded
// entirely (where appearance is authored — primitives and the plot engine).
function isExcluded(relPath) {
  return (
    relPath.startsWith("components/ui/") ||
    relPath.startsWith("print/ui/") ||
    relPath.startsWith("print/plot/")
  );
}
```

- [ ] **Step 3: Verify the guard still runs clean**

Run: `npm run lint:design`
Expected: exits 0 (no `src/print/plot/` files exist yet; the edit is forward-compatible).

- [ ] **Step 4: Commit**

```bash
git add package.json package-lock.json scripts/check-design.mjs
git commit -m "build(plot): promote d3-scale/d3-shape to direct deps + exempt print/plot from design guard"
```

---

## Task 2: Projection — d3-scale wrapper

**Files:**
- Create: `src/print/plot/projection.ts`
- Test: `test/print-plot/projection.test.ts`

- [ ] **Step 1: Write the failing test**

Create `test/print-plot/projection.test.ts`:
```ts
import { describe, it, expect } from "vitest";
import {
  makeAxis,
  makeProjection,
  positiveExtent,
} from "../../src/print/plot/projection";

describe("makeAxis", () => {
  it("linear to/invert round-trips", () => {
    const a = makeAxis([0, 10], [0, 100], "linear");
    expect(a.to(5)).toBeCloseTo(50);
    expect(a.invert(50)).toBeCloseTo(5);
  });

  it("log maps decades evenly and inverts", () => {
    const a = makeAxis([0.01, 1], [0, 200], "log");
    expect(a.to(0.01)).toBeCloseTo(0);
    expect(a.to(1)).toBeCloseTo(200);
    expect(a.to(0.1)).toBeCloseTo(100);
    expect(a.invert(100)).toBeCloseTo(0.1);
  });

  it("ticks returns values inside the domain", () => {
    const a = makeAxis([0.01, 1], [0, 200], "log");
    const t = a.ticks(5);
    expect(t.length).toBeGreaterThan(0);
    for (const v of t) {
      expect(v).toBeGreaterThanOrEqual(0.01);
      expect(v).toBeLessThanOrEqual(1);
    }
  });
});

describe("makeProjection", () => {
  it("inverts the y range so the domain max sits at the top (px 0)", () => {
    const p = makeProjection({
      xDomain: [0.01, 1],
      yDomain: [1, 1000],
      plotWidth: 300,
      plotHeight: 200,
      xType: "log",
      yType: "log",
    });
    expect(p.y.to(1000)).toBeCloseTo(0);
    expect(p.y.to(1)).toBeCloseTo(200);
    expect(p.x.to(0.01)).toBeCloseTo(0);
    expect(p.x.to(1)).toBeCloseTo(300);
  });
});

describe("positiveExtent", () => {
  it("ignores non-finite and non-positive values", () => {
    expect(positiveExtent([0, -1, 2, NaN, 8])).toEqual([2, 8]);
  });

  it("falls back when there is no positive data", () => {
    expect(positiveExtent([0, -5, NaN])).toEqual([1, 10]);
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npx vitest run test/print-plot/projection.test.ts`
Expected: FAIL — cannot resolve `../../src/print/plot/projection`.

- [ ] **Step 3: Write the implementation**

Create `src/print/plot/projection.ts`:
```ts
import { scaleLinear, scaleLog } from "d3-scale";

export type ScaleType = "log" | "linear";

/** One axis: data↔pixel mapping plus "nice" tick values. Backed by d3-scale. */
export interface Axis1D {
  /** data value → pixel. */
  to(value: number): number;
  /** pixel → data value. */
  invert(px: number): number;
  /** "Nice" tick values across the domain (default ~6). */
  ticks(count?: number): number[];
  readonly domain: readonly [number, number];
  readonly range: readonly [number, number];
  readonly type: ScaleType;
}

export function makeAxis(
  domain: [number, number],
  range: [number, number],
  type: ScaleType,
): Axis1D {
  const scale =
    type === "log" ? scaleLog<number, number>() : scaleLinear<number, number>();
  scale.domain(domain);
  scale.range(range);
  return {
    to: (v: number) => scale(v),
    invert: (px: number) => scale.invert(px),
    ticks: (count = 6) => scale.ticks(count),
    domain,
    range,
    type,
  };
}

export interface Projection {
  x: Axis1D;
  y: Axis1D;
}

/**
 * Build x (left→right) and y (top→bottom *inverted*: domain max at px 0) axes
 * for a plot area of `plotWidth × plotHeight` pixels.
 */
export function makeProjection(opts: {
  xDomain: [number, number];
  yDomain: [number, number];
  plotWidth: number;
  plotHeight: number;
  xType?: ScaleType;
  yType?: ScaleType;
}): Projection {
  const {
    xDomain,
    yDomain,
    plotWidth,
    plotHeight,
    xType = "log",
    yType = "log",
  } = opts;
  return {
    x: makeAxis(xDomain, [0, plotWidth], xType),
    y: makeAxis(yDomain, [plotHeight, 0], yType),
  };
}

/**
 * Tight [min, max] over the positive, finite values — the domain a log scale
 * needs. Falls back to [1, 10] when there is no positive data, and pads a
 * degenerate single-value extent so the scale is non-empty.
 */
export function positiveExtent(
  values: number[],
  fallback: [number, number] = [1, 10],
): [number, number] {
  let min = Infinity;
  let max = -Infinity;
  for (const v of values) {
    if (!Number.isFinite(v) || v <= 0) continue;
    if (v < min) min = v;
    if (v > max) max = v;
  }
  if (min === Infinity) return fallback;
  if (min === max) return [min * 0.9, max * 1.1];
  return [min, max];
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npx vitest run test/print-plot/projection.test.ts`
Expected: PASS (4 tests / 6 assertions green).

- [ ] **Step 5: Commit**

```bash
git add src/print/plot/projection.ts test/print-plot/projection.test.ts
git commit -m "feat(plot): owned d3-scale projection (makeAxis/makeProjection/positiveExtent)"
```

---

## Task 3: Interaction — ported pure behaviours

**Files:**
- Create: `src/print/plot/interaction.ts`
- Test: `test/print-plot/interaction.test.ts`

- [ ] **Step 1: Write the failing test**

Create `test/print-plot/interaction.test.ts`:
```ts
import { describe, it, expect } from "vitest";
import {
  hitTestPeaks,
  zoomXDomain,
  PEAK_HIT_PX,
} from "../../src/print/plot/interaction";

const toPx = (q: number) => q * 100; // q=1 → 100px

describe("hitTestPeaks", () => {
  it("returns the nearest peak within tolerance", () => {
    const peaks = [
      { id: 1, q: 1 },
      { id: 2, q: 2 },
    ];
    expect(hitTestPeaks(peaks, 105, toPx, PEAK_HIT_PX)?.id).toBe(1);
  });

  it("returns null beyond tolerance", () => {
    const peaks = [{ id: 1, q: 1 }];
    expect(hitTestPeaks(peaks, 130, toPx, PEAK_HIT_PX)).toBeNull();
  });

  it("skips optimistic placeholders (id < 0)", () => {
    const peaks = [{ id: -1, q: 1 }];
    expect(hitTestPeaks(peaks, 100, toPx, PEAK_HIT_PX)).toBeNull();
  });

  it("the later peak wins an exact tie (<=)", () => {
    const peaks = [
      { id: 1, q: 1 },
      { id: 2, q: 1 },
    ];
    expect(hitTestPeaks(peaks, 100, toPx, PEAK_HIT_PX)?.id).toBe(2);
  });
});

describe("zoomXDomain — log", () => {
  it("zooms in about the cursor, narrowing the domain", () => {
    const next = zoomXDomain({
      cursorQ: 0.1,
      deltaY: -100,
      current: [0.01, 1],
      extent: [0.01, 1],
      type: "log",
    });
    expect(next).not.toBeNull();
    expect(next![0]).toBeGreaterThanOrEqual(0.01);
    expect(next![1]).toBeLessThanOrEqual(1);
    expect(next![1] - next![0]).toBeLessThan(1 - 0.01);
  });

  it("clamps to the full extent when zooming out", () => {
    const next = zoomXDomain({
      cursorQ: 0.1,
      deltaY: 1000,
      current: [0.05, 0.2],
      extent: [0.01, 1],
      type: "log",
    });
    expect(next![0]).toBeGreaterThanOrEqual(0.01);
    expect(next![1]).toBeLessThanOrEqual(1);
  });

  it("returns null when the span would get too small", () => {
    const tiny = zoomXDomain({
      cursorQ: 0.1,
      deltaY: -100000,
      current: [0.0999, 0.1001],
      extent: [0.01, 1],
      type: "log",
    });
    expect(tiny).toBeNull();
  });
});

describe("zoomXDomain — linear", () => {
  it("narrows about the cursor in linear space", () => {
    const next = zoomXDomain({
      cursorQ: 5,
      deltaY: -100,
      current: [0, 10],
      extent: [0, 10],
      type: "linear",
    });
    expect(next![1] - next![0]).toBeLessThan(10);
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npx vitest run test/print-plot/interaction.test.ts`
Expected: FAIL — cannot resolve `../../src/print/plot/interaction`.

- [ ] **Step 3: Write the implementation**

Create `src/print/plot/interaction.ts`. The bodies of `hitTestPeaks` and `zoomXDomain` are ported verbatim from `PlotSurface.tsx` (hit-test `:104-123`, wheel math `:258-284`) — do not "improve" the math; behaviour parity is the point.
```ts
/**
 * Pure interaction math ported from the legacy PlotSurface so the greenfield
 * engine reproduces its feel exactly. No DOM, no d3 — caller supplies the q→px
 * projection.
 */

/** Click within this many pixels of a peak to act on it. (PlotSurface #180.) */
export const PEAK_HIT_PX = 10;
/** Drags shorter than this are treated as clicks, not brushes. */
export const BRUSH_DEADZONE_PX = 4;
/** Refuse zooms narrower than this fraction of the full extent. */
export const MIN_SPAN_FRAC = 1e-4;

export interface HitPeak {
  id: number;
  q: number;
}

/**
 * Nearest clickable peak (id ≥ 0) to a pixel x, within `tolerancePx`, or null.
 * On an exact tie the later peak wins (`<=`), matching PlotSurface:104-123.
 */
export function hitTestPeaks<T extends HitPeak>(
  peaks: T[],
  clickX: number,
  toPx: (q: number) => number,
  tolerancePx: number,
): T | null {
  let best: T | null = null;
  let bestDist = tolerancePx;
  for (const p of peaks) {
    if (p.id < 0) continue;
    const px = toPx(p.q);
    if (!Number.isFinite(px)) continue;
    const d = Math.abs(px - clickX);
    if (d <= bestDist) {
      best = p;
      bestDist = d;
    }
  }
  return best;
}

/**
 * Wheel-zoom the x domain about the cursor (log-aware). Returns the new
 * [min, max] clamped to `extent`, or null when the zoom would collapse the
 * span below MIN_SPAN_FRAC of the extent. Ported from PlotSurface:258-284.
 */
export function zoomXDomain(opts: {
  cursorQ: number;
  deltaY: number;
  current: [number, number];
  extent: [number, number];
  type: "log" | "linear";
}): [number, number] | null {
  const { cursorQ, deltaY, current, extent, type } = opts;
  const [curMin, curMax] = current;
  const [q0, qN] = extent;
  const factor = Math.exp(deltaY * 0.001);
  let newMin: number;
  let newMax: number;
  if (type === "log") {
    const logMin = Math.log(curMin);
    const logMax = Math.log(curMax);
    const logCur = Math.log(Math.max(cursorQ, 1e-6));
    newMin = Math.max(q0, Math.exp(logCur - (logCur - logMin) * factor));
    newMax = Math.min(qN, Math.exp(logCur + (logMax - logCur) * factor));
  } else {
    newMin = Math.max(q0, cursorQ - (cursorQ - curMin) * factor);
    newMax = Math.min(qN, cursorQ + (curMax - cursorQ) * factor);
  }
  if (newMax - newMin < (qN - q0) * MIN_SPAN_FRAC) return null;
  return [newMin, newMax];
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npx vitest run test/print-plot/interaction.test.ts`
Expected: PASS (6 tests green).

- [ ] **Step 5: Commit**

```bash
git add src/print/plot/interaction.ts test/print-plot/interaction.test.ts
git commit -m "feat(plot): port hitTestPeaks + log-aware zoomXDomain as pure fns"
```

---

## Task 4: PlotFrame — SVG root + ResizeObserver + gestures

**Files:**
- Create: `src/print/plot/PlotFrame.tsx`
- Test: `test/print-plot/PlotFrame.test.tsx`

**Design note:** `PlotFrame` is projection-agnostic. It measures its own width (skipped when an explicit `width` is supplied — used by tests and fixed-size stories), computes the plot area from `margins`, and hands a `PlotDims` to a `render` callback. It emits gestures as **container-relative pixels** (`onWheelPx`/`onClickPx`/`onDoubleClickPx`); `TracePlot` translates them to q. The wheel listener is attached via `addEventListener(..., { passive: false })` in an effect because React's synthetic `onWheel` is passive and cannot `preventDefault()` (this is why the legacy code bound wheel imperatively).

- [ ] **Step 1: Write the failing test**

Create `test/print-plot/PlotFrame.test.tsx`:
```tsx
import { render } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { PlotFrame } from "../../src/print/plot/PlotFrame";

describe("PlotFrame", () => {
  const margins = { top: 10, right: 10, bottom: 20, left: 30 };

  it("renders an svg and calls render with the computed plot dims", () => {
    const renderSpy = vi.fn(() => null);
    render(
      <PlotFrame
        height={200}
        margins={margins}
        width={400}
        data-testid="frame"
        render={renderSpy}
      />,
    );
    expect(renderSpy).toHaveBeenCalledWith(
      expect.objectContaining({
        width: 400,
        height: 200,
        plotWidth: 360, // 400 - 30 - 10
        plotHeight: 170, // 200 - 10 - 20
        margins,
      }),
    );
  });

  it("renders a testid'd svg element", () => {
    const { container } = render(
      <PlotFrame
        height={200}
        margins={margins}
        width={400}
        data-testid="frame"
        render={() => null}
      />,
    );
    expect(container.querySelector('svg[data-testid="frame"]')).toBeTruthy();
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npx vitest run test/print-plot/PlotFrame.test.tsx`
Expected: FAIL — cannot resolve `../../src/print/plot/PlotFrame`.

- [ ] **Step 3: Write the implementation**

Create `src/print/plot/PlotFrame.tsx`:
```tsx
import { useEffect, useRef, useState, type ReactNode } from "react";

export interface Margins {
  top: number;
  right: number;
  bottom: number;
  left: number;
}

export interface PlotDims {
  width: number;
  height: number;
  plotWidth: number;
  plotHeight: number;
  margins: Margins;
}

export interface PlotFrameProps {
  height: number;
  margins: Margins;
  /** Controlled width; when omitted the frame measures its container. */
  width?: number;
  /** Width used before the first measurement (and in non-DOM tests). */
  defaultWidth?: number;
  className?: string;
  "data-testid"?: string;
  /** Container-relative pixel gestures. q-translation happens in the caller. */
  onWheelPx?: (deltaY: number, px: number, py: number) => void;
  onClickPx?: (px: number, py: number, altKey: boolean) => void;
  onDoubleClickPx?: () => void;
  /** Render the plot body, given the resolved pixel dimensions. */
  render: (dims: PlotDims) => ReactNode;
}

export function PlotFrame({
  height,
  margins,
  width,
  defaultWidth = 640,
  className,
  "data-testid": testid,
  onWheelPx,
  onClickPx,
  onDoubleClickPx,
  render,
}: PlotFrameProps): JSX.Element {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const [measured, setMeasured] = useState<number | null>(null);

  // Measure container width unless an explicit width is supplied.
  useEffect(() => {
    if (width !== undefined) return;
    const el = containerRef.current;
    if (!el || typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver((entries) => {
      const w = entries[0]?.contentRect.width ?? 0;
      if (w > 0) setMeasured(w);
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, [width]);

  // React's onWheel is passive and cannot preventDefault(); bind non-passive.
  useEffect(() => {
    const el = containerRef.current;
    if (!el || !onWheelPx) return;
    function handle(ev: WheelEvent): void {
      ev.preventDefault();
      const rect = el!.getBoundingClientRect();
      onWheelPx!(ev.deltaY, ev.clientX - rect.left, ev.clientY - rect.top);
    }
    el.addEventListener("wheel", handle, { passive: false });
    return () => el.removeEventListener("wheel", handle);
  }, [onWheelPx]);

  const w = width ?? measured ?? defaultWidth;
  const plotWidth = Math.max(0, w - margins.left - margins.right);
  const plotHeight = Math.max(0, height - margins.top - margins.bottom);
  const dims: PlotDims = { width: w, height, plotWidth, plotHeight, margins };

  function handleClick(ev: React.MouseEvent): void {
    if (!onClickPx) return;
    const rect = containerRef.current!.getBoundingClientRect();
    onClickPx(ev.clientX - rect.left, ev.clientY - rect.top, ev.altKey);
  }

  return (
    <div
      ref={containerRef}
      className={className}
      style={{ width: width !== undefined ? width : "100%" }}
    >
      <svg
        width={w}
        height={height}
        viewBox={`0 0 ${w} ${height}`}
        data-testid={testid}
        onClick={handleClick}
        {...(onDoubleClickPx ? { onDoubleClick: onDoubleClickPx } : {})}
      >
        <g transform={`translate(${margins.left},${margins.top})`}>
          {render(dims)}
        </g>
      </svg>
    </div>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npx vitest run test/print-plot/PlotFrame.test.tsx`
Expected: PASS (2 tests green).

- [ ] **Step 5: Commit**

```bash
git add src/print/plot/PlotFrame.tsx test/print-plot/PlotFrame.test.tsx
git commit -m "feat(plot): PlotFrame SVG root with ResizeObserver width + non-passive wheel"
```

---

## Task 5: Axis — ticks as JSX

**Files:**
- Create: `src/print/plot/Axis.tsx`
- Test: `test/print-plot/Axis.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/print-plot/Axis.test.tsx`:
```tsx
import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeAxis } from "../../src/print/plot/projection";
import { Axis } from "../../src/print/plot/Axis";

describe("Axis", () => {
  it("renders one tick label per tick value (bottom)", () => {
    const a = makeAxis([0.01, 1], [0, 300], "log");
    const ticks = a.ticks(6);
    const { container } = render(
      <svg>
        <Axis
          axis={a}
          orientation="bottom"
          plotWidth={300}
          plotHeight={150}
        />
      </svg>,
    );
    const texts = container.querySelectorAll('[data-role="axis-bottom"] text');
    expect(texts.length).toBe(ticks.length);
  });

  it("renders a label when provided (left)", () => {
    const a = makeAxis([1, 1000], [150, 0], "log");
    const { container } = render(
      <svg>
        <Axis
          axis={a}
          orientation="left"
          plotWidth={300}
          plotHeight={150}
          label="I"
        />
      </svg>,
    );
    const labels = [...container.querySelectorAll("text")].filter(
      (t) => t.textContent === "I",
    );
    expect(labels.length).toBe(1);
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npx vitest run test/print-plot/Axis.test.tsx`
Expected: FAIL — cannot resolve `../../src/print/plot/Axis`.

- [ ] **Step 3: Write the implementation**

Create `src/print/plot/Axis.tsx`. Mono tick labels honour the "Mono-Means-Measurement" DESIGN.md invariant; colours/fonts use CSS tokens (the file is design-guard-exempt, so inline `style` with var tokens is allowed and is the natural fit for SVG):
```tsx
import { type Axis1D } from "./projection";
import { formatAxis } from "../../lib/plot/formatAxis";

export interface AxisProps {
  axis: Axis1D;
  orientation: "bottom" | "left";
  plotWidth: number;
  plotHeight: number;
  label?: string;
}

const TICK_LABEL_STYLE: React.CSSProperties = {
  fontFamily: "var(--font-mono)",
  fontSize: 10,
  fill: "var(--color-ink-soft)",
};

const AXIS_LABEL_STYLE: React.CSSProperties = {
  fontFamily: "var(--font-sans)",
  fontSize: 11,
  fill: "var(--color-ink-soft)",
};

export function Axis({
  axis,
  orientation,
  plotWidth,
  plotHeight,
  label,
}: AxisProps): JSX.Element {
  const ticks = axis.ticks(6);
  const isBottom = orientation === "bottom";
  return (
    <g data-role={`axis-${orientation}`} aria-hidden="true">
      {ticks.map((t) => {
        const p = axis.to(t);
        if (!Number.isFinite(p)) return null;
        return isBottom ? (
          <g key={t} transform={`translate(${p},${plotHeight})`}>
            <line y2={5} stroke="var(--color-ink-faint)" strokeWidth={1} />
            <text y={16} textAnchor="middle" style={TICK_LABEL_STYLE}>
              {formatAxis(t)}
            </text>
          </g>
        ) : (
          <g key={t} transform={`translate(0,${p})`}>
            <line x2={-5} stroke="var(--color-ink-faint)" strokeWidth={1} />
            <text x={-8} dy="0.32em" textAnchor="end" style={TICK_LABEL_STYLE}>
              {formatAxis(t)}
            </text>
          </g>
        );
      })}
      {label != null &&
        (isBottom ? (
          <text
            x={plotWidth / 2}
            y={36}
            textAnchor="middle"
            style={AXIS_LABEL_STYLE}
          >
            {label}
          </text>
        ) : (
          <text
            transform={`translate(-40,${plotHeight / 2}) rotate(-90)`}
            textAnchor="middle"
            style={AXIS_LABEL_STYLE}
          >
            {label}
          </text>
        ))}
    </g>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npx vitest run test/print-plot/Axis.test.tsx`
Expected: PASS (2 tests green).

- [ ] **Step 5: Commit**

```bash
git add src/print/plot/Axis.tsx test/print-plot/Axis.test.tsx
git commit -m "feat(plot): Axis — d3 ticks rendered as mono-labelled JSX"
```

---

## Task 6: TraceLine mark — line + σ-band

**Files:**
- Create: `src/print/plot/marks/TraceLine.tsx`
- Test: `test/print-plot/TraceLine.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/print-plot/TraceLine.test.tsx`:
```tsx
import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeProjection } from "../../src/print/plot/projection";
import { TraceLine } from "../../src/print/plot/marks/TraceLine";

describe("TraceLine", () => {
  const proj = makeProjection({
    xDomain: [0.1, 0.3],
    yDomain: [1, 30],
    plotWidth: 300,
    plotHeight: 200,
    xType: "log",
    yType: "log",
  });

  it("draws a band path and a line path (line starts with M)", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 2, 1.5] };
    const { container } = render(
      <svg>
        <TraceLine trace={trace} projection={proj} color="oklch(0.5 0.1 30)" />
      </svg>,
    );
    const paths = container.querySelectorAll('[data-role="trace-line"] path');
    expect(paths.length).toBe(2); // band + line
    expect(paths[1]!.getAttribute("d")).toMatch(/^M/);
  });

  it("omits the band when band={false}", () => {
    const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 2, 1.5] };
    const { container } = render(
      <svg>
        <TraceLine
          trace={trace}
          projection={proj}
          color="oklch(0.5 0.1 30)"
          band={false}
        />
      </svg>,
    );
    const paths = container.querySelectorAll('[data-role="trace-line"] path');
    expect(paths.length).toBe(1);
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npx vitest run test/print-plot/TraceLine.test.tsx`
Expected: FAIL — cannot resolve `../../src/print/plot/marks/TraceLine`.

- [ ] **Step 3: Write the implementation**

Create `src/print/plot/marks/TraceLine.tsx`:
```tsx
import { line as d3line, area as d3area } from "d3-shape";
import { type Projection } from "../projection";
import { type Trace } from "../../../api";

export interface TraceLineProps {
  trace: Trace;
  projection: Projection;
  /** Resolved stroke/fill colour (e.g. phaseColor output). */
  color: string;
  /** Draw the ±σ band behind the line (default true; needs trace.sigma). */
  band?: boolean;
}

export function TraceLine({
  trace,
  projection,
  color,
  band = true,
}: TraceLineProps): JSX.Element {
  const { x, y } = projection;
  const n = Math.min(trace.q.length, trace.I.length);
  const idx = Array.from({ length: n }, (_, i) => i);
  const valid = (i: number): boolean =>
    Number.isFinite(trace.q[i]!) &&
    trace.q[i]! > 0 &&
    Number.isFinite(trace.I[i]!) &&
    trace.I[i]! > 0;

  const linePath =
    d3line<number>()
      .defined(valid)
      .x((i) => x.to(trace.q[i]!))
      .y((i) => y.to(trace.I[i]!))(idx) ?? "";

  const sigma = trace.sigma;
  const bandPath =
    band && sigma
      ? d3area<number>()
          .defined((i) => valid(i) && Number.isFinite(sigma[i]!))
          .x((i) => x.to(trace.q[i]!))
          .y0((i) => y.to(Math.max(trace.I[i]! - sigma[i]!, 1e-9)))
          .y1((i) => y.to(trace.I[i]! + sigma[i]!))(idx) ?? ""
      : "";

  return (
    <g data-role="trace-line">
      {bandPath ? (
        <path d={bandPath} fill={color} opacity={0.15} stroke="none" />
      ) : null}
      <path
        d={linePath}
        fill="none"
        stroke={color}
        strokeWidth={1.25}
        strokeLinejoin="round"
      />
    </g>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npx vitest run test/print-plot/TraceLine.test.tsx`
Expected: PASS (2 tests green).

- [ ] **Step 5: Commit**

```bash
git add src/print/plot/marks/TraceLine.tsx test/print-plot/TraceLine.test.tsx
git commit -m "feat(plot): TraceLine mark — d3-shape line + sigma band"
```

---

## Task 7: PlotPeaks mark — composes PeakGlyph

**Files:**
- Create: `src/print/plot/marks/PlotPeaks.tsx`
- Test: `test/print-plot/PlotPeaks.test.tsx`

**Composition contract (verified):** `PeakGlyph` emits a `<g data-role="peak-glyph">` and must live inside a parent `<svg>`; it takes a `descriptor` built by `peakGlyph(opts)` and a pre-resolved colour. `PlotPeaks` maps each peak's `(q, intensity)` to plot pixels via the projection and renders one `<PeakGlyph>` per real peak. Peaks with `id < 0` (optimistic placeholders) are skipped; peaks with null intensity anchor to `baselineI`.

- [ ] **Step 1: Write the failing test**

Create `test/print-plot/PlotPeaks.test.tsx`:
```tsx
import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { makeProjection } from "../../src/print/plot/projection";
import { PlotPeaks, type PlotPeak } from "../../src/print/plot/marks/PlotPeaks";

describe("PlotPeaks", () => {
  const proj = makeProjection({
    xDomain: [0.1, 0.3],
    yDomain: [1, 30],
    plotWidth: 300,
    plotHeight: 200,
    xType: "log",
    yType: "log",
  });

  it("renders one glyph per real peak and skips id < 0", () => {
    const peaks: PlotPeak[] = [
      { id: 1, q: 0.15, source: "auto", intensity: 10 },
      { id: -2, q: 0.2, source: "auto", intensity: 5 },
    ];
    const { container } = render(
      <svg>
        <PlotPeaks peaks={peaks} projection={proj} color="oklch(0.5 0.1 30)" />
      </svg>,
    );
    expect(
      container.querySelectorAll('[data-role="peak-glyph"]').length,
    ).toBe(1);
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npx vitest run test/print-plot/PlotPeaks.test.tsx`
Expected: FAIL — cannot resolve `../../src/print/plot/marks/PlotPeaks`.

- [ ] **Step 3: Write the implementation**

Create `src/print/plot/marks/PlotPeaks.tsx`:
```tsx
import { type Projection } from "../projection";
import { PeakGlyph } from "../../ui/PeakGlyph";
import { peakGlyph } from "../../ui/peakMark";

export interface PlotPeak {
  id: number;
  q: number;
  intensity?: number | null;
  source: "auto" | "manual";
  excluded?: boolean;
  predictedAbsent?: boolean;
  hot?: boolean;
}

export interface PlotPeaksProps {
  peaks: PlotPeak[];
  projection: Projection;
  /** Resolved colour for the glyphs (phase colour, usually). */
  color: string;
  /** y data-value for peaks lacking an intensity (anchors them to baseline). */
  baselineI?: number;
  /** Paper colour threaded to PeakGlyph's halo (export-parity). */
  paperColor?: string;
}

export function PlotPeaks({
  peaks,
  projection,
  color,
  baselineI,
  paperColor,
}: PlotPeaksProps): JSX.Element {
  const { x, y } = projection;
  return (
    <g data-role="plot-peaks">
      {peaks.map((p) => {
        if (p.id < 0) return null;
        const px = x.to(p.q);
        if (!Number.isFinite(px)) return null;
        const iVal = p.intensity ?? baselineI;
        const py =
          iVal != null && Number.isFinite(iVal) ? y.to(iVal) : y.range[0];
        const descriptor = peakGlyph({
          source: p.source,
          color,
          ...(p.predictedAbsent ? { predictedAbsent: true } : {}),
          ...(p.excluded ? { excluded: true } : {}),
          ...(p.hot ? { hot: true } : {}),
        });
        return (
          <PeakGlyph
            key={p.id}
            descriptor={descriptor}
            x={px}
            y={py}
            dataPeakId={p.id}
            {...(paperColor ? { haloStroke: paperColor } : {})}
          />
        );
      })}
    </g>
  );
}
```

> Note: `y.range[0]` is the bottom-of-plot pixel (the inverted y range starts at `plotHeight`), so null-intensity peaks sit on the baseline.

- [ ] **Step 4: Run the test to verify it passes**

Run: `npx vitest run test/print-plot/PlotPeaks.test.tsx`
Expected: PASS (1 test green).

- [ ] **Step 5: Verify the `peakGlyph` opts shape against live source**

Before relying on the opts, confirm the signature (the descriptor builder owns shape selection):
Run: `grep -n "export function peakGlyph" src/print/ui/peakMark.ts`
Expected: a signature accepting `{ source, color, predictedAbsent?, excluded?, hot? }`. If a field name differs, fix the call in `PlotPeaks.tsx` and re-run Step 4.

- [ ] **Step 6: Commit**

```bash
git add src/print/plot/marks/PlotPeaks.tsx test/print-plot/PlotPeaks.test.tsx
git commit -m "feat(plot): PlotPeaks mark — composes print/ui PeakGlyph in plot space"
```

---

## Task 8: TracePlot — the engine + barrel

**Files:**
- Create: `src/print/plot/TracePlot.tsx`
- Create: `src/print/plot/index.ts`
- Test: `test/print-plot/TracePlot.test.tsx`

**Design note:** `TracePlot` owns the projection and the gesture→q translation. Because the projection depends on the measured plot width (owned by `PlotFrame`), it is built inside `PlotFrame`'s `render` callback and stashed in a ref so the stable gesture handlers read the latest projection without stale closures. Click translation subtracts `margins.left` (gestures arrive container-relative; the projection's x range is plot-area-relative).

- [ ] **Step 1: Write the failing test**

Create `test/print-plot/TracePlot.test.tsx`:
```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { TracePlot, type TraceModel } from "../../src/print/plot/TracePlot";

const model: TraceModel = {
  trace: { q: [0.1, 0.2, 0.3], I: [10, 20, 15], sigma: [1, 1, 1] },
  peaks: [{ id: 1, q: 0.2, source: "auto", intensity: 20 }],
  phase: "Ia3d",
};

describe("TracePlot", () => {
  it("renders axes, a trace line, and peak glyphs for a single trace", () => {
    const { container } = render(
      <TracePlot traces={[model]} width={500} height={300} data-testid="tp" />,
    );
    expect(container.querySelector('svg[data-testid="tp"]')).toBeTruthy();
    expect(
      container.querySelector('[data-role="trace-line"] path'),
    ).toBeTruthy();
    expect(
      container.querySelectorAll('[data-role="peak-glyph"]').length,
    ).toBe(1);
    expect(container.querySelector('[data-role="axis-bottom"]')).toBeTruthy();
    expect(container.querySelector('[data-role="axis-left"]')).toBeTruthy();
  });

  it("renders no axes when axes={false} (mini level)", () => {
    const { container } = render(
      <TracePlot traces={[model]} width={120} height={48} axes={false} />,
    );
    expect(container.querySelector('[data-role="axis-bottom"]')).toBeNull();
    expect(
      container.querySelector('[data-role="trace-line"] path'),
    ).toBeTruthy();
  });

  it("calls onAddPeak when clicking empty space", () => {
    const onAddPeak = vi.fn();
    const { container } = render(
      <TracePlot
        traces={[model]}
        width={500}
        height={300}
        interaction={{ onXDomain: () => {}, onAddPeak }}
      />,
    );
    // Click far from the single peak (q=0.2 ≈ px 274 in plot space). px=148.
    fireEvent.click(container.querySelector("svg")!, {
      clientX: 200,
      clientY: 100,
    });
    expect(onAddPeak).toHaveBeenCalledTimes(1);
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npx vitest run test/print-plot/TracePlot.test.tsx`
Expected: FAIL — cannot resolve `../../src/print/plot/TracePlot`.

- [ ] **Step 3: Write the implementation**

Create `src/print/plot/TracePlot.tsx`:
```tsx
import { useCallback, useRef } from "react";
import { PlotFrame, type Margins, type PlotDims } from "./PlotFrame";
import {
  makeProjection,
  positiveExtent,
  type Projection,
  type ScaleType,
} from "./projection";
import { Axis } from "./Axis";
import { TraceLine } from "./marks/TraceLine";
import { PlotPeaks, type PlotPeak } from "./marks/PlotPeaks";
import { hitTestPeaks, zoomXDomain, PEAK_HIT_PX } from "./interaction";
import { phaseColor } from "../../phases";
import { type Trace } from "../../api";

export interface TraceModel {
  trace: Trace;
  peaks: PlotPeak[];
  phase: string | null;
}

export interface PlotContext {
  projection: Projection;
  dims: PlotDims;
  hitTest: (peaks: PlotPeak[], px: number, tol?: number) => PlotPeak | null;
}

export interface TracePlotInteraction {
  onXDomain: (d: [number, number] | null) => void;
  onAddPeak?: (q: number) => void;
  onClickPeak?: (peakId: number, altKey: boolean) => void;
  onReset?: () => void;
  hitTolerancePx?: number;
}

export interface TracePlotProps {
  /** 1 = hero/mini; >1 overlays in shared scales. */
  traces: TraceModel[];
  /** Visible x window; null = full data extent. */
  xDomain?: [number, number] | null;
  xType?: ScaleType;
  yType?: ScaleType;
  axes?: boolean;
  xLabel?: string;
  yLabel?: string;
  interaction?: TracePlotInteraction | false;
  overlay?: (ctx: PlotContext) => React.ReactNode;
  height?: number;
  width?: number;
  /** Paper colour for peak-glyph halos. */
  paperColor?: string;
  className?: string;
  "data-testid"?: string;
}

const UNINDEXED_COLOR = "var(--color-ink-faint)";

interface PlotState {
  projection: Projection;
  dims: PlotDims;
  curXDomain: [number, number];
  xExtent: [number, number];
  xType: ScaleType;
  interaction: TracePlotInteraction | null;
  peaks: PlotPeak[];
}

export function TracePlot(props: TracePlotProps): JSX.Element {
  const {
    traces,
    xDomain = null,
    xType = "log",
    yType = "log",
    axes = true,
    xLabel = "q (Å⁻¹)",
    yLabel = "I",
    interaction = false,
    overlay,
    height = 320,
    width,
    paperColor,
    className,
    "data-testid": testid,
  } = props;

  const xExtent = positiveExtent(traces.flatMap((t) => t.trace.q));
  const yExtent = positiveExtent(traces.flatMap((t) => t.trace.I));
  const curXDomain = xDomain ?? xExtent;

  const margins: Margins = axes
    ? { top: 12, right: 14, bottom: 40, left: 52 }
    : { top: 4, right: 4, bottom: 4, left: 4 };

  const stateRef = useRef<PlotState | null>(null);

  const handleWheelPx = useCallback((deltaY: number, px: number) => {
    const s = stateRef.current;
    if (!s || !s.interaction) return;
    const cursorQ = s.projection.x.invert(px - s.dims.margins.left);
    if (!Number.isFinite(cursorQ)) return;
    const next = zoomXDomain({
      cursorQ,
      deltaY,
      current: s.curXDomain,
      extent: s.xExtent,
      type: s.xType,
    });
    if (next) s.interaction.onXDomain(next);
  }, []);

  const handleClickPx = useCallback(
    (px: number, _py: number, altKey: boolean) => {
      const s = stateRef.current;
      if (!s || !s.interaction) return;
      const plotPx = px - s.dims.margins.left;
      const q = s.projection.x.invert(plotPx);
      if (!Number.isFinite(q)) return;
      const tol = s.interaction.hitTolerancePx ?? PEAK_HIT_PX;
      const hit = hitTestPeaks(
        s.peaks,
        plotPx,
        (qq) => s.projection.x.to(qq),
        tol,
      );
      if (hit && s.interaction.onClickPeak) {
        s.interaction.onClickPeak(hit.id, altKey);
      } else if (s.interaction.onAddPeak) {
        s.interaction.onAddPeak(q);
      }
    },
    [],
  );

  const handleDblClick = useCallback(() => {
    const s = stateRef.current;
    if (!s || !s.interaction) return;
    if (s.interaction.onReset) s.interaction.onReset();
    else s.interaction.onXDomain(null);
  }, []);

  return (
    <PlotFrame
      height={height}
      margins={margins}
      {...(width !== undefined ? { width } : {})}
      {...(className ? { className } : {})}
      {...(testid ? { "data-testid": testid } : {})}
      {...(interaction
        ? {
            onWheelPx: handleWheelPx,
            onClickPx: handleClickPx,
            onDoubleClickPx: handleDblClick,
          }
        : {})}
      render={(dims) => {
        const projection = makeProjection({
          xDomain: curXDomain,
          yDomain: yExtent,
          plotWidth: dims.plotWidth,
          plotHeight: dims.plotHeight,
          xType,
          yType,
        });
        const allPeaks = traces.flatMap((t) => t.peaks);
        stateRef.current = {
          projection,
          dims,
          curXDomain,
          xExtent,
          xType,
          interaction: interaction || null,
          peaks: allPeaks,
        };
        const ctx: PlotContext = {
          projection,
          dims,
          hitTest: (ps, px, tol) =>
            hitTestPeaks(ps, px, (q) => projection.x.to(q), tol ?? PEAK_HIT_PX),
        };
        return (
          <>
            {axes ? (
              <>
                <Axis
                  axis={projection.x}
                  orientation="bottom"
                  plotWidth={dims.plotWidth}
                  plotHeight={dims.plotHeight}
                  label={xLabel}
                />
                <Axis
                  axis={projection.y}
                  orientation="left"
                  plotWidth={dims.plotWidth}
                  plotHeight={dims.plotHeight}
                  label={yLabel}
                />
              </>
            ) : null}
            {traces.map((t, i) => (
              <TraceLine
                key={`line-${i}`}
                trace={t.trace}
                projection={projection}
                color={t.phase ? phaseColor(t.phase) : UNINDEXED_COLOR}
              />
            ))}
            {traces.map((t, i) => (
              <PlotPeaks
                key={`peaks-${i}`}
                peaks={t.peaks}
                projection={projection}
                color={t.phase ? phaseColor(t.phase) : UNINDEXED_COLOR}
                baselineI={yExtent[0]}
                {...(paperColor ? { paperColor } : {})}
              />
            ))}
            {overlay ? overlay(ctx) : null}
          </>
        );
      }}
    />
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npx vitest run test/print-plot/TracePlot.test.tsx`
Expected: PASS (3 tests green).

- [ ] **Step 5: Write the barrel**

Create `src/print/plot/index.ts`:
```ts
export { TracePlot } from "./TracePlot";
export type {
  TracePlotProps,
  TraceModel,
  PlotContext,
  TracePlotInteraction,
} from "./TracePlot";
export type { PlotPeak } from "./marks/PlotPeaks";
export { TraceLine } from "./marks/TraceLine";
export { PlotPeaks } from "./marks/PlotPeaks";
export { Axis } from "./Axis";
export { PlotFrame } from "./PlotFrame";
export type { Margins, PlotDims } from "./PlotFrame";
export {
  makeProjection,
  makeAxis,
  positiveExtent,
  type Projection,
  type Axis1D,
  type ScaleType,
} from "./projection";
export {
  hitTestPeaks,
  zoomXDomain,
  PEAK_HIT_PX,
  BRUSH_DEADZONE_PX,
  MIN_SPAN_FRAC,
} from "./interaction";
```

- [ ] **Step 6: Typecheck the whole engine**

Run: `npx tsc --noEmit -p tsconfig.build.json`
Expected: exits 0 (no type errors). If `exactOptionalPropertyTypes` flags an optional-prop assignment, convert it to a conditional spread (`...(x ? { x } : {})`).

- [ ] **Step 7: Commit**

```bash
git add src/print/plot/TracePlot.tsx src/print/plot/index.ts test/print-plot/TracePlot.test.tsx
git commit -m "feat(plot): TracePlot engine — projection + marks + axes + gesture→q"
```

---

## Task 9: Storybook hero + mini proof against real fixtures

**Files:**
- Create: `src/print/plot/TracePlot.stories.tsx`

**Goal of this task:** prove the engine carries the hardest single-trace case (the Focus hero) and the cheapest (mini) against the captured real data, and confirm the design guard + build stay green.

- [ ] **Step 1: Write the stories**

Create `src/print/plot/TracePlot.stories.tsx`. The model adapter maps a `SeriesMember`'s snapshot peaks + measured trace into a `TraceModel`:
```tsx
import type { Meta, StoryObj } from "@storybook/react-vite";
import { TracePlot, type TraceModel } from "./TracePlot";
import type { PlotPeak } from "./marks/PlotPeaks";
import { realTraces } from "../fixtures/realTraces";
import { realMembers, transitionSeries } from "../fixtures/realSeriesMembers";
import type { SeriesMember } from "../../api";

function modelFor(member: SeriesMember): TraceModel {
  const trace = realTraces[member.exposure_id as number]!;
  const peaks: PlotPeak[] = (member.snapshot?.effective_peaks ?? []).map(
    (p) => ({
      id: p.id,
      q: p.q,
      intensity: p.intensity,
      source: p.source,
    }),
  );
  return {
    trace,
    peaks,
    phase: member.snapshot?.confirmed_index?.phase ?? null,
  };
}

const heroModel = modelFor(realMembers[0]!); // exp 65, Ia3d

const meta = {
  title: "plot/TracePlot",
  component: TracePlot,
  parameters: { layout: "padded" },
} satisfies Meta<typeof TracePlot>;

export default meta;
type Story = StoryObj<typeof meta>;

// Focus hero: axes on, full interaction (zoom + click-to-add/select peaks).
export const Hero: Story = {
  args: {
    traces: [heroModel],
    height: 360,
    interaction: {
      onXDomain: (d) => console.log("xDomain", d),
      onAddPeak: (q) => console.log("addPeak", q),
      onClickPeak: (id, alt) => console.log("clickPeak", id, alt),
    },
    paperColor: "var(--color-paper)",
  },
};

// Mini: axes + interaction off, small — cheap enough for a masonry of cards.
export const Mini: Story = {
  args: {
    traces: [heroModel],
    height: 64,
    width: 200,
    axes: false,
    interaction: false,
  },
};

// Overlay: the three-member Sample-9 transition in shared scales (a preview of
// the Plan #2 band layout — here simply overlaid).
export const TransitionOverlay: Story = {
  args: {
    traces: transitionSeries.map(modelFor),
    height: 360,
  },
};
```

- [ ] **Step 2: Verify the fixtures expose what the stories import**

Run:
```bash
grep -n "export const realTraces" src/print/fixtures/realTraces.ts
grep -n "export const realMembers\|export const transitionSeries" src/print/fixtures/realSeriesMembers.ts
```
Expected: `realTraces` (a `Record<number, Trace>`), `realMembers`, and `transitionSeries` all present. If `exposure_id` is nullable on `SeriesMember` and `tsc` complains about the index, the synthetic empty-state members are not used here (only `realMembers`/`transitionSeries`, which have real exposure ids) — the `as number` assertion in `modelFor` is intentional and safe for these fixtures.

- [ ] **Step 3: Typecheck (stories are compiled by the build)**

Run: `npx tsc --noEmit -p tsconfig.build.json`
Expected: exits 0.

- [ ] **Step 4: Run the design guard**

Run: `npm run lint:design`
Expected: exits 0 — `src/print/plot/**` is exempt (Task 1), so its inline SVG appearance does not trip the guard.

- [ ] **Step 5: Build Storybook to confirm the stories compile**

Run: `npm run build-storybook`
Expected: exits 0; the `plot/TracePlot` stories (`Hero`, `Mini`, `TransitionOverlay`) are in the built catalog. (Optional manual check: `npm run storybook` and eyeball the hero — log trace + σ band + Ia3d-coloured peak glyphs + mono axes; scroll-wheel zooms about the cursor; clicking empty space logs `addPeak`, clicking a glyph logs `clickPeak`; double-click resets.)

- [ ] **Step 6: Run the full unit suite + build gate**

Run: `npx vitest run test/print-plot/` then `npm run build`
Expected: all `test/print-plot/` specs green; `npm run build` (lint:design + tsc + vite build) exits 0.

- [ ] **Step 7: Commit**

```bash
git add src/print/plot/TracePlot.stories.tsx
git commit -m "feat(plot): Storybook hero/mini/overlay proof on real fixtures"
```

---

## Self-Review (run after implementing all tasks)

**1. Spec coverage** — map each spec section to a task:
- Own-the-projection (`projection(domain,range,type)`) → Task 2. ✓
- Keep/port behaviours (zoom-about-cursor, hitTestPeaks, dblclick reset) → Task 3 (zoom + hit-test) + Task 8 (dblclick wired). ✓ *(Note: the spec's "q-link relative tolerance" does not exist in source — it is pixel `PEAK_HIT_PX`; reflected here, fix the spec.)*
- Drop the machinery (Plot host, scale-readback, two-SVG stack, resize dance) → not ported; the engine has none. ✓
- Shared spine (phaseColor, PeakGlyph vocabulary, log-x default, mono labels) → Task 5 (mono Axis), Task 7 (PeakGlyph), Task 8 (phaseColor, log defaults). ✓
- Data model `TraceModel` + measured feeder → Task 8 (`TraceModel`), Task 9 (`modelFor`). The modelled `buildMiniWaterfall` feeder is a Plan #2 concern (band layout). ✓
- File layout (`src/print/plot/` + internal primitives + design-guard allowlist) → Tasks 1–8. ✓
- Component surface (`traces`, `axes`, `interaction`, `overlay`, `size`/`className`) → Task 8 (`layout` prop dropped from the sketch as YAGNI for single-trace; overlay/axes/interaction/width/className present). ✓
- Testing (unit projection/hit-test/zoom + Storybook stories) → Tasks 2,3,9. ✓
- Migration steps 3–4 (consumers, export, drop Plot) → **Plan #2** (out of scope here, by design). ✓

**2. Placeholder scan** — no "TBD"/"add error handling"/"similar to Task N"; every code step shows complete code. ✓

**3. Type consistency** — `Projection`/`Axis1D`/`ScaleType` (Task 2) used identically in Tasks 6–8; `PlotPeak` defined in Task 7 and imported by Task 8; `PlotDims`/`Margins` defined in Task 4 and consumed in Task 8; `peakGlyph` opts verified against live source in Task 7 Step 5. ✓

**4. Open questions deferred to Plan #2** (from the spec): line-path decimation threshold for large drag-zooms (the hero trace is ~956 pts — fine undecimated; decimation only matters once many overlaid traces land), and figure-export sequencing.

---

## Execution Handoff

Plan complete and saved to `docs/superpowers/plans/2026-05-31-greenfield-trace-plot-engine.md`. Two execution options:

**1. Subagent-Driven (recommended)** — dispatch a fresh subagent per task, two-stage review (spec compliance, then code quality) between tasks, fast iteration.

**2. Inline Execution** — execute tasks in this session with batched checkpoints for review.

Which approach?
