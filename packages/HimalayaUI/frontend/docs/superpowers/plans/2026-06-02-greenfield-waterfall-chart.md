# WaterfallChart Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the Series-surface stacked-offset `WaterfallChart` renderer by composing N single-trace `TracePlot`s, after a small `TracePlot`/`TraceLine` refactor.

**Architecture:** Each member trace is its own axis-less `TracePlot` box (structural per-row partitioning replaces the legacy `computeYBands`/`applyNormalization` math), stacked bottom-up under one shared bottom `Axis`. Peak anchors come free from the existing peaks layer (`confirmed_index.peak_ids`); row hot/dim and the sample label live in a thin wrapper. Pure renderer — props in, `onHoverRow`/`onHoverQ` out.

**Tech Stack:** React 18 + TS strict (`exactOptionalPropertyTypes` ON → conditional spread for optional/`data-*` props, never explicit `undefined`), Vite, Storybook CSF3 (`@storybook/react-vite`, auto-globs `src/print/**/*.stories.tsx`), Vitest + RTL under JSDOM (assert via `data-role`/`data-*`/text, never class strings — JSDOM computes no layout). d3-shape (via the existing marks). Spec: `docs/superpowers/specs/2026-06-02-greenfield-waterfall-chart-design.md`.

**Working dir:** all commands run from `packages/HimalayaUI/frontend/`. The Vitest pre-tool hook injects `--run`, so `npm test -- <pattern>` is one-shot. Commit trailer on every commit: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

---

## File structure

| File | Responsibility |
|---|---|
| `src/print/plot/marks/TraceLine.tsx` (modify) | Accept a `color` prop for the line stroke (default keeps ink-soft). |
| `src/print/plot/TracePlot.tsx` (modify) | `traces: TraceModel[]` → single `trace: TraceModel`; thread phase colour to the line; add `yHeadroom` top pad. |
| `src/print/plot/TracePlot.stories.tsx` (modify) | Update to the single-trace API; delete the multi-trace `TransitionOverlay` story (superseded by the waterfall). |
| `src/print/waterfall/waterfallModel.ts` (create) | Pure view-model: `WaterfallAnchor`/`WaterfallRow`, `toWaterfallRows`, `waterfallQDomain`. |
| `src/print/waterfall/WaterfallChart.tsx` (create) | The renderer: stacks N single-trace `TracePlot`s + shared axis + labels + row hover. |
| `src/print/waterfall/waterfall.fixtures.ts` (create) | Plain module deriving `WaterfallRow[]` from `realMembers` + `realTraces`. |
| `src/print/waterfall/WaterfallChart.stories.tsx` (create) | Storybook stories (Wide/Narrow/Hover/States). |
| `src/print/waterfall/index.ts` (create) | Barrel. |
| `test/print-plot/TraceLine.test.tsx` (create) | TraceLine colour test. |
| `test/print-waterfall/waterfallModel.test.ts` (create) | Model unit tests. |
| `test/print-waterfall/WaterfallChart.test.tsx` (create) | Renderer unit tests. |

---

## Task 1: TraceLine `color` prop

The waterfall needs phase-coloured lines; today `TraceLine` hard-codes `stroke="var(--color-ink-soft)"`. Add an optional `color` (default preserves current behaviour).

**Files:**
- Modify: `src/print/plot/marks/TraceLine.tsx`
- Test: `test/print-plot/TraceLine.test.tsx` (create)

- [ ] **Step 1: Write the failing test**

Create `test/print-plot/TraceLine.test.tsx`:

```tsx
import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { TraceLine } from "../../src/print/plot/marks/TraceLine";
import { makeProjection } from "../../src/print/plot/projection";
import type { Trace } from "../../src/api";

const trace: Trace = { q: [0.02, 0.05, 0.1], I: [100, 50, 10], sigma: [1, 1, 1] };
const projection = makeProjection({
  xDomain: [0.02, 0.1], yDomain: [10, 100],
  plotWidth: 200, plotHeight: 100, xType: "log", yType: "log",
});

// The line is the second <path> when a band is present; with band={false} it is
// the only path. Use band={false} so the line path is unambiguous.
function linePath(container: HTMLElement): SVGPathElement {
  const paths = container.querySelectorAll('[data-role="trace-line"] path');
  return paths[paths.length - 1] as SVGPathElement;
}

describe("TraceLine color", () => {
  it("defaults the line stroke to ink-soft", () => {
    const { container } = render(
      <svg><TraceLine trace={trace} projection={projection} band={false} /></svg>,
    );
    expect(linePath(container).getAttribute("stroke")).toBe("var(--color-ink-soft)");
  });

  it("uses the provided color for the line stroke", () => {
    const { container } = render(
      <svg><TraceLine trace={trace} projection={projection} band={false} color="oklch(0.570 0.150 58)" /></svg>,
    );
    expect(linePath(container).getAttribute("stroke")).toBe("oklch(0.570 0.150 58)");
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npm test -- print-plot/TraceLine`
Expected: the second test FAILS (stroke is `var(--color-ink-soft)`, not the passed colour).

- [ ] **Step 3: Implement the `color` prop**

In `src/print/plot/marks/TraceLine.tsx`, add `color` to the props interface and default it, then use it on the line `<path>`:

```tsx
export interface TraceLineProps {
  trace: Trace;
  projection: Projection;
  /** Draw the ±σ band behind the line (default true; needs trace.sigma). */
  band?: boolean;
  /** Line stroke colour (default the neutral ink-soft). */
  color?: string;
}

export function TraceLine({
  trace,
  projection,
  band = true,
  color = "var(--color-ink-soft)",
}: TraceLineProps): JSX.Element {
```

Then change the line `<path>` stroke from the literal to `color`:

```tsx
      <path
        d={linePath}
        fill="none"
        stroke={color}
        strokeWidth={1.8}
        strokeLinejoin="round"
      />
```

(Leave the band `<path>` fill as-is.)

- [ ] **Step 4: Run the test to verify it passes**

Run: `npm test -- print-plot/TraceLine`
Expected: both tests PASS.

- [ ] **Step 5: tsc + commit**

Run: `npx tsc --noEmit -p tsconfig.build.json` → exit 0.

```bash
git add src/print/plot/marks/TraceLine.tsx test/print-plot/TraceLine.test.tsx
git commit -m "feat(plot): TraceLine accepts a line color prop (default ink-soft)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: TracePlot → single trace + phase-coloured line + `yHeadroom`

Narrow `TracePlot` to one trace (the only way it is used), colour the line by phase, and add an optional top headroom so peaks don't slam the ceiling when stacked. Update the stories.

**Files:**
- Modify: `src/print/plot/TracePlot.tsx`
- Modify: `src/print/plot/TracePlot.stories.tsx`
- Test: `test/print-plot/TracePlot.test.tsx` (existing — update to new API), plus a new headroom/colour assertion.

- [ ] **Step 1: Inspect the existing TracePlot tests**

Run: `ls test/print-plot/ && grep -rln "traces" test/print-plot/`
This lists the test files referencing the old `traces` array. You will update each `traces={[m]}` / `traces: [m]` to `trace={m}` / `trace: m` in Step 4. Read them now so the refactor keeps their intent.

- [ ] **Step 2: Write the new failing assertions**

Append to `test/print-plot/TracePlot.test.tsx` (adjust the import block if the helpers already exist there):

```tsx
import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { TracePlot } from "../../src/print/plot/TracePlot";
import type { TraceModel } from "../../src/print/plot/TracePlot";

const model: TraceModel = {
  trace: { q: [0.02, 0.05, 0.1], I: [100, 50, 10], sigma: [1, 1, 1] },
  peaks: [{ id: 1, q: 0.05, intensity: 50, source: "auto" }],
  phase: "Ia3d",
};

describe("TracePlot single-trace API", () => {
  it("colours the trace line by phase (not ink-soft)", () => {
    const { container } = render(<TracePlot trace={model} height={200} />);
    const paths = container.querySelectorAll('[data-role="trace-line"] path');
    const line = paths[paths.length - 1] as SVGPathElement;
    // phaseColor("Ia3d") is a non-ink-soft oklch token.
    expect(line.getAttribute("stroke")).not.toBe("var(--color-ink-soft)");
    expect(line.getAttribute("stroke")).toContain("oklch");
  });

  it("falls back to ink-faint for a null-phase trace", () => {
    const { container } = render(
      <TracePlot trace={{ ...model, phase: null }} height={200} />,
    );
    const paths = container.querySelectorAll('[data-role="trace-line"] path');
    const line = paths[paths.length - 1] as SVGPathElement;
    expect(line.getAttribute("stroke")).toBe("var(--color-ink-faint)");
  });
});
```

- [ ] **Step 3: Run to verify it fails**

Run: `npm test -- print-plot/TracePlot`
Expected: FAIL — `trace` is not a prop yet (TS/runtime error or the line is ink-soft).

- [ ] **Step 4: Refactor `TracePlot.tsx`**

Apply these edits to `src/print/plot/TracePlot.tsx`:

(a) Props interface — replace the `traces` field:

```tsx
export interface TracePlotProps {
  /** The single trace to render (hero / mini / one waterfall row). */
  trace: TraceModel;
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
  /** Multiply the y-domain top by (1 + yHeadroom) so peaks keep headroom below
   *  the ceiling — used by the stacked waterfall. Default 0 (no change). */
  yHeadroom?: number;
  show?: { peaks?: boolean; labels?: boolean; band?: boolean };
  highlightPeakIds?: ReadonlySet<number>;
}
```

(b) Destructure — replace `traces` with `trace`, add `yHeadroom`:

```tsx
  const {
    trace,
    xDomain = null,
    xType = "log",
    yType = "log",
    axes = true,
    xLabel = "q (Å⁻¹)",
    yLabel = "intensity (a.u.)",
    interaction = false,
    overlay,
    height = 320,
    width,
    paperColor,
    className,
    "data-testid": testid,
    yHeadroom = 0,
    show,
    highlightPeakIds,
  } = props;
```

(c) Extents — replace the `flatMap` extents and apply headroom:

```tsx
  const xExtent = positiveExtent(trace.trace.q);
  const rawYExtent = positiveExtent(trace.trace.I);
  const yExtent: [number, number] = [rawYExtent[0], rawYExtent[1] * (1 + yHeadroom)];
  const curXDomain = xDomain ?? xExtent;
```

(d) `stateRef` peaks — replace `allPeaks`:

```tsx
        const allPeaks = trace.peaks;
```

(e) Trace line — replace the `traces.map(...)` line block with a single coloured line:

```tsx
            <TraceLine
              trace={trace.trace}
              projection={projection}
              band={layers.band}
              color={trace.phase ? phaseColor(trace.phase) : UNINDEXED_COLOR}
            />
```

(f) Peaks layer — replace the `traces.map(...)` peaks block with a single layer:

```tsx
            {layers.peaks ? (() => {
              const peaksWithHover = hoverId == null
                ? trace.peaks
                : trace.peaks.map((p) => (p.id === hoverId ? { ...p, hot: true } : p));
              return (
                <PlotPeaks
                  peaks={peaksWithHover}
                  projection={projection}
                  color={trace.phase ? phaseColor(trace.phase) : UNINDEXED_COLOR}
                  baselineI={yExtent[0]}
                  {...(paperColor ? { paperColor } : {})}
                  {...(interaction ? { onPeakFocus: setHoverId } : {})}
                  {...(highlightPeakIds ? { highlightPeakIds } : {})}
                />
              );
            })() : null}
```

(g) Labels layer — replace the `traces.map(...)` labels block with a single layer:

```tsx
            {layers.labels ? (
              <PlotLabels
                peaks={trace.peaks}
                projection={projection}
                color={trace.phase ? phaseColor(trace.phase) : UNINDEXED_COLOR}
                baselineI={yExtent[0]}
                {...(highlightPeakIds ? { highlightPeakIds } : {})}
              />
            ) : null}
```

(h) q-readout — replace the `traces.flatMap(...).find(...)` with a single-trace find:

```tsx
              const hovered = hoverId == null || !layers.peaks
                ? null
                : trace.peaks.find((p) => p.id === hoverId) ?? null;
```

- [ ] **Step 5: Update the stories to the single-trace API**

In `src/print/plot/TracePlot.stories.tsx`:

- Delete the `TransitionOverlay` story (lines defining `export const TransitionOverlay`) — it is the only multi-trace usage and is superseded by `WaterfallChart`.
- Remove `transitionSeries` from the `realSeriesMembers` import (line 6) since it is now unused.
- Replace every `traces: [heroModel]` → `trace: heroModel`, `traces: [annotatedModel]` → `trace: annotatedModel`, and the two JSX `traces={[heroModel]}` / `traces={[annotatedModel]}` → `trace={heroModel}` / `trace={annotatedModel}`.

- [ ] **Step 6: Update the existing TracePlot tests to the new API**

In `test/print-plot/TracePlot.test.tsx` (and any sibling print-plot test using TracePlot), change each `traces={[m]}` → `trace={m}` and `traces: [m]` → `trace: m`. Any test that passed a 2+ element `traces` array must be split or dropped (the API is single-trace now; if such a test asserted overlay behaviour, delete it — the waterfall covers stacking).

- [ ] **Step 7: Run the full print-plot suite**

Run: `npm test -- print-plot`
Expected: all PASS (including the two new colour assertions).

- [ ] **Step 8: tsc + commit**

Run: `npx tsc --noEmit -p tsconfig.build.json` → exit 0.

```bash
git add src/print/plot/TracePlot.tsx src/print/plot/TracePlot.stories.tsx test/print-plot/
git commit -m "refactor(plot): TracePlot single-trace + phase-coloured line + yHeadroom

Narrows traces[] -> trace (only ever used single-trace), colours the line
by phase, adds an optional top headroom for stacked use. Drops the
multi-trace TransitionOverlay story (superseded by WaterfallChart).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `waterfallModel.ts` — pure view-model

**Files:**
- Create: `src/print/waterfall/waterfallModel.ts`
- Test: `test/print-waterfall/waterfallModel.test.ts`

- [ ] **Step 1: Write the failing test**

Create `test/print-waterfall/waterfallModel.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { toWaterfallRows, waterfallQDomain } from "../../src/print/waterfall/waterfallModel";
import { realMembers, formFactorMember, unindexedMember } from "../../src/print/fixtures/realSeriesMembers";
import { realTraces } from "../../src/print/fixtures/realTraces";

describe("toWaterfallRows", () => {
  it("returns one row per member, preserving display order", () => {
    const rows = toWaterfallRows(realMembers, realTraces);
    expect(rows.length).toBe(realMembers.length);
    const orders = realMembers.map((m) => m.display_order);
    // rows are in the same order as the (already display-ordered) input
    expect(rows.map((_, i) => realMembers[i]!.display_order)).toEqual(orders);
  });

  it("derives the dominant phase from confirmed_phases / confirmed_index", () => {
    const rows = toWaterfallRows([realMembers[0]!], realTraces);
    expect(rows[0]!.phase).toBe("Ia3d"); // exp 65
  });

  it("emits an anchor per confirmed_index.peak_id, resolved to its effective-peak q", () => {
    const rows = toWaterfallRows([realMembers[0]!], realTraces);
    const ci = realMembers[0]!.snapshot.confirmed_index!;
    expect(rows[0]!.anchors.map((a) => a.id).sort()).toEqual([...ci.peak_ids].sort());
    // every anchor q matches the effective_peak with that id
    const byId = new Map(realMembers[0]!.snapshot.effective_peaks.map((p) => [p.id, p.q]));
    for (const a of rows[0]!.anchors) expect(a.q).toBe(byId.get(a.id));
  });

  it("yields zero anchors and a null phase for a form-factor member", () => {
    const rows = toWaterfallRows([formFactorMember], realTraces);
    expect(rows[0]!.anchors).toEqual([]);
    expect(rows[0]!.phase).toBeNull();
    expect(rows[0]!.state).toBe("form_factor");
  });

  it("yields zero anchors for an unindexed (null-assignment) member", () => {
    const rows = toWaterfallRows([unindexedMember], realTraces);
    expect(rows[0]!.anchors).toEqual([]);
    expect(rows[0]!.phase).toBeNull();
  });

  it("uses an empty trace when no measured trace is found (no throw)", () => {
    const orphan = { ...realMembers[0]!, exposure_id: 999999 };
    const rows = toWaterfallRows([orphan], realTraces);
    expect(rows[0]!.trace.q).toEqual([]);
    expect(rows[0]!.trace.I).toEqual([]);
  });
});

describe("waterfallQDomain", () => {
  it("returns the padded positive q-extent across all rows", () => {
    const rows = toWaterfallRows(realMembers, realTraces);
    const [lo, hi] = waterfallQDomain(rows);
    const allQ = rows.flatMap((r) => r.trace.q).filter((q) => q > 0);
    expect(lo).toBeLessThanOrEqual(Math.min(...allQ));
    expect(hi).toBeGreaterThanOrEqual(Math.max(...allQ));
    expect(lo).toBeGreaterThan(0);
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npm test -- print-waterfall/waterfallModel`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `waterfallModel.ts`**

Create `src/print/waterfall/waterfallModel.ts`:

```ts
import type { SeriesMember, Trace } from "../../api";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

/** One indexed-peak bead position on a member's trace. */
export interface WaterfallAnchor {
  /** Effective-peak id (from peak_ids) — maps 1:1 to PlotPeak.id. */
  id: number;
  q: number;
  /** Member dominant phase (drives bead colour); null = unindexed. */
  phase: string | null;
}

/** One member row, fully resolved for rendering (low→high variable order). */
export interface WaterfallRow {
  key: string;
  label: string;
  trace: Trace;
  phase: string | null;
  state: "indexed" | "form_factor" | "null";
  anchors: WaterfallAnchor[];
  bandHeight: number;
  yOffset: number;
}

/** Dominant phase: first confirmed_phases entry, else confirmed_index.phase, else null. */
function dominantPhase(member: SeriesMember): string | null {
  const snap = member.snapshot;
  if (snap.assignment_state === "form_factor" || snap.assignment_state === "null") return null;
  const cp = snap.confirmed_phases;
  if (cp && cp.length > 0) return cp[0]!.phase;
  return snap.confirmed_index?.phase ?? null;
}

/** Map API members + a trace lookup → rows in input (display) order. */
export function toWaterfallRows(
  members: SeriesMember[],
  tracesById: Record<number, Trace>,
): WaterfallRow[] {
  return members.map((member) => {
    const snap = member.snapshot;
    const state: WaterfallRow["state"] = snap.assignment_state ?? "indexed";
    const phase = dominantPhase(member);
    const trace = (member.exposure_id != null && tracesById[member.exposure_id]) || EMPTY_TRACE;

    const qById = new Map(snap.effective_peaks.map((p) => [p.id, p.q]));
    const anchors: WaterfallAnchor[] =
      state === "indexed" && snap.confirmed_index
        ? snap.confirmed_index.peak_ids
            .filter((id) => qById.has(id))
            .map((id) => ({ id, q: qById.get(id)!, phase }))
        : [];

    return {
      key: String(member.id ?? member.exposure_id ?? member.display_order),
      label: member.label_override ?? (member.exposure_id != null ? `exp ${member.exposure_id}` : ""),
      trace,
      phase,
      state,
      anchors,
      bandHeight: member.band_height ?? 1,
      yOffset: member.y_offset ?? 0,
    };
  });
}

/** Positive q-extent across all rows, padded ×0.97 / ×1.03. */
export function waterfallQDomain(rows: WaterfallRow[]): [number, number] {
  let lo = Infinity;
  let hi = 0;
  for (const r of rows) {
    for (const q of r.trace.q) {
      if (q > 0 && Number.isFinite(q)) {
        if (q < lo) lo = q;
        if (q > hi) hi = q;
      }
    }
  }
  if (!Number.isFinite(lo) || hi <= 0) return [0.01, 1]; // no positive data fallback
  return [lo * 0.97, hi * 1.03];
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `npm test -- print-waterfall/waterfallModel`
Expected: all PASS.

- [ ] **Step 5: tsc + commit**

Run: `npx tsc --noEmit -p tsconfig.build.json` → exit 0.

```bash
git add src/print/waterfall/waterfallModel.ts test/print-waterfall/waterfallModel.test.ts
git commit -m "feat(waterfall): pure view-model (toWaterfallRows + waterfallQDomain)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: `WaterfallChart.tsx` — the renderer

Stacks N single-trace `TracePlot`s bottom-up under one shared `Axis`, with row hot/dim and a sample label per row. The component measures no layout (JSDOM-safe): it positions rows by `bandHeight` ratios into a fixed total height and renders each member's box at a computed pixel offset inside one outer SVG-less flex column; the shared q-axis is one `Axis` inside its own small SVG at the bottom.

To keep it JSDOM-testable and avoid width measurement, use a **fixed internal coordinate height** per the row stack and let the outer container scale via CSS (`width: 100%`, `max-width`). Each member `TracePlot` gets an explicit `width` derived from `maxWidth` and a `height` derived from its `bandHeight` share of the total.

**Files:**
- Create: `src/print/waterfall/WaterfallChart.tsx`
- Test: `test/print-waterfall/WaterfallChart.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/print-waterfall/WaterfallChart.test.tsx`:

```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { WaterfallChart } from "../../src/print/waterfall/WaterfallChart";
import type { WaterfallRow } from "../../src/print/waterfall/waterfallModel";

const ROWS: WaterfallRow[] = [
  {
    key: "a", label: "1:0", phase: "Ia3d", state: "indexed",
    trace: { q: [0.02, 0.05, 0.1], I: [100, 60, 12], sigma: [1, 1, 1] },
    anchors: [{ id: 1, q: 0.05, phase: "Ia3d" }],
    bandHeight: 1, yOffset: 0,
  },
  {
    key: "b", label: "1:1", phase: "Im3m", state: "indexed",
    trace: { q: [0.02, 0.06, 0.12], I: [90, 55, 10], sigma: [1, 1, 1] },
    anchors: [{ id: 2, q: 0.06, phase: "Im3m" }],
    bandHeight: 1, yOffset: 0,
  },
];

describe("WaterfallChart", () => {
  it("renders one row group per member", () => {
    const { container } = render(<WaterfallChart rows={ROWS} />);
    expect(container.querySelectorAll('[data-role="wf-row"]').length).toBe(2);
  });

  it("stacks display-order-0 at the BOTTOM (largest top offset)", () => {
    const { container } = render(<WaterfallChart rows={ROWS} maxWidth={800} />);
    const groups = Array.from(container.querySelectorAll('[data-role="wf-row"]'));
    const topA = Number((groups.find((g) => g.getAttribute("data-key") === "a") as HTMLElement).style.top.replace("px", ""));
    const topB = Number((groups.find((g) => g.getAttribute("data-key") === "b") as HTMLElement).style.top.replace("px", ""));
    expect(topA).toBeGreaterThan(topB); // row "a" (order 0) sits lower on screen
  });

  it("renders a phase-coloured trace line per row (not ink-soft)", () => {
    const { container } = render(<WaterfallChart rows={ROWS} />);
    const lines = container.querySelectorAll('[data-role="trace-line"] path');
    expect(lines.length).toBeGreaterThanOrEqual(2);
    // at least one line is phase-coloured
    const strokes = Array.from(lines).map((p) => p.getAttribute("stroke"));
    expect(strokes.some((s) => s?.includes("oklch"))).toBe(true);
  });

  it("renders a bead per anchor (the peaks layer)", () => {
    const { container } = render(<WaterfallChart rows={ROWS} />);
    // PlotPeaks delegates each peak to PeakGlyph (data-role="peak-glyph").
    expect(container.querySelectorAll('[data-role="peak-glyph"]').length).toBeGreaterThanOrEqual(2);
  });

  it("renders exactly one shared bottom axis", () => {
    const { container } = render(<WaterfallChart rows={ROWS} />);
    expect(container.querySelectorAll('[data-role="axis-bottom"]').length).toBe(1);
  });

  it("renders the sample label per row", () => {
    const { getByText } = render(<WaterfallChart rows={ROWS} />);
    expect(getByText("1:0")).toBeTruthy();
    expect(getByText("1:1")).toBeTruthy();
  });

  it("marks the hovered row hot and dims the others; fires onHoverRow", () => {
    const spy = vi.fn();
    const { container } = render(<WaterfallChart rows={ROWS} onHoverRow={spy} />);
    const rowA = container.querySelector('[data-role="wf-row"][data-key="a"]')!;
    fireEvent.mouseEnter(rowA);
    expect(spy).toHaveBeenCalledWith("a");
    const rowB = container.querySelector('[data-role="wf-row"][data-key="b"]')!;
    expect(rowA.getAttribute("data-hot")).toBe("true");
    expect(rowB.getAttribute("data-dim")).toBe("true");
    fireEvent.mouseLeave(rowA);
    expect(spy).toHaveBeenCalledWith(undefined);
  });

  it("respects a controlled hoveredKey", () => {
    const { container } = render(<WaterfallChart rows={ROWS} hoveredKey="b" />);
    expect(container.querySelector('[data-role="wf-row"][data-key="b"]')!.getAttribute("data-hot")).toBe("true");
    expect(container.querySelector('[data-role="wf-row"][data-key="a"]')!.getAttribute("data-dim")).toBe("true");
  });
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `npm test -- print-waterfall/WaterfallChart`
Expected: FAIL — module not found.

- [ ] **Step 3: Implement `WaterfallChart.tsx`**

Create `src/print/waterfall/WaterfallChart.tsx`:

```tsx
import { useState } from "react";
import { TracePlot, type TraceModel } from "../plot/TracePlot";
import { Axis } from "../plot/Axis";
import { makeAxis } from "../plot/projection";
import type { PlotPeak } from "../plot/marks/PlotPeaks";
import { waterfallQDomain, type WaterfallRow } from "./waterfallModel";

const AXIS_H = 44;          // height of the shared bottom axis strip
const LABEL_W = 56;         // right-margin label gutter
const TOTAL_H = 420;        // fixed internal stack height (CSS scales width)
const Y_HEADROOM = 0.12;    // each box keeps ~12% headroom above its peak

export interface WaterfallChartProps {
  /** Rows low→high (display order); rendered bottom-up. */
  rows: WaterfallRow[];
  xType?: "log" | "linear";
  /** Controlled hot row; falls back to internal hover when omitted. */
  hoveredKey?: string;
  onHoverRow?: (key?: string) => void;
  onHoverQ?: (q?: number) => void;
  /** Fit-to-width ceiling (px). Default 1080 (mockup plate). */
  maxWidth?: number;
  /** Floor below which it stops shrinking (px). Default 560. */
  minWidth?: number;
  /** PLACEMENT ONLY. */
  className?: string;
}

function anchorsToPeaks(row: WaterfallRow): PlotPeak[] {
  return row.anchors.map((a) => ({ id: a.id, q: a.q, source: "auto" as const }));
}

export function WaterfallChart({
  rows,
  xType = "log",
  hoveredKey,
  onHoverRow,
  onHoverQ,
  maxWidth = 1080,
  minWidth = 560,
  className = "",
}: WaterfallChartProps): JSX.Element {
  const [internalHot, setInternalHot] = useState<string | undefined>(undefined);
  const hot = hoveredKey ?? internalHot;

  const qDomain = waterfallQDomain(rows);
  const plotW = Math.max(minWidth, maxWidth) - LABEL_W;

  // Total band weight; each row's pixel height is its share of TOTAL_H.
  const totalWeight = rows.reduce((s, r) => s + Math.max(0, r.bandHeight), 0) || rows.length;

  // Bottom-up: display-order-0 (rows[0]) sits at the bottom. Walk rows in
  // reverse for top→bottom pixel placement.
  let cumulative = 0;
  const placed = [...rows].reverse().map((row) => {
    const h = (Math.max(0, row.bandHeight) / totalWeight) * TOTAL_H;
    const top = cumulative + row.yOffset;
    cumulative += h;
    return { row, top, h };
  });

  const setHot = (key?: string): void => {
    if (hoveredKey === undefined) setInternalHot(key);
    onHoverRow?.(key);
  };

  return (
    <div
      className={`relative w-full ${className}`}
      style={{ maxWidth }}
      data-testid="waterfall"
    >
      {/* Row stack */}
      <div className="relative" style={{ height: TOTAL_H }}>
        {placed.map(({ row, top, h }) => {
          const model: TraceModel = { trace: row.trace, peaks: anchorsToPeaks(row), phase: row.phase };
          const isHot = hot === row.key;
          const isDim = hot !== undefined && !isHot;
          return (
            <div
              key={row.key}
              data-role="wf-row"
              data-key={row.key}
              {...(isHot ? { "data-hot": "true" } : {})}
              {...(isDim ? { "data-dim": "true" } : {})}
              className="absolute left-0 flex items-center"
              style={{ top, height: h, width: "100%", opacity: isDim ? 0.32 : 1 }}
              onMouseEnter={() => setHot(row.key)}
              onMouseLeave={() => setHot(undefined)}
            >
              <TracePlot
                trace={model}
                axes={false}
                xType={xType}
                xDomain={qDomain}
                yHeadroom={Y_HEADROOM}
                height={h}
                width={plotW}
                paperColor="var(--color-plate)"
                show={{ peaks: true, labels: false, band: false }}
                interaction={{
                  onXDomain: () => {},
                  onAddPeak: () => {},
                  onClickPeak: () => {},
                }}
              />
              <div
                data-role="wf-label"
                className={`font-mono text-meta leading-none pl-2 ${isHot ? "text-ink" : "text-ink-soft"}`}
                style={{ width: LABEL_W }}
              >
                {row.label}
              </div>
            </div>
          );
        })}
      </div>

      {/* Shared bottom q-axis */}
      <svg width={plotW} height={AXIS_H} role="img" aria-label="scattering vector q" data-testid="wf-axis">
        <Axis
          axis={makeAxis(qDomain, [0, plotW], xType)}
          orientation="bottom"
          plotWidth={plotW}
          plotHeight={0}
          label="q (Å⁻¹)"
        />
      </svg>
    </div>
  );
}
```

Notes for the implementer:
- `onHoverQ` is threaded for API parity but the per-box `TracePlot` owns its own peak hover/q-readout; if a test for `onHoverQ` is added later, wire it through the box `interaction`. It is intentionally unused now — keep it in the props so the wrapper's contract is stable. (If `tsc`/lint flags the unused prop, prefix with a `void onHoverQ;` no-op or omit it from the destructure until used. Prefer omitting from destructure: do **not** destructure `onHoverQ` if unused, to keep the build clean — but keep it in the Props interface.)
- The per-bead role is `data-role="peak-glyph"` (PlotPeaks wraps with `data-role="plot-peaks"` and delegates each peak to PeakGlyph, which renders `data-role="peak-glyph"`) — already reflected in the Step 1 test.
- `Axis` props are exactly `{ axis, orientation, plotWidth, plotHeight, label? }` (verified). The `orientation="bottom"` ticks render below `plotHeight`; with `plotHeight={0}` the tick labels sit at the top of the 44px strip. If labels clip, increase `AXIS_H` — do not move the axis into the row stack.

- [ ] **Step 4: Run the renderer test**

Run: `npm test -- print-waterfall/WaterfallChart`
Expected: all PASS. (The peak-bead selector `data-role="peak-glyph"` is already correct.)

- [ ] **Step 5: Design-guard + tsc**

Run: `npm run lint:design`
Expected: clean (proves WaterfallChart is placement-only — no inline appearance literals). If it flags a literal, move that appearance into the plot layer or a `ui/` primitive; do **not** add a `print/waterfall/` exemption unless a raw SVG appearance literal is genuinely unavoidable (and if so, mirror the `print/comb/` exemption in `scripts/check-design.mjs` + `test/check-design.test.ts` as a separate committed step).

Run: `npx tsc --noEmit -p tsconfig.build.json` → exit 0.

- [ ] **Step 6: Commit**

```bash
git add src/print/waterfall/WaterfallChart.tsx test/print-waterfall/WaterfallChart.test.tsx
git commit -m "feat(waterfall): WaterfallChart renderer (stacked single-trace TracePlots)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Fixtures, stories, barrel + full build

**Files:**
- Create: `src/print/waterfall/waterfall.fixtures.ts`
- Create: `src/print/waterfall/WaterfallChart.stories.tsx`
- Create: `src/print/waterfall/index.ts`

- [ ] **Step 1: Create the fixture module**

Create `src/print/waterfall/waterfall.fixtures.ts` (plain module — NOT a stories file, so it is not globbed as CSF):

```ts
import { toWaterfallRows, type WaterfallRow } from "./waterfallModel";
import {
  realMembers,
  transitionSeries,
  formFactorMember,
  unindexedMember,
} from "../fixtures/realSeriesMembers";
import { realTraces } from "../fixtures/realTraces";

/** The full real Sample-9 set (Ia3d → Im3m+Ia3d → Im3m + stress members). */
export const FULL: WaterfallRow[] = toWaterfallRows(realMembers, realTraces);

/** The three-member transition only — the canonical hero story. */
export const TRANSITION: WaterfallRow[] = toWaterfallRows(transitionSeries, realTraces);

/** Mixed states: an indexed member, a form-factor member, an unindexed member. */
export const MIXED_STATES: WaterfallRow[] = toWaterfallRows(
  [realMembers[0]!, formFactorMember, unindexedMember],
  realTraces,
);
```

- [ ] **Step 2: Create the barrel**

Create `src/print/waterfall/index.ts`:

```ts
export { WaterfallChart } from "./WaterfallChart";
export type { WaterfallChartProps } from "./WaterfallChart";
export {
  toWaterfallRows, waterfallQDomain,
  type WaterfallRow, type WaterfallAnchor,
} from "./waterfallModel";
```

- [ ] **Step 3: Create the stories**

Create `src/print/waterfall/WaterfallChart.stories.tsx`:

```tsx
import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { WaterfallChart } from "./WaterfallChart";
import { FULL, TRANSITION, MIXED_STATES } from "./waterfall.fixtures";

const meta: Meta<typeof WaterfallChart> = {
  title: "waterfall/WaterfallChart",
  component: WaterfallChart,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof WaterfallChart>;

const frame = "rounded border border-hair bg-plate p-4";

// Wide: the full real Sample-9 set at the plate max width.
export const Wide: Story = {
  render: () => (
    <div className={frame} style={{ width: 980 }}>
      <WaterfallChart rows={FULL} maxWidth={920} />
    </div>
  ),
};

// Narrow: fit-to-width down to the floor (no horizontal scroll).
export const Narrow: Story = {
  render: () => (
    <div className={frame} style={{ width: 600 }}>
      <WaterfallChart rows={TRANSITION} maxWidth={560} minWidth={560} />
    </div>
  ),
};

// Linear q.
export const LinearQ: Story = {
  render: () => (
    <div className={frame} style={{ width: 980 }}>
      <WaterfallChart rows={TRANSITION} maxWidth={920} xType="linear" />
    </div>
  ),
};

// Mixed states: indexed + form-factor (no anchors, neutral line) + unindexed.
export const MixedStates: Story = {
  render: () => (
    <div className={frame} style={{ width: 980 }}>
      <WaterfallChart rows={MIXED_STATES} maxWidth={920} />
    </div>
  ),
};

// Controlled hover: a parent drives the hot row (e.g. from a member list).
export const ControlledHover: Story = {
  render: () => {
    // eslint-disable-next-line react-hooks/rules-of-hooks
    const [hot, setHot] = useState<string | undefined>(TRANSITION[0]?.key);
    return (
      <div className={frame} style={{ width: 980 }}>
        <WaterfallChart rows={TRANSITION} maxWidth={920} hoveredKey={hot} onHoverRow={setHot} />
      </div>
    );
  },
};
```

- [ ] **Step 4: Run the waterfall suite + design guard + tsc**

Run: `npm test -- print-waterfall`
Expected: all PASS (model + chart).

Run: `npm run lint:design`
Expected: clean.

Run: `npx tsc --noEmit -p tsconfig.build.json`
Expected: exit 0.

- [ ] **Step 5: Full build + Storybook build**

Run: `npm run build`
Expected: exit 0 (`lint:design` + `tsc` + `vite build`).

Run: `npm run build-storybook`
Expected: exit 0 (proves the new stories compile + glob cleanly).

- [ ] **Step 6: Commit**

```bash
git add src/print/waterfall/waterfall.fixtures.ts src/print/waterfall/WaterfallChart.stories.tsx src/print/waterfall/index.ts
git commit -m "feat(waterfall): real-data fixtures, stories, barrel

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification (after all tasks)

- [ ] `npm test -- print-plot` and `npm test -- print-waterfall` both green.
- [ ] `npm run build` exit 0; `npm run build-storybook` exit 0.
- [ ] Manual fidelity: `npm run storybook` → `waterfall/WaterfallChart` (Wide) vs the `#waterfall` region of `docs/redesign-mockups/2026-05-29-series-plot.html` — stacked bottom-up, phase-coloured lines, anchor beads, one bottom q-axis, right-margin labels, row hot/dim on hover.
- [ ] Dispatch a final `frontend-reviewer` pass over the whole batch.

## Notes / deferrals (carried from the spec)

- **Migration track** and **heatmap** are explicitly out of scope (track blocked on per-peak reflection identity the member snapshot doesn't carry; heatmap is YAGNI + competes with phase colour).
- Per-peak coexistence anchor colouring deferred — anchors colour by the member's dominant phase.
- The vertical `PhaseStrip` companion already exists; placing it beside the waterfall is the tier-2 Series plate's job, not this renderer's.
