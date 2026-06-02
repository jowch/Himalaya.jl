# Greenfield Comb + Residual Charts Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build two pure SVG renderers — `CombChart` (reflection-comb teeth) and `ResidualChart` (per-phase Δq/q diagnostic) — in a new `src/print/comb/` directory, structurally parallel via a shared `CombScaffold`.

**Architecture:** A shared `CombScaffold` owns the rowed layout (pinned left gutter + horizontally-scrollable always-log q-pane + bottom q-axis) and hands a projection context to a render-prop. `CombChart` draws teeth/caret/ring glyphs into it; `ResidualChart` draws a per-row tolerance band + residual points into it. Both consume a pure `combModel` (view-model types + row assembly + q-domain). No Zustand/Query/URL/API — props in, `onHoverQ` out.

**Tech Stack:** React 18 + TypeScript (strict, `exactOptionalPropertyTypes` ON), Vite, `d3-scale` (via the existing `src/print/plot/projection.ts` helpers), Storybook CSF3 (`@storybook/react-vite`), Vitest + React Testing Library (JSDOM).

**Spec:** [docs/superpowers/specs/2026-06-02-greenfield-comb-residual-charts-design.md](../specs/2026-06-02-greenfield-comb-residual-charts-design.md)

---

## Conventions (read before starting)

- **Working dir:** all `npm`/`git` commands run from `packages/HimalayaUI/frontend/`. (The repo is a git worktree at `.claude/worktrees/greenfield-ui-rebuild/`. Do NOT push; the branch stays unmerged — the whole "The Print" rebuild lands as one piece.)
- **Commit trailer (every commit):**
  ```
  Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
  ```
- **`src/print/comb/` is design-guard-exempt** once Task 1 lands — appearance literals (`var(--color-…)`, `fontSize`, inline numbers) are allowed there, exactly as in `src/print/detector/`. Test files live in `test/print-comb/` and assert via `data-role` / `data-*` / text — **never class strings**.
- **`exactOptionalPropertyTypes` is ON:** never pass an explicit `undefined` for an optional prop. Use conditional spread: `{...(maxWidth !== undefined ? { maxWidth } : {})}`.
- **Test runner:** `npm test -- print-comb/<Name>` (a Vitest pre-tool hook injects `--run`; no watch). Type/build gate is `npx tsc --noEmit -p tsconfig.build.json` (the build tsconfig; the *default* root `tsc` has unrelated pre-existing errors in `test/` and is NOT the gate).
- **Per-task cadence:** write the failing test → run it red → implement → run it green → `npm run lint:design` clean → commit. Final task adds the full `npm run build` + `npm run build-storybook` gate.

## File Structure

```
src/print/comb/
  combModel.ts        # CombTooth, CombSeries, CombRow types + assembleRows() + combQDomain()
  CombScaffold.tsx    # shared rowed layout: pinned gutter + scrollable log-q SVG pane + bottom q-axis; render-prop ctx
  CombChart.tsx       # teeth renderer (Tooth / AbsentCaret / leftover Ring)
  ResidualChart.tsx   # per-row zero-line + tolerance band + residual points (clamp+mark overflow)
  CombChart.stories.tsx
  ResidualChart.stories.tsx
  index.ts            # barrel
test/print-comb/
  combModel.test.ts
  CombScaffold.test.tsx
  CombChart.test.tsx
  ResidualChart.test.tsx
scripts/check-design.mjs   # MODIFY: add "print/comb/" to isExcluded()
```

---

### Task 1: Register `print/comb/` with the design guard

**Files:**
- Modify: `scripts/check-design.mjs:134-144` (the `isExcluded` function + its comment)

- [ ] **Step 1: Add the exclusion prefix**

In `scripts/check-design.mjs`, update the comment and `isExcluded` to add the comb dir:

```js
// src/components/ui/**, src/print/ui/**, src/print/plot/**, src/print/detector/**,
// and src/print/comb/** are excluded (appearance authored — primitives, the trace-plot
// engine, the detector rendering layer, and the comb/residual rendering layer).
function isExcluded(relPath) {
  return (
    relPath.startsWith("components/ui/") ||
    relPath.startsWith("print/ui/") ||
    relPath.startsWith("print/plot/") ||
    relPath.startsWith("print/detector/") ||
    relPath.startsWith("print/comb/")
  );
}
```

- [ ] **Step 2: Verify the exclusion is live**

The module exports `scanContent` (pure). Confirm an appearance literal under `print/comb/` is now exempt:

```bash
node -e "import('./scripts/check-design.mjs').then(m => { const v = m.scanContent('print/comb/x.tsx', 'const a = \"text-[10px]\"'); if (v.length) { console.error('NOT excluded', v); process.exit(1); } console.log('print/comb excluded OK'); })"
```
Expected: `print/comb excluded OK`

- [ ] **Step 3: Confirm the existing design lint still passes, and update any check-design unit test**

```bash
npm run lint:design
grep -rl "isExcluded\|scanContent\|print/detector" test/ 2>/dev/null
```
Expected: `lint:design` exits 0. The `grep` lists any test that drives the design guard. If one exists and asserts the exclusion list (e.g. it checks `print/detector/` is excluded), add a parallel `print/comb/` assertion and run that single file with `npm test -- <path>`; if the grep finds nothing, there is no such test and you can skip it.

- [ ] **Step 4: Commit**

```bash
git add scripts/check-design.mjs
git commit -m "Design guard: exempt src/print/comb/ (comb/residual renderer layer)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: `combModel.ts` — view-model + row assembly + q-domain

**Files:**
- Create: `src/print/comb/combModel.ts`
- Test: `test/print-comb/combModel.test.ts`

- [ ] **Step 1: Write the failing test**

Create `test/print-comb/combModel.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { assembleRows, combQDomain, type CombSeries } from "../../src/print/comb/combModel";

const pn3m: CombSeries = {
  phase: "Pn3m",
  color: "var(--a)",
  teeth: [
    { q: 0.10, label: "√2", observed: true },
    { q: 0.1225, label: "√3", observed: false },
  ],
};
const im3m: CombSeries = {
  phase: "Im3m",
  color: "var(--b)",
  teeth: [{ q: 0.12, label: "√2", observed: true }],
};

describe("assembleRows", () => {
  it("orders assigned rows, then the preview row, then the leftover row", () => {
    const rows = assembleRows([pn3m, im3m], im3m, [0.2]);
    expect(rows.map((r) => r.kind)).toEqual(["assigned", "assigned", "preview", "leftover"]);
  });

  it("omits the preview row when there is no hovered series", () => {
    const rows = assembleRows([pn3m], undefined, [0.2]);
    expect(rows.some((r) => r.kind === "preview")).toBe(false);
  });

  it("omits the leftover row when there are no leftover peaks", () => {
    const rows = assembleRows([pn3m], undefined, []);
    expect(rows.some((r) => r.kind === "leftover")).toBe(false);
  });

  it("carries the leftover q-values on the leftover row", () => {
    const rows = assembleRows([pn3m], undefined, [0.2, 0.25]);
    const leftover = rows.find((r) => r.kind === "leftover");
    expect(leftover && leftover.kind === "leftover" && leftover.qs).toEqual([0.2, 0.25]);
  });
});

describe("combQDomain", () => {
  it("spans every q (teeth + leftover) with ~10% pad and ignores non-positive", () => {
    const rows = assembleRows([pn3m], undefined, [0.2, -1, 0]);
    const [lo, hi] = combQDomain(rows);
    expect(lo).toBeCloseTo(0.10 * 0.9, 6);
    expect(hi).toBeCloseTo(0.20 * 1.1, 6);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

```bash
npm test -- print-comb/combModel
```
Expected: FAIL — cannot resolve `../../src/print/comb/combModel` (module does not exist yet).

- [ ] **Step 3: Implement `combModel.ts`**

Create `src/print/comb/combModel.ts`:

```ts
import { positiveExtent } from "../plot/projection";

/** One predicted Bragg reflection in a phase's comb. */
export interface CombTooth {
  /** Predicted reflection q (Å⁻¹). */
  q: number;
  /** hkl or √N label shown above the tooth, e.g. "√6". */
  label: string;
  /** A claimed observed peak exists → solid tooth; otherwise a predicted-absent caret. */
  observed: boolean;
  /** Fractional residual (q_obs − q_pred)/q_pred for the claimed peak. Used by ResidualChart. */
  residual?: number;
}

/** One assigned (or hovered-preview) phase's comb row. `color` is a token ref
 *  (phaseColor(phase)) supplied by the composite. */
export interface CombSeries {
  phase: string;
  color: string;
  /** Comb gutter line 2, e.g. "a = 197 Å" (from lattice_d). */
  latticeLabel?: string;
  /** Residual gutter line 2 value, e.g. 0.998 (from r_squared). */
  rSquared?: number;
  teeth: CombTooth[];
}

/** A laid-out row. `assigned`/`preview` carry a series; `leftover` carries raw q-values. */
export type CombRow =
  | { kind: "assigned"; series: CombSeries }
  | { kind: "preview"; series: CombSeries }
  | { kind: "leftover"; qs: number[] };

/** Ordered rows: assigned phases, then the optional hovered-preview row, then the
 *  optional leftover (unindexed-peaks) row. ResidualChart passes `leftover: []`
 *  so it never gets a leftover row. */
export function assembleRows(
  assigned: CombSeries[],
  hovered: CombSeries | undefined,
  leftover: number[],
): CombRow[] {
  const rows: CombRow[] = assigned.map((series) => ({ kind: "assigned", series }));
  if (hovered) rows.push({ kind: "preview", series: hovered });
  if (leftover.length > 0) rows.push({ kind: "leftover", qs: leftover });
  return rows;
}

/** Padded log-q domain over every q in the rows (teeth + leftover). ~10% pad each
 *  side (mockup convention). Non-positive q are ignored by positiveExtent. */
export function combQDomain(rows: CombRow[]): [number, number] {
  const qs: number[] = [];
  for (const row of rows) {
    if (row.kind === "leftover") qs.push(...row.qs);
    else for (const t of row.series.teeth) qs.push(t.q);
  }
  const [lo, hi] = positiveExtent(qs, [0.01, 0.2]);
  return [lo * 0.9, hi * 1.1];
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
npm test -- print-comb/combModel
```
Expected: PASS (5 tests).

- [ ] **Step 5: Commit**

```bash
git add src/print/comb/combModel.ts test/print-comb/combModel.test.ts
git commit -m "comb: combModel — view-model types, row assembly, padded log-q domain

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: `CombScaffold.tsx` — shared rowed layout

**Files:**
- Create: `src/print/comb/CombScaffold.tsx`
- Test: `test/print-comb/CombScaffold.test.tsx`

The scaffold renders a flex row: a fixed-width **pinned HTML gutter** (one block per row, phase title + sub-line) beside an **`overflow-x-auto` pane** holding the SVG. It computes the always-log q-axis, lays out rows at a constant `ROW_H`, draws faint q-gridlines + a bottom q-axis (via `axisTicks`) + per-row baselines (dashed for `preview`), and calls `children(ctx)` to inject chart-specific glyphs. The SVG width has a floor (`MIN_PX_PER_DECADE × decades`, min `MIN_PLOT_W`) so teeth never crowd; when the container is narrower, the pane scrolls.

- [ ] **Step 1: Write the failing test**

Create `test/print-comb/CombScaffold.test.tsx`:

```tsx
import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { CombScaffold, MIN_PLOT_W, type ScaffoldRow } from "../../src/print/comb/CombScaffold";

const rows: ScaffoldRow[] = [
  { gutterTitle: "Pn3m", gutterSub: "a = 197 Å" },
  { gutterTitle: "Im3m", gutterSub: "a = 142 Å", preview: true },
];

describe("CombScaffold", () => {
  it("renders one pinned gutter block per row, with the title and sub-line", () => {
    const { container, getByText } = render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} ariaLabel="comb">{() => null}</CombScaffold>,
    );
    expect(container.querySelectorAll('[data-role="gutter-row"]').length).toBe(2);
    expect(getByText("Pn3m")).toBeTruthy();
    expect(getByText("a = 197 Å")).toBeTruthy();
  });

  it("draws a baseline per row; the preview row's baseline is dashed", () => {
    const { container } = render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} ariaLabel="comb">{() => null}</CombScaffold>,
    );
    const baselines = container.querySelectorAll('[data-role="row-baseline"]');
    expect(baselines.length).toBe(2);
    const preview = container.querySelector('[data-role="row-baseline"][data-preview="true"]')!;
    expect(preview.getAttribute("stroke-dasharray")).toBeTruthy();
  });

  it("draws labelled q-axis ticks (major/mid), skipping minor labels", () => {
    const { container } = render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} ariaLabel="comb">{() => null}</CombScaffold>,
    );
    expect(container.querySelectorAll('[data-role="q-tick"]').length).toBeGreaterThan(0);
  });

  it("enforces a min plot width so teeth never crowd (narrow container → wide SVG, pane scrolls)", () => {
    const { container } = render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} maxWidth={320} ariaLabel="comb">{() => null}</CombScaffold>,
    );
    const svg = container.querySelector("svg")!;
    expect(Number(svg.getAttribute("data-plot-w"))).toBeGreaterThanOrEqual(MIN_PLOT_W);
  });

  it("invokes the render-prop with a working log-q projection and a baselineY mapper", () => {
    let seen: { x0: number; x1: number; rowCount: number } | null = null;
    render(
      <CombScaffold rows={rows} xDomain={[0.05, 0.25]} ariaLabel="comb">
        {(ctx) => { seen = { x0: ctx.x.to(0.05), x1: ctx.x.to(0.25), rowCount: ctx.rowCount }; return null; }}
      </CombScaffold>,
    );
    expect(seen!.x1).toBeGreaterThan(seen!.x0); // q increases left→right
    expect(seen!.rowCount).toBe(2);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

```bash
npm test -- print-comb/CombScaffold
```
Expected: FAIL — cannot resolve `../../src/print/comb/CombScaffold`.

- [ ] **Step 3: Implement `CombScaffold.tsx`**

Create `src/print/comb/CombScaffold.tsx`:

```tsx
import { type ReactNode } from "react";
import { makeAxis, axisTicks, type Axis1D } from "../plot/projection";

export const GUTTER_W = 96;
export const ROW_H = 48;
export const TOP_PAD = 10;
export const AXIS_H = 28;
export const PLOT_PAD_X = 14;
export const MIN_PX_PER_DECADE = 240;
export const MIN_PLOT_W = 320;
const DEFAULT_MAX_W = 760;

export interface ScaffoldRow {
  gutterTitle: string;
  gutterSub?: string;
  preview?: boolean;
}

/** What the scaffold hands the render-prop so a chart can place glyphs. */
export interface ScaffoldCtx {
  /** Log-q axis mapping data → pixel across the plot area. */
  x: Axis1D;
  plotW: number;
  padX: number;
  rowH: number;
  topPad: number;
  rowCount: number;
  /** y of row i's baseline (teeth/rings/zero-line sit here; teeth grow upward). */
  baselineY: (i: number) => number;
  rowsBottom: number;
}

interface Props {
  rows: ScaffoldRow[];
  xDomain: [number, number];
  maxWidth?: number;
  ariaLabel: string;
  children: (ctx: ScaffoldCtx) => ReactNode;
}

/** Compact q tick label: 0.05, 0.1, 0.12, 0.2 … */
function formatTick(q: number): string {
  return String(Number(q.toPrecision(2)));
}

export function CombScaffold({ rows, xDomain, maxWidth, ariaLabel, children }: Props): JSX.Element {
  const maxW = maxWidth ?? DEFAULT_MAX_W;
  const [lo, hi] = xDomain;
  // Always log. Floor the plot width so √N teeth never crowd below legibility;
  // when the container is narrower than gutter+floor, the pane scrolls.
  const decades = Math.max(0.5, Math.log10(hi / lo));
  const floorPlotW = Math.max(MIN_PLOT_W, Math.ceil(decades * MIN_PX_PER_DECADE));
  const fitPlotW = Math.max(0, maxW - GUTTER_W - 2 * PLOT_PAD_X);
  const plotW = Math.max(floorPlotW, fitPlotW);
  const svgW = plotW + 2 * PLOT_PAD_X;
  const rowsBottom = TOP_PAD + rows.length * ROW_H;
  const svgH = rowsBottom + AXIS_H;

  const x = makeAxis([lo, hi], [PLOT_PAD_X, PLOT_PAD_X + plotW], "log");
  const baselineY = (i: number): number => TOP_PAD + i * ROW_H + ROW_H - 14;
  const ctx: ScaffoldCtx = {
    x, plotW, padX: PLOT_PAD_X, rowH: ROW_H, topPad: TOP_PAD,
    rowCount: rows.length, baselineY, rowsBottom,
  };
  const ticks = axisTicks(x, 6);

  return (
    <div className="flex items-stretch" style={{ maxWidth: maxW }} data-testid="comb-scaffold">
      {/* Pinned gutter (HTML, never scrolled) — vertically aligned to the SVG rows. */}
      <div className="shrink-0 select-none" style={{ width: GUTTER_W }}>
        {rows.map((r, i) => (
          <div
            key={i}
            data-role="gutter-row"
            data-preview={r.preview ? "true" : undefined}
            className="flex flex-col justify-end pr-2"
            style={{ height: ROW_H, marginTop: i === 0 ? TOP_PAD : 0, paddingBottom: 8 }}
          >
            <div className={`font-mono text-data leading-none ${r.preview ? "text-ink-faint" : "text-ink"}`}>
              {r.gutterTitle}
            </div>
            {r.gutterSub ? (
              <div className="font-mono text-meta text-ink-faint leading-none mt-0.5">{r.gutterSub}</div>
            ) : null}
          </div>
        ))}
      </div>

      {/* Scrollable q-pane. */}
      <div className="overflow-x-auto" style={{ maxWidth: maxW - GUTTER_W }}>
        <svg width={svgW} height={svgH} role="img" aria-label={ariaLabel} data-plot-w={plotW}>
          {ticks.map((t, i) => {
            const px = x.to(t.value);
            return (
              <g key={i} data-role="q-tick" data-tick-kind={t.kind}>
                <line
                  x1={px} y1={TOP_PAD} x2={px} y2={rowsBottom}
                  stroke="var(--color-hair)" strokeWidth={1}
                  opacity={t.kind === "major" ? 0.9 : 0.5}
                />
                {t.kind !== "minor" ? (
                  <text
                    x={px} y={rowsBottom + 16} textAnchor="middle"
                    className="font-mono" fontSize={9.5} fill="var(--color-ink-faint)"
                  >
                    {formatTick(t.value)}
                  </text>
                ) : null}
              </g>
            );
          })}
          {rows.map((r, i) => (
            <line
              key={i}
              data-role="row-baseline"
              data-preview={r.preview ? "true" : undefined}
              x1={PLOT_PAD_X} y1={baselineY(i)} x2={PLOT_PAD_X + plotW} y2={baselineY(i)}
              stroke="var(--color-hair-strong)" strokeWidth={1}
              {...(r.preview ? { strokeDasharray: "3 3" } : {})}
              opacity={0.8}
            />
          ))}
          {children(ctx)}
        </svg>
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

```bash
npm test -- print-comb/CombScaffold
```
Expected: PASS (5 tests).

- [ ] **Step 5: Lint + commit**

```bash
npm run lint:design
git add src/print/comb/CombScaffold.tsx test/print-comb/CombScaffold.test.tsx
git commit -m "comb: CombScaffold — pinned gutter + scrollable log-q pane + bottom axis

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: `CombChart.tsx` — teeth renderer

**Files:**
- Create: `src/print/comb/CombChart.tsx`
- Test: `test/print-comb/CombChart.test.tsx`

Draws, per row: solid **teeth** (stem + filled cap + hkl label) for observed reflections, dashed-stem **carets** for predicted-absent, and **rings** for leftover peaks. Hover/q-link emphasis is wider stroke + a ring in the glyph's **own color** (never the terracotta accent). Hovering any glyph fires `onHoverQ(q)`.

- [ ] **Step 1: Write the failing test**

Create `test/print-comb/CombChart.test.tsx`:

```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { CombChart } from "../../src/print/comb/CombChart";
import type { CombSeries } from "../../src/print/comb/combModel";

const PN3M: CombSeries = {
  phase: "Pn3m",
  color: "oklch(0.570 0.150 58)",
  latticeLabel: "a = 197 Å",
  teeth: [
    { q: 0.0712, label: "√2", observed: true },
    { q: 0.0872, label: "√3", observed: true },
    { q: 0.1126, label: "√5", observed: false }, // predicted-absent
  ],
};

describe("CombChart", () => {
  it("renders one comb-row per assigned series plus a leftover row", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[0.2]} />);
    expect(container.querySelectorAll('[data-role="comb-row"]').length).toBe(2);
    expect(container.querySelector('[data-row-kind="leftover"]')).toBeTruthy();
  });

  it("draws a solid tooth (stem + cap + label) for observed reflections", () => {
    const { container, getByText } = render(<CombChart assigned={[PN3M]} leftover={[]} />);
    expect(container.querySelectorAll('[data-role="tooth"]').length).toBe(2);
    expect(container.querySelector('[data-role="tooth-cap"]')).toBeTruthy();
    expect(getByText("√2")).toBeTruthy();
    // Tooth strokes the phase colour, not the accent.
    expect(container.querySelector('[data-role="tooth-stem"]')!.getAttribute("stroke")).toBe("oklch(0.570 0.150 58)");
  });

  it("draws a dashed caret (not a tooth) for a predicted-absent reflection", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[]} />);
    expect(container.querySelectorAll('[data-role="caret"]').length).toBe(1);
    const stem = container.querySelector('[data-role="caret-stem"]')!;
    expect(stem.getAttribute("stroke-dasharray")).toBeTruthy();
  });

  it("draws a ring per leftover peak", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[0.2, 0.25]} />);
    expect(container.querySelectorAll('[data-role="leftover-ring"]').length).toBe(2);
  });

  it("renders a dashed preview row when a hovered candidate is supplied", () => {
    const { container } = render(<CombChart assigned={[PN3M]} hovered={PN3M} leftover={[]} />);
    expect(container.querySelector('[data-role="row-baseline"][data-preview="true"]')).toBeTruthy();
  });

  it("a hovered tooth goes hot, keeping its OWN colour (no accent recolor) + a ring", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[]} hoveredQ={0.0712} />);
    const hot = container.querySelector('[data-role="tooth"][data-hot="true"]')!;
    expect(hot.getAttribute("data-q")).toBe("0.0712");
    const stem = hot.querySelector('[data-role="tooth-stem"]')!;
    expect(stem.getAttribute("stroke")).toBe("oklch(0.570 0.150 58)");
    expect(stem.getAttribute("stroke")).not.toBe("var(--color-accent)");
    expect(hot.querySelector('[data-role="tooth-ring"]')).toBeTruthy();
  });

  it("fires onHoverQ(q) on tooth enter and onHoverQ(undefined) on leave", () => {
    const spy = vi.fn();
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[]} onHoverQ={spy} />);
    const tooth = container.querySelector('[data-role="tooth"]')!;
    fireEvent.mouseEnter(tooth);
    expect(spy).toHaveBeenCalledWith(0.0712);
    fireEvent.mouseLeave(tooth);
    expect(spy).toHaveBeenCalledWith(undefined);
  });

  it("enforces the min plot width at a narrow maxWidth (pane scrolls)", () => {
    const { container } = render(<CombChart assigned={[PN3M]} leftover={[]} maxWidth={320} />);
    expect(Number(container.querySelector("svg")!.getAttribute("data-plot-w"))).toBeGreaterThanOrEqual(320);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

```bash
npm test -- print-comb/CombChart
```
Expected: FAIL — cannot resolve `../../src/print/comb/CombChart`.

- [ ] **Step 3: Implement `CombChart.tsx`**

Create `src/print/comb/CombChart.tsx`:

```tsx
import { CombScaffold, type ScaffoldRow, type ScaffoldCtx } from "./CombScaffold";
import { assembleRows, combQDomain, type CombSeries, type CombRow } from "./combModel";

const TOL = 1e-6;
const LEFTOVER_COLOR = "var(--color-ink-faint)";

interface Props {
  assigned: CombSeries[];
  hovered?: CombSeries;
  leftover: number[];
  /** Incoming q-link; the matching glyph across all rows lights hot. */
  hoveredQ?: number;
  /** Fired on glyph enter (q) / leave (undefined). Omit → inert glyphs. */
  onHoverQ?: (q?: number) => void;
  maxWidth?: number;
  className?: string;
}

export function CombChart({ assigned, hovered, leftover, hoveredQ, onHoverQ, maxWidth, className }: Props): JSX.Element {
  const rows = assembleRows(assigned, hovered, leftover);
  const xDomain = combQDomain(rows);
  const scaffoldRows: ScaffoldRow[] = rows.map(rowToGutter);
  return (
    <div className={className}>
      <CombScaffold
        rows={scaffoldRows}
        xDomain={xDomain}
        ariaLabel="reflection comb"
        {...(maxWidth !== undefined ? { maxWidth } : {})}
      >
        {(ctx) => rows.map((row, i) => renderRow(row, i, ctx, hoveredQ, onHoverQ))}
      </CombScaffold>
    </div>
  );
}

function rowToGutter(row: CombRow): ScaffoldRow {
  if (row.kind === "leftover") return { gutterTitle: "leftover", gutterSub: `${row.qs.length} unindexed` };
  const base: ScaffoldRow = {
    gutterTitle: row.series.phase,
    ...(row.series.latticeLabel ? { gutterSub: row.series.latticeLabel } : {}),
  };
  return row.kind === "preview" ? { ...base, preview: true } : base;
}

function renderRow(
  row: CombRow,
  i: number,
  ctx: ScaffoldCtx,
  hoveredQ: number | undefined,
  onHoverQ?: (q?: number) => void,
): JSX.Element {
  const y0 = ctx.baselineY(i);
  const toothH = ctx.rowH - 26;
  const hov = (q: number): { onMouseEnter: () => void; onMouseLeave: () => void } | Record<string, never> =>
    onHoverQ ? { onMouseEnter: () => onHoverQ(q), onMouseLeave: () => onHoverQ(undefined) } : {};
  const isHot = (q: number): boolean => hoveredQ !== undefined && Math.abs(q - hoveredQ) <= TOL;

  if (row.kind === "leftover") {
    return (
      <g key={i} data-role="comb-row" data-row-kind="leftover">
        {row.qs.map((q, j) => {
          const hot = isHot(q);
          return (
            <g
              key={j}
              data-role="leftover-ring"
              data-q={q}
              {...(hot ? { "data-hot": "true" } : {})}
              style={{ cursor: onHoverQ ? "pointer" : "default" }}
              {...hov(q)}
            >
              <circle
                cx={ctx.x.to(q)} cy={y0} r={hot ? 5 : 3.6}
                fill="none" stroke={LEFTOVER_COLOR} strokeWidth={hot ? 2.4 : 1.6}
              />
            </g>
          );
        })}
      </g>
    );
  }

  const { color } = row.series;
  return (
    <g key={i} data-role="comb-row" data-row-kind={row.kind}>
      {row.series.teeth.map((t, j) => {
        const px = ctx.x.to(t.q);
        const hot = isHot(t.q);
        if (t.observed) {
          const apexY = y0 - toothH;
          return (
            <g
              key={j}
              data-role="tooth"
              data-q={t.q}
              {...(hot ? { "data-hot": "true" } : {})}
              style={{ cursor: onHoverQ ? "pointer" : "default" }}
              {...hov(t.q)}
            >
              <line data-role="tooth-stem" x1={px} y1={y0} x2={px} y2={apexY} stroke={color} strokeWidth={hot ? 2.8 : 2} />
              <circle data-role="tooth-cap" cx={px} cy={apexY} r={hot ? 5 : 2.6} fill={color} />
              {hot ? (
                <circle data-role="tooth-ring" cx={px} cy={apexY} r={8} fill="none" stroke={color} strokeWidth={1.4} opacity={0.6} />
              ) : null}
              <text
                data-role="tooth-mlabel" x={px} y={apexY - 7} textAnchor="middle"
                className="font-mono" fontSize={9.5} fontWeight={700} fill={color}
              >
                {t.label}
              </text>
            </g>
          );
        }
        // Predicted-absent: dashed faint stem (~72%) + hollow caret.
        const ah = toothH * 0.72;
        const apexY = y0 - ah;
        const w = hot ? 2 : 1.3;
        return (
          <g
            key={j}
            data-role="caret"
            data-q={t.q}
            {...(hot ? { "data-hot": "true" } : {})}
            style={{ cursor: onHoverQ ? "pointer" : "default" }}
            {...hov(t.q)}
          >
            <line
              data-role="caret-stem" x1={px} y1={y0} x2={px} y2={apexY}
              stroke={LEFTOVER_COLOR} strokeWidth={w} strokeDasharray="1.5 1.8"
            />
            <path
              data-role="caret-cap" d={`M${px - 3.2} ${apexY} L${px} ${apexY - 4} L${px + 3.2} ${apexY}`}
              fill="none" stroke={LEFTOVER_COLOR} strokeWidth={w}
            />
          </g>
        );
      })}
    </g>
  );
}
```

- [ ] **Step 4: Write the shared fixtures module**

Create `src/print/comb/comb.fixtures.ts`. (Fixtures must NOT live in a `*.stories.tsx` file — Storybook CSF treats every named export there as a story, so data objects would pollute/break the story list. Both story files import from here.)

```ts
import { phaseColor } from "../../phases";
import type { CombSeries } from "./combModel";

// A Pn3m series (amber) with one predicted-absent reflection, and an Im3m series
// (sage) whose √8 residual overflows the ±3% residual domain. Plus leftover peaks.
export const PN3M: CombSeries = {
  phase: "Pn3m", color: phaseColor("Pn3m"), latticeLabel: "a = 197 Å", rSquared: 0.998,
  teeth: [
    { q: 0.0712, label: "√2", observed: true, residual: 0.004 },
    { q: 0.0872, label: "√3", observed: true, residual: -0.006 },
    { q: 0.1007, label: "√4", observed: true, residual: 0.001 },
    { q: 0.1126, label: "√5", observed: false },
    { q: 0.1234, label: "√6", observed: true, residual: 0.012 },
    { q: 0.1333, label: "√7", observed: false },
  ],
};
export const IM3M: CombSeries = {
  phase: "Im3m", color: phaseColor("Im3m"), latticeLabel: "a = 142 Å", rSquared: 0.991,
  teeth: [
    { q: 0.0626, label: "√2", observed: true, residual: -0.003 },
    { q: 0.0885, label: "√4", observed: true, residual: 0.007 },
    { q: 0.1084, label: "√6", observed: false },
    { q: 0.1252, label: "√8", observed: true, residual: 0.041 }, // overflow demo
  ],
};
export const LEFTOVER: number[] = [0.156, 0.205];
```

- [ ] **Step 5: Write the story**

Create `src/print/comb/CombChart.stories.tsx`:

```tsx
import type { Meta, StoryObj } from "@storybook/react-vite";
import { useState } from "react";
import { CombChart } from "./CombChart";
import { ResidualChart } from "./ResidualChart";
import { PN3M, IM3M, LEFTOVER } from "./comb.fixtures";

const meta: Meta<typeof CombChart> = { title: "comb/CombChart", component: CombChart };
export default meta;
type Story = StoryObj<typeof CombChart>;

const frame = "rounded border border-hair bg-plate p-3";

export const Wide: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><CombChart assigned={[PN3M, IM3M]} leftover={LEFTOVER} maxWidth={760} /></div>,
};

export const NarrowScrolling: Story = {
  render: () => <div className={frame} style={{ width: 360 }}><CombChart assigned={[PN3M, IM3M]} leftover={LEFTOVER} maxWidth={360} /></div>,
};

export const WithPreviewRow: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><CombChart assigned={[PN3M]} hovered={IM3M} leftover={LEFTOVER} maxWidth={760} /></div>,
};

// Comb + residual sharing a hoveredQ — proves the q-link and the tracked layout.
export const LinkedWithResidual: Story = {
  render: () => {
    const [hoveredQ, setHoveredQ] = useState<number | undefined>(undefined);
    const linkProps = { ...(hoveredQ !== undefined ? { hoveredQ } : {}), onHoverQ: setHoveredQ };
    return (
      <div className="flex flex-col gap-3" style={{ width: 760 }}>
        <div className={frame}><CombChart assigned={[PN3M, IM3M]} leftover={LEFTOVER} maxWidth={760} {...linkProps} /></div>
        <div className={frame}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={760} {...linkProps} /></div>
      </div>
    );
  },
};
```

- [ ] **Step 6: Run tests + lint**

```bash
npm test -- print-comb/CombChart && npm run lint:design
```
Expected: 8 CombChart tests PASS; `lint:design` exits 0. Do NOT run the `tsc -p tsconfig.build.json` gate yet — `CombChart.stories.tsx` imports `ResidualChart`, which lands in Task 5, so the type gate can only pass once Task 5 is done (it runs there).

- [ ] **Step 7: Commit**

```bash
git add src/print/comb/CombChart.tsx src/print/comb/comb.fixtures.ts src/print/comb/CombChart.stories.tsx test/print-comb/CombChart.test.tsx
git commit -m "comb: CombChart — teeth / absent-caret / leftover-ring renderer + fixtures + stories

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: `ResidualChart.tsx` — per-row Δq/q diagnostic

**Files:**
- Create: `src/print/comb/ResidualChart.tsx`
- Test: `test/print-comb/ResidualChart.test.tsx`

Tracks CombChart's scaffold (assigned + preview rows, **no leftover**). Per row: a faint ±2.2% tolerance band centred on the row baseline (the zero line), and a residual dot per observed reflection at `(q, −residual·scale)`. The y-domain is fixed symmetric (±3%); a residual beyond it **clamps to the edge and is marked** with an open dot. The gutter sub-line shows `R²`.

- [ ] **Step 1: Write the failing test**

Create `test/print-comb/ResidualChart.test.tsx`:

```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ResidualChart } from "../../src/print/comb/ResidualChart";
import type { CombSeries } from "../../src/print/comb/combModel";

const PN3M: CombSeries = {
  phase: "Pn3m",
  color: "oklch(0.570 0.150 58)",
  rSquared: 0.998,
  teeth: [
    { q: 0.0712, label: "√2", observed: true, residual: 0.010 },  // positive → above baseline
    { q: 0.0872, label: "√3", observed: true, residual: -0.008 }, // negative → below baseline
    { q: 0.1126, label: "√5", observed: false },                  // no residual → no point
    { q: 0.1234, label: "√6", observed: true, residual: 0.041 },  // beyond ±3% → overflow
  ],
};

describe("ResidualChart", () => {
  it("renders one resid-row per assigned/preview series and NO leftover row", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} hovered={PN3M} />);
    expect(container.querySelectorAll('[data-role="resid-row"]').length).toBe(2);
    expect(container.querySelector('[data-row-kind="leftover"]')).toBeNull();
  });

  it("shows the R² in the gutter", () => {
    const { getByText } = render(<ResidualChart assigned={[PN3M]} />);
    expect(getByText(/R².*0\.998/)).toBeTruthy();
  });

  it("draws a tolerance band per row", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} />);
    expect(container.querySelectorAll('[data-role="tolerance-band"]').length).toBe(1);
  });

  it("plots a point only for observed reflections with a residual", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} />);
    expect(container.querySelectorAll('[data-role="resid-point"]').length).toBe(3); // √2 √3 √6, not √5
  });

  it("places a positive residual above the baseline and a negative one below", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} />);
    const baseY = Number(container.querySelector('[data-role="row-baseline"]')!.getAttribute("y1"));
    const pos = container.querySelector('[data-role="resid-point"][data-q="0.0712"] circle')!;
    const neg = container.querySelector('[data-role="resid-point"][data-q="0.0872"] circle')!;
    expect(Number(pos.getAttribute("cy"))).toBeLessThan(baseY);    // above = smaller y
    expect(Number(neg.getAttribute("cy"))).toBeGreaterThan(baseY); // below = larger y
  });

  it("marks an out-of-domain residual as overflow (open dot)", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} />);
    const overflow = container.querySelector('[data-role="resid-point"][data-overflow="true"] circle')!;
    expect(overflow.getAttribute("data-q") ?? overflow.parentElement!.getAttribute("data-q")).toBe("0.1234");
    expect(overflow.getAttribute("fill")).toBe("none");
  });

  it("a hovered point goes hot keeping its own colour (no accent)", () => {
    const { container } = render(<ResidualChart assigned={[PN3M]} hoveredQ={0.0712} />);
    const hot = container.querySelector('[data-role="resid-point"][data-hot="true"]')!;
    expect(hot.getAttribute("data-q")).toBe("0.0712");
    const circle = hot.querySelector("circle")!;
    expect(circle.getAttribute("stroke")).toBe("oklch(0.570 0.150 58)");
    expect(circle.getAttribute("stroke")).not.toBe("var(--color-accent)");
  });

  it("fires onHoverQ on point enter/leave", () => {
    const spy = vi.fn();
    const { container } = render(<ResidualChart assigned={[PN3M]} onHoverQ={spy} />);
    const point = container.querySelector('[data-role="resid-point"]')!;
    fireEvent.mouseEnter(point);
    expect(spy).toHaveBeenCalledWith(0.0712);
    fireEvent.mouseLeave(point);
    expect(spy).toHaveBeenCalledWith(undefined);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

```bash
npm test -- print-comb/ResidualChart
```
Expected: FAIL — cannot resolve `../../src/print/comb/ResidualChart`.

- [ ] **Step 3: Implement `ResidualChart.tsx`**

Create `src/print/comb/ResidualChart.tsx`:

```tsx
import { CombScaffold, ROW_H, type ScaffoldRow, type ScaffoldCtx } from "./CombScaffold";
import { assembleRows, combQDomain, type CombSeries, type CombRow } from "./combModel";

const TOL = 1e-6;
const RESID_DOMAIN = 0.03; // fixed symmetric y-domain (±3%)
const BAND = 0.022;        // tolerance band drawn at ±2.2%
const HALF_SPAN = ROW_H / 2 - 9; // px from baseline to the row's ±RESID_DOMAIN edge

interface Props {
  assigned: CombSeries[];
  hovered?: CombSeries;
  hoveredQ?: number;
  onHoverQ?: (q?: number) => void;
  maxWidth?: number;
  className?: string;
}

export function ResidualChart({ assigned, hovered, hoveredQ, onHoverQ, maxWidth, className }: Props): JSX.Element {
  // No leftover row in the residual view (unindexed peaks have no predicted q).
  const rows = assembleRows(assigned, hovered, []);
  const xDomain = combQDomain(rows);
  const scaffoldRows: ScaffoldRow[] = rows.map(rowToGutter);
  return (
    <div className={className}>
      <CombScaffold
        rows={scaffoldRows}
        xDomain={xDomain}
        ariaLabel="indexing residuals"
        {...(maxWidth !== undefined ? { maxWidth } : {})}
      >
        {(ctx) => rows.map((row, i) => renderResidRow(row, i, ctx, hoveredQ, onHoverQ))}
      </CombScaffold>
    </div>
  );
}

function rowToGutter(row: CombRow): ScaffoldRow {
  if (row.kind === "leftover") return { gutterTitle: "leftover" }; // unreachable (no leftover passed)
  const sub = row.series.rSquared !== undefined ? `R² ${row.series.rSquared.toFixed(3)}` : undefined;
  const base: ScaffoldRow = { gutterTitle: row.series.phase, ...(sub ? { gutterSub: sub } : {}) };
  return row.kind === "preview" ? { ...base, preview: true } : base;
}

function renderResidRow(
  row: CombRow,
  i: number,
  ctx: ScaffoldCtx,
  hoveredQ: number | undefined,
  onHoverQ?: (q?: number) => void,
): JSX.Element {
  if (row.kind === "leftover") return <g key={i} data-role="resid-row" data-row-kind="leftover" />; // unreachable
  const { color } = row.series;
  const y0 = ctx.baselineY(i);
  const bandPx = (BAND / RESID_DOMAIN) * HALF_SPAN;
  // Positive residual (over-predicted) sits ABOVE the baseline (smaller y).
  const yFor = (res: number): number => {
    const clamped = Math.max(-RESID_DOMAIN, Math.min(RESID_DOMAIN, res));
    return y0 - (clamped / RESID_DOMAIN) * HALF_SPAN;
  };
  return (
    <g key={i} data-role="resid-row" data-row-kind={row.kind}>
      <rect
        data-role="tolerance-band"
        x={ctx.padX} y={y0 - bandPx} width={ctx.plotW} height={bandPx * 2}
        fill="var(--color-ink)" opacity={0.05}
      />
      {row.series.teeth
        .filter((t) => t.observed && t.residual !== undefined)
        .map((t, j) => {
          const res = t.residual as number;
          const overflow = Math.abs(res) > RESID_DOMAIN;
          const hot = hoveredQ !== undefined && Math.abs(t.q - hoveredQ) <= TOL;
          const cx = ctx.x.to(t.q);
          const cy = yFor(res);
          const r = hot ? 4 : 2.6;
          return (
            <g
              key={j}
              data-role="resid-point"
              data-q={t.q}
              {...(overflow ? { "data-overflow": "true" } : {})}
              {...(hot ? { "data-hot": "true" } : {})}
              style={{ cursor: onHoverQ ? "pointer" : "default" }}
              {...(onHoverQ ? { onMouseEnter: () => onHoverQ(t.q), onMouseLeave: () => onHoverQ(undefined) } : {})}
            >
              <circle
                cx={cx} cy={cy} r={r}
                fill={overflow ? "none" : color}
                stroke={color}
                strokeWidth={overflow ? 1.6 : hot ? 1.4 : 0}
              />
              {hot ? <circle cx={cx} cy={cy} r={r + 3} fill="none" stroke={color} strokeWidth={1.4} opacity={0.6} /> : null}
            </g>
          );
        })}
    </g>
  );
}
```

- [ ] **Step 4: Write the story**

Create `src/print/comb/ResidualChart.stories.tsx`:

```tsx
import type { Meta, StoryObj } from "@storybook/react-vite";
import { ResidualChart } from "./ResidualChart";
import { PN3M, IM3M } from "./comb.fixtures";

const meta: Meta<typeof ResidualChart> = { title: "comb/ResidualChart", component: ResidualChart };
export default meta;
type Story = StoryObj<typeof ResidualChart>;

const frame = "rounded border border-hair bg-plate p-3";

export const Wide: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={760} /></div>,
};

export const NarrowScrolling: Story = {
  render: () => <div className={frame} style={{ width: 360 }}><ResidualChart assigned={[PN3M, IM3M]} maxWidth={360} /></div>,
};

// IM3M's √8 residual (0.041) exceeds the ±3% domain → an open overflow dot at the edge.
export const OverflowSignal: Story = {
  render: () => <div className={frame} style={{ width: 760 }}><ResidualChart assigned={[IM3M]} maxWidth={760} /></div>,
};
```

- [ ] **Step 5: Run tests + lint + typecheck**

```bash
npm test -- print-comb/ResidualChart && npm run lint:design && npx tsc --noEmit -p tsconfig.build.json
```
Expected: 8 ResidualChart tests PASS; `lint:design` exits 0; `tsc -p tsconfig.build.json` exits 0 (the CombChart story's `ResidualChart` import now resolves).

- [ ] **Step 6: Commit**

```bash
git add src/print/comb/ResidualChart.tsx src/print/comb/ResidualChart.stories.tsx test/print-comb/ResidualChart.test.tsx
git commit -m "comb: ResidualChart — per-row Δq/q points, tolerance band, R² gutter, overflow signal

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: Barrel + full build/storybook gate

**Files:**
- Create: `src/print/comb/index.ts`

- [ ] **Step 1: Create the barrel**

Create `src/print/comb/index.ts`:

```ts
export { CombChart } from "./CombChart";
export { ResidualChart } from "./ResidualChart";
export { CombScaffold } from "./CombScaffold";
export {
  assembleRows, combQDomain,
  type CombTooth, type CombSeries, type CombRow,
} from "./combModel";
```

- [ ] **Step 2: Run the whole comb suite**

```bash
npm test -- print-comb
```
Expected: all four files pass (5 + 5 + 8 + 8 = 26 tests).

- [ ] **Step 3: Full build + storybook gate**

```bash
npm run build && npm run build-storybook
```
Expected: both exit 0. `npm run build` = `lint:design && tsc -p tsconfig.build.json && vite build`. `build-storybook` compiles all `src/print/**/*.stories.tsx` including the two new comb stories.

- [ ] **Step 4: Commit**

```bash
git add src/print/comb/index.ts
git commit -m "comb: barrel export for CombChart / ResidualChart / CombScaffold / combModel

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Verification (whole batch)

- `npm test -- print-comb` — 26 tests green.
- `npm run lint:design` — exit 0 (proves `print/comb/` is correctly registered as exempt and nothing leaks appearance outside it).
- `npx tsc --noEmit -p tsconfig.build.json` — exit 0.
- `npm run build` — exit 0.
- `npm run build-storybook` — exit 0.
- Manual fidelity: `npm run storybook` → `comb/CombChart` (Wide / NarrowScrolling / WithPreviewRow / LinkedWithResidual) and `comb/ResidualChart` (Wide / NarrowScrolling / OverflowSignal). Confirm: teeth stay legible and the pane scrolls (gutter pinned) at 360px; the preview row baseline is dashed; hover lights the matching tooth/point in its own phase colour across both charts; the residual overflow dot reads as a clear open marker.

## Notes for the implementer

- **Do NOT push or merge.** This branch (`worktree-greenfield-ui-rebuild`) stays unmerged — the whole "The Print" rebuild lands as one piece.
- The comb's glyph vocabulary (tooth / caret / ring) is intentionally **its own**, distinct from the trace-plot `PeakGlyph` (triangle/diamond) and from the existing `CombLegend`. Do not wire in `PeakGlyph` or `CombLegend` — legend placement + CombLegend reconciliation is deferred to the future `CombsPanel` composite.
- Hover is **never** the terracotta accent — emphasis is wider stroke + a same-colour ring. This matches the detector rings and trace-plot engine.
- Fixtures live in `src/print/comb/comb.fixtures.ts` (a plain module, NOT a `*.stories.tsx` — CSF would treat data exports as stories). Both story files import `PN3M` / `IM3M` / `LEFTOVER` from there.
```
