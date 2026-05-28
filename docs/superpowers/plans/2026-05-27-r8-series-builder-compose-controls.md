# R8 — Series Builder Compose Controls Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the series-builder "compose the figure" controls (trace offset slider, log/linear scale toggle, floating offset dock in full-bleed), the figure-as-printed-plate framing (plate container, kicker tag-row, fig-sub, auto caption), and the autogroup read-rail card — matching `docs/redesign-mockups/series-builder.html`.

**Architecture:** The on-screen plot is `MultiTracePlot`, which already supports `xType: "log" | "linear"`. The waterfall vertical spacing is driven by `applyNormalization`'s `workingBandFraction` (default 0.7) — a larger fraction = taller traces filling more of each band = the "more offset" feel. We thread an `offset` (0.4–1.4) → `workingBandFraction` and a `scaleMode` ("log"/"linear") → `xType` from the page down through `MultiTracePlot` → `MemberTraceLayer` → `applyNormalization`. New rail UI (`OffsetSlider`, `ScaleToggle`, `AutogroupCard`, `OffsetDock`) lives in the rail / page. The plate container, kicker tag-row, fig-sub and caption are page-level JSX using existing Print tokens.

**Tech Stack:** React 18 + TypeScript strict, TailwindCSS 4 ("The Print" tokens), Vitest + RTL, Observable Plot, Playwright (live verification).

**Findings covered:** B-F (offset slider + scale toggle), B-G (floating offset dock), B-H (fig-tags + fig-sub + caption), B-I (autogroup card), B-J (figure-as-plate container), B-K/B-L (P3: edit-rail `bg-plate`, rail recessed margin).

**Out of scope (do NOT touch):** heatmap enablement + track-reflections (B-D/B-E → #208); leave the disabled heatmap button as-is. Export pipeline band-fill math (separate normalization; export already gated). `phases.ts` (R0b done). MultiTracePlot tooltip/selection chrome (R0c done).

---

## File structure

| File | Responsibility | Action |
|---|---|---|
| `src/lib/comparison/normalization.ts` | already accepts `workingBandFraction`; no change needed | read-only |
| `src/components/MemberTraceLayer.tsx` | thread `workingBandFraction` into `applyNormalization` | modify |
| `src/components/MultiTracePlot.tsx` | accept `workingBandFraction` prop, forward to layer | modify |
| `src/components/OffsetSlider.tsx` | rail offset slider (0.4–1.4×) | create |
| `src/components/ScaleToggle.tsx` | log/linear segmented toggle | create |
| `src/components/OffsetDock.tsx` | floating offset dock (full-bleed) | create |
| `src/components/AutogroupCard.tsx` | read-rail "Auto-grouped" card | create |
| `src/components/SeriesBuilderRail.tsx` | mount offset slider, scale toggle, autogroup card; `bg-plate` edit input; recessed margin | modify |
| `src/pages/SeriesBuilderPage.tsx` | offset/scaleMode state; plate container; kicker tags; fig-sub; caption; offset dock; thread offset/xType to plot | modify |
| `test/OffsetSlider.test.tsx` | slider behavior | create |
| `test/ScaleToggle.test.tsx` | toggle behavior | create |
| `test/OffsetDock.test.tsx` | dock behavior | create |
| `test/AutogroupCard.test.tsx` | autogroup card render | create |
| `test/SeriesBuilderPage.test.tsx` | plate, caption, offset/scale wiring to plot | modify |

**Offset → workingBandFraction mapping (load-bearing):** slider range 0.4–1.4 (mockup default 1.2). Map linearly to a workingBandFraction in `[0.45, 0.95]`:
`fraction = 0.45 + (offset - 0.4) / (1.4 - 0.4) * (0.95 - 0.45)`.
At offset 1.2 → fraction ≈ 0.85 (taller, denser stack, near mockup default); at 0.4 → 0.45 (loose); at 1.4 → 0.95 (dense). Implemented as a pure helper `offsetToBandFraction(offset)` in `MultiTracePlot.tsx` (exported for the test).

---

### Task 1: Thread `workingBandFraction` through MemberTraceLayer → applyNormalization

**Files:**
- Modify: `src/components/MemberTraceLayer.tsx`
- Test: `test/MultiTracePlot.test.tsx` (existing; we add a layer-level unit is overkill — covered by Task 2 plot prop test). This task is a refactor verified by `npm run build` + existing tests.

- [ ] **Step 1: Add `workingBandFraction` to MemberMarksProps and pass to applyNormalization**

In `src/components/MemberTraceLayer.tsx`, add to the `MemberMarksProps` interface (after `sampleIdFor?`):

```typescript
  /**
   * Fraction of each total band occupied by the working band (R8 offset
   * slider). Forwarded to `applyNormalization`; when omitted the library
   * default (DEFAULT_WORKING_BAND_FRACTION = 0.7) applies. A larger value
   * makes traces taller within their band — the "tighter waterfall" the
   * builder's offset slider composes.
   */
  workingBandFraction?: number;
```

In `buildMemberPeakRows`, destructure it and pass it through:

```typescript
  const { member, trace, yBand, peakDisplay, highlightedIndexId, workingBandFraction } = props;
```

and change the `applyNormalization` call:

```typescript
  const linePoints = applyNormalization(
    { q: trace.q, I: trace.I },
    reference,
    yBand,
    workingBandFraction,
  );
```

(`applyNormalization`'s 4th param defaults to `DEFAULT_WORKING_BAND_FRACTION` when `undefined` — passing `undefined` is safe.)

- [ ] **Step 2: Run build to verify types**

Run: `npm run build`
Expected: PASS (tsc --noEmit + vite build clean).

- [ ] **Step 3: Run existing layer/plot tests**

Run: `npm test -- MultiTracePlot`
Expected: PASS (no behavior change — default fraction preserved).

- [ ] **Step 4: Commit**

```bash
git add src/components/MemberTraceLayer.tsx
git commit -m "R8: thread workingBandFraction through MemberTraceLayer (#231)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: MultiTracePlot accepts `workingBandFraction` + `offsetToBandFraction` helper

**Files:**
- Modify: `src/components/MultiTracePlot.tsx`
- Test: `test/MultiTracePlot.test.tsx`

- [ ] **Step 1: Write the failing test**

Add to `test/MultiTracePlot.test.tsx` (inside the top-level describe; import `offsetToBandFraction` from the component):

```typescript
import { offsetToBandFraction } from "../src/components/MultiTracePlot";

describe("offsetToBandFraction", () => {
  it("maps the slider range 0.4..1.4 into a band fraction 0.45..0.95", () => {
    expect(offsetToBandFraction(0.4)).toBeCloseTo(0.45, 5);
    expect(offsetToBandFraction(1.4)).toBeCloseTo(0.95, 5);
  });
  it("maps the mockup default 1.2 to ~0.85", () => {
    expect(offsetToBandFraction(1.2)).toBeCloseTo(0.85, 2);
  });
  it("clamps out-of-range offsets to the fraction bounds", () => {
    expect(offsetToBandFraction(0)).toBeCloseTo(0.45, 5);
    expect(offsetToBandFraction(99)).toBeCloseTo(0.95, 5);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- MultiTracePlot`
Expected: FAIL — `offsetToBandFraction is not a function` / no export.

- [ ] **Step 3: Add the helper + prop + forwarding**

In `src/components/MultiTracePlot.tsx`, add near the top (after the margin consts, before `MultiTracePlotProps`):

```typescript
/**
 * Map the builder's trace-offset slider (0.4..1.4, mockup default 1.2) to the
 * working-band fraction consumed by `applyNormalization`. A larger offset =
 * taller traces filling more of each band = a tighter waterfall. Clamped so
 * the slider can never collapse a band (min 0.45) or overflow it (max 0.95).
 */
export function offsetToBandFraction(offset: number): number {
  const t = Math.min(1, Math.max(0, (offset - 0.4) / (1.4 - 0.4)));
  return 0.45 + t * (0.95 - 0.45);
}
```

Add to `MultiTracePlotProps` (after `sampleIdFor?`):

```typescript
  /**
   * Working-band fraction for the waterfall offset slider (R8, #231). When
   * omitted the layer falls back to `DEFAULT_WORKING_BAND_FRACTION`. The page
   * derives this from its offset slider via `offsetToBandFraction`.
   */
  workingBandFraction?: number;
```

In the `MultiTracePlot` props destructure, add `workingBandFraction,`. Then in `renderPlot`, inside the `layerProps` object literal (where `showPeakTicks` etc. are set), add:

```typescript
        ...(workingBandFraction !== undefined ? { workingBandFraction } : {}),
```

Finally, add `workingBandFraction` to the `renderPlot` `useCallback` dependency array (find the dep array at the end of the `useCallback(() => {...}, [...])` and add it).

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- MultiTracePlot`
Expected: PASS.

- [ ] **Step 5: Run build**

Run: `npm run build`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add src/components/MultiTracePlot.tsx test/MultiTracePlot.test.tsx
git commit -m "R8: MultiTracePlot accepts workingBandFraction + offsetToBandFraction (#231)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: OffsetSlider component

**Files:**
- Create: `src/components/OffsetSlider.tsx`
- Test: `test/OffsetSlider.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/OffsetSlider.test.tsx`:

```typescript
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { OffsetSlider } from "../src/components/OffsetSlider";

describe("OffsetSlider", () => {
  it("renders the current value to two decimals", () => {
    render(<OffsetSlider value={1.2} onChange={vi.fn()} />);
    expect(screen.getByTestId("offset-slider-value")).toHaveTextContent("1.20");
  });

  it("calls onChange with the parsed numeric value on input", () => {
    const onChange = vi.fn();
    render(<OffsetSlider value={1.2} onChange={onChange} />);
    fireEvent.change(screen.getByTestId("offset-slider"), { target: { value: "0.8" } });
    expect(onChange).toHaveBeenCalledWith(0.8);
  });

  it("exposes the slider range and step matching the mockup", () => {
    render(<OffsetSlider value={1.2} onChange={vi.fn()} />);
    const input = screen.getByTestId("offset-slider") as HTMLInputElement;
    expect(input.min).toBe("0.4");
    expect(input.max).toBe("1.4");
    expect(input.step).toBe("0.05");
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- OffsetSlider`
Expected: FAIL — cannot find module `OffsetSlider`.

- [ ] **Step 3: Implement OffsetSlider**

Create `src/components/OffsetSlider.tsx`:

```typescript
interface OffsetSliderProps {
  /** Current offset multiplier (0.4..1.4). */
  value: number;
  onChange: (next: number) => void;
}

/**
 * OffsetSlider — the series-builder "trace offset" control (R8 / B-F). Scales
 * the vertical spacing of the waterfall stack. The page maps this 0.4..1.4
 * range to a working-band fraction via `offsetToBandFraction`. Mockup:
 * `series-builder.html` `.slider-row` / `#offset`.
 */
export function OffsetSlider({ value, onChange }: OffsetSliderProps): JSX.Element {
  return (
    <div className="flex flex-col gap-1.5" data-testid="offset-slider-row">
      <div className="flex items-baseline justify-between">
        <span className="text-xs font-semibold text-ink">Trace offset</span>
        <span className="text-xs font-semibold text-ink" data-testid="offset-slider-value">
          {value.toFixed(2)}&times;
        </span>
      </div>
      <input
        type="range"
        data-testid="offset-slider"
        aria-label="Trace offset"
        min="0.4"
        max="1.4"
        step="0.05"
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        className="h-[3px] w-full cursor-pointer appearance-none rounded-full bg-hair-strong accent-print-accent"
      />
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- OffsetSlider`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/components/OffsetSlider.tsx test/OffsetSlider.test.tsx
git commit -m "R8: OffsetSlider compose control (#231)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: ScaleToggle component

**Files:**
- Create: `src/components/ScaleToggle.tsx`
- Test: `test/ScaleToggle.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/ScaleToggle.test.tsx`:

```typescript
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScaleToggle } from "../src/components/ScaleToggle";

describe("ScaleToggle", () => {
  it("marks the active scale via aria-pressed", () => {
    render(<ScaleToggle value="log" onChange={vi.fn()} />);
    expect(screen.getByTestId("scale-log")).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByTestId("scale-linear")).toHaveAttribute("aria-pressed", "false");
  });

  it("calls onChange('linear') when linear is clicked", () => {
    const onChange = vi.fn();
    render(<ScaleToggle value="log" onChange={onChange} />);
    fireEvent.click(screen.getByTestId("scale-linear"));
    expect(onChange).toHaveBeenCalledWith("linear");
  });

  it("calls onChange('log') when log is clicked", () => {
    const onChange = vi.fn();
    render(<ScaleToggle value="linear" onChange={onChange} />);
    fireEvent.click(screen.getByTestId("scale-log"));
    expect(onChange).toHaveBeenCalledWith("log");
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- ScaleToggle`
Expected: FAIL — cannot find module `ScaleToggle`.

- [ ] **Step 3: Implement ScaleToggle**

Create `src/components/ScaleToggle.tsx`:

```typescript
export type ScaleMode = "log" | "linear";

interface ScaleToggleProps {
  value: ScaleMode;
  onChange: (next: ScaleMode) => void;
}

/**
 * ScaleToggle — log/linear q-axis segmented control (R8 / B-F). Drives
 * MultiTracePlot's `xType`. Active segment uses the ink-fill the mockup's
 * `.seg button.on` defines (`background:var(--ink);color:var(--paper)`),
 * matching the sibling RepresentationToggle's active state.
 */
export function ScaleToggle({ value, onChange }: ScaleToggleProps): JSX.Element {
  return (
    <div
      data-testid="scale-toggle"
      role="group"
      aria-label="q-axis scale"
      className="inline-flex overflow-hidden rounded border border-hair-strong"
    >
      <button
        type="button"
        data-testid="scale-log"
        aria-pressed={value === "log"}
        onClick={() => onChange("log")}
        className={`flex-1 px-3 py-1.5 text-xs ${value === "log" ? "bg-ink text-paper" : "text-ink-faint"}`}
      >
        log q
      </button>
      <button
        type="button"
        data-testid="scale-linear"
        aria-pressed={value === "linear"}
        onClick={() => onChange("linear")}
        className={`flex-1 px-3 py-1.5 text-xs ${value === "linear" ? "bg-ink text-paper" : "text-ink-faint"}`}
      >
        linear q
      </button>
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- ScaleToggle`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/components/ScaleToggle.tsx test/ScaleToggle.test.tsx
git commit -m "R8: ScaleToggle log/linear control (#231)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 5: AutogroupCard component

**Files:**
- Create: `src/components/AutogroupCard.tsx`
- Test: `test/AutogroupCard.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/AutogroupCard.test.tsx`:

```typescript
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { AutogroupCard } from "../src/components/AutogroupCard";

describe("AutogroupCard", () => {
  it("reads N samples and the ordering variable into the body", () => {
    render(
      <AutogroupCard sampleCount={6} orderingVariable="LL37 : lipid ratio" onConfirm={vi.fn()} onAdjust={vi.fn()} />,
    );
    expect(screen.getByTestId("autogroup-card")).toHaveTextContent("6 samples");
    expect(screen.getByTestId("autogroup-card")).toHaveTextContent("LL37 : lipid ratio");
  });

  it("falls back to a placeholder when ordering variable is null", () => {
    render(<AutogroupCard sampleCount={3} orderingVariable={null} onConfirm={vi.fn()} onAdjust={vi.fn()} />);
    expect(screen.getByTestId("autogroup-card")).toHaveTextContent("their names");
  });

  it("fires onConfirm and onAdjust", () => {
    const onConfirm = vi.fn();
    const onAdjust = vi.fn();
    render(<AutogroupCard sampleCount={6} orderingVariable="dose" onConfirm={onConfirm} onAdjust={onAdjust} />);
    fireEvent.click(screen.getByTestId("autogroup-confirm"));
    fireEvent.click(screen.getByTestId("autogroup-adjust"));
    expect(onConfirm).toHaveBeenCalled();
    expect(onAdjust).toHaveBeenCalled();
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- AutogroupCard`
Expected: FAIL — cannot find module `AutogroupCard`.

- [ ] **Step 3: Implement AutogroupCard**

Create `src/components/AutogroupCard.tsx`:

```typescript
interface AutogroupCardProps {
  sampleCount: number;
  orderingVariable: string | null;
  onConfirm: () => void;
  onAdjust: () => void;
}

/**
 * AutogroupCard — the confident-expert "Auto-grouped" note in the series
 * builder read rail (R8 / B-I). Mockup `series-builder.html` `.autogroup`:
 * a spark glyph + "Auto-grouped" title, a sentence naming how many samples
 * Himalaya read as one series and the ordering variable, and Confirm / Adjust
 * link buttons. Plate surface, hair border.
 */
export function AutogroupCard({
  sampleCount, orderingVariable, onConfirm, onAdjust,
}: AutogroupCardProps): JSX.Element {
  return (
    <div
      data-testid="autogroup-card"
      className="rounded-md border border-hair bg-plate p-3.5"
    >
      <div className="mb-1.5 flex items-center gap-1.5">
        <svg className="h-[15px] w-[15px] shrink-0" viewBox="0 0 16 16" fill="none" aria-hidden="true">
          <path
            d="M8 1.4l1.7 4.2 4.5.3-3.5 2.9 1.2 4.4L8 10.9 4.1 13.2l1.2-4.4L1.8 5.9l4.5-.3z"
            fill="var(--color-print-accent)"
          />
        </svg>
        <span className="text-xs font-bold text-ink">Auto-grouped</span>
      </div>
      <p className="text-xs leading-relaxed text-ink-soft">
        Himalaya read <b className="font-semibold text-ink">{sampleCount} samples</b> as one
        series from{" "}
        {orderingVariable !== null ? (
          <>
            their names, ordered by{" "}
            <b className="font-semibold text-ink">{orderingVariable}</b>.
          </>
        ) : (
          <>their names.</>
        )}
      </p>
      <div className="mt-2 flex gap-3">
        <button
          type="button"
          data-testid="autogroup-confirm"
          onClick={onConfirm}
          className="text-xs font-semibold text-print-accent hover:underline"
        >
          Confirm series
        </button>
        <button
          type="button"
          data-testid="autogroup-adjust"
          onClick={onAdjust}
          className="text-xs font-semibold text-ink-faint hover:underline"
        >
          Adjust
        </button>
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- AutogroupCard`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/components/AutogroupCard.tsx test/AutogroupCard.test.tsx
git commit -m "R8: AutogroupCard read-rail note (#231)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 6: OffsetDock component (floating, full-bleed)

**Files:**
- Create: `src/components/OffsetDock.tsx`
- Test: `test/OffsetDock.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/OffsetDock.test.tsx`:

```typescript
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { OffsetDock } from "../src/components/OffsetDock";

describe("OffsetDock", () => {
  it("does not render when hidden", () => {
    render(<OffsetDock show={false} value={1.2} onChange={vi.fn()} />);
    expect(screen.queryByTestId("offset-dock")).not.toBeInTheDocument();
  });

  it("renders the dock + value when shown", () => {
    render(<OffsetDock show value={1.2} onChange={vi.fn()} />);
    expect(screen.getByTestId("offset-dock")).toBeInTheDocument();
    expect(screen.getByTestId("offset-dock-value")).toHaveTextContent("1.20");
  });

  it("calls onChange with the parsed value on input", () => {
    const onChange = vi.fn();
    render(<OffsetDock show value={1.2} onChange={onChange} />);
    fireEvent.change(screen.getByTestId("offset-dock-slider"), { target: { value: "0.6" } });
    expect(onChange).toHaveBeenCalledWith(0.6);
  });
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- OffsetDock`
Expected: FAIL — cannot find module `OffsetDock`.

- [ ] **Step 3: Implement OffsetDock**

Create `src/components/OffsetDock.tsx`:

```typescript
interface OffsetDockProps {
  /** Only renders when true — shown only in full-bleed (rail collapsed). */
  show: boolean;
  value: number;
  onChange: (next: number) => void;
}

/**
 * OffsetDock — floating offset control kept reachable when the rail is
 * collapsed into full-bleed (R8 / B-G). Mockup `series-builder.html` `.dock`:
 * a plate-surfaced pill pinned bottom-right with an "Offset" label, the
 * slider, and the live value. Mirrors the rail OffsetSlider's range/step so
 * the page can keep both in sync.
 */
export function OffsetDock({ show, value, onChange }: OffsetDockProps): JSX.Element | null {
  if (!show) return null;
  return (
    <div
      data-testid="offset-dock"
      className="fixed bottom-6 right-6 z-10 flex items-center gap-3.5 rounded-xl border border-hair-strong bg-plate px-4 py-3 shadow-[0_8px_26px_-10px_rgba(60,52,40,.34)]"
    >
      <span className="text-[10px] font-bold uppercase tracking-wide text-ink-faint">Offset</span>
      <input
        type="range"
        data-testid="offset-dock-slider"
        aria-label="Trace offset"
        min="0.4"
        max="1.4"
        step="0.05"
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        className="w-32 cursor-pointer appearance-none rounded-full bg-hair-strong accent-print-accent"
      />
      <span
        data-testid="offset-dock-value"
        className="min-w-[40px] text-xs font-bold text-ink"
      >
        {value.toFixed(2)}&times;
      </span>
    </div>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- OffsetDock`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/components/OffsetDock.tsx test/OffsetDock.test.tsx
git commit -m "R8: OffsetDock floating full-bleed control (#231)

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 7: Wire offset/scale + autogroup into SeriesBuilderRail

**Files:**
- Modify: `src/components/SeriesBuilderRail.tsx`
- Test: `test/SeriesBuilderRail.test.tsx`

- [ ] **Step 1: Write the failing tests**

Add to `test/SeriesBuilderRail.test.tsx` — first extend the `setup` defaults so the new required props are present, then add cases. Replace the `setup` function's `props` object with:

```typescript
  const props = {
    collapsed: false,
    onToggleCollapsed: vi.fn(),
    representation: "waterfall" as const,
    onRepresentationChange: vi.fn(),
    orderingVariable: "LL37 : lipid ratio",
    offset: 1.2,
    onOffsetChange: vi.fn(),
    scaleMode: "log" as const,
    onScaleModeChange: vi.fn(),
    sampleCount: 6,
    onConfirmSeries: vi.fn(),
    onAdjustSeries: vi.fn(),
    exportControls: <div data-testid="mock-export" />,
    ...over,
  };
```

Add cases inside the `describe`:

```typescript
  it("renders the offset slider, scale toggle, and autogroup card", () => {
    setup();
    expect(screen.getByTestId("offset-slider")).toBeInTheDocument();
    expect(screen.getByTestId("scale-toggle")).toBeInTheDocument();
    expect(screen.getByTestId("autogroup-card")).toBeInTheDocument();
  });

  it("forwards offset slider input to onOffsetChange", () => {
    const props = setup();
    fireEvent.change(screen.getByTestId("offset-slider"), { target: { value: "0.8" } });
    expect(props.onOffsetChange).toHaveBeenCalledWith(0.8);
  });

  it("forwards scale toggle to onScaleModeChange", () => {
    const props = setup();
    fireEvent.click(screen.getByTestId("scale-linear"));
    expect(props.onScaleModeChange).toHaveBeenCalledWith("linear");
  });

  it("hides the autogroup card and compose controls in edit mode", () => {
    setup({ editControls: <div data-testid="mock-edit" /> });
    expect(screen.queryByTestId("autogroup-card")).not.toBeInTheDocument();
    expect(screen.getByTestId("mock-edit")).toBeInTheDocument();
  });
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- SeriesBuilderRail`
Expected: FAIL — offset-slider / scale-toggle / autogroup-card not found.

- [ ] **Step 3: Implement the rail changes**

Rewrite `src/components/SeriesBuilderRail.tsx`:

```typescript
import type { ReactNode } from "react";
import { RepresentationToggle, type Representation } from "./RepresentationToggle";
import { OffsetSlider } from "./OffsetSlider";
import { ScaleToggle, type ScaleMode } from "./ScaleToggle";
import { AutogroupCard } from "./AutogroupCard";

interface SeriesBuilderRailProps {
  collapsed: boolean;
  onToggleCollapsed: () => void;
  representation: Representation;
  onRepresentationChange: (next: Representation) => void;
  orderingVariable: string | null;
  /** Trace-offset slider (B-F). 0.4..1.4. */
  offset: number;
  onOffsetChange: (next: number) => void;
  /** log/linear q-axis scale (B-F). */
  scaleMode: ScaleMode;
  onScaleModeChange: (next: ScaleMode) => void;
  /** Autogroup card (B-I): how many samples Himalaya read as this series. */
  sampleCount: number;
  onConfirmSeries: () => void;
  onAdjustSeries: () => void;
  /** Parent injects <FigureExportControls/> (it owns the export spec thunk). */
  exportControls: ReactNode;
  /**
   * I3.5b — recipe-edit controls. Injected only in edit mode; when present
   * the static ordering-variable line + the autogroup card are replaced by
   * the recipe editor (the compose controls stay, since offset/scale shape
   * the figure in both modes).
   */
  editControls?: ReactNode;
}

/**
 * SeriesBuilderRail — the editing margin of the series builder. Read mode
 * shows the autogroup card (B-I), the static ordering-variable line, the
 * representation toggle, the compose controls (offset slider + scale toggle,
 * B-F), and the injected export controls. Edit mode (I3.5b) swaps the
 * autogroup + ordering line for the recipe editor. The rail collapses to
 * full-bleed; the floating OffsetDock (rendered by the page) keeps offset
 * reachable while collapsed (B-G).
 *
 * B-K/B-L: edit inputs use `bg-plate`; the rail itself is the recessed
 * `bg-paper-sunk` margin against the bright plate (matching the mockup
 * `.rail{background:var(--paper-sunk)}`).
 */
export function SeriesBuilderRail({
  collapsed, onToggleCollapsed, representation, onRepresentationChange,
  orderingVariable, offset, onOffsetChange, scaleMode, onScaleModeChange,
  sampleCount, onConfirmSeries, onAdjustSeries, exportControls, editControls,
}: SeriesBuilderRailProps): JSX.Element {
  if (collapsed) {
    return (
      <button
        type="button"
        data-testid="rail-restore"
        onClick={onToggleCollapsed}
        title="Show the editing rail"
        className="flex items-center gap-1.5 border-l border-hair px-2 text-xs text-ink-faint"
      >
        <span aria-hidden="true">‹</span> Compose
      </button>
    );
  }
  return (
    <aside
      data-testid="series-builder-rail"
      className="flex w-[336px] shrink-0 flex-col gap-5 overflow-y-auto border-l border-hair bg-paper-sunk p-4"
    >
      <div className="flex items-center justify-between">
        <span className="text-xs font-semibold uppercase tracking-wide text-ink">Compose</span>
        <button
          type="button"
          data-testid="rail-collapse-toggle"
          onClick={onToggleCollapsed}
          title="Collapse the rail — full-bleed"
          className="rounded px-1.5 text-ink-faint hover:text-ink"
        >
          ›
        </button>
      </div>

      {editControls !== undefined ? (
        <section className="flex flex-col gap-1.5 [&_input]:bg-plate" data-testid="rail-edit">
          {editControls}
        </section>
      ) : (
        <>
          <AutogroupCard
            sampleCount={sampleCount}
            orderingVariable={orderingVariable}
            onConfirm={onConfirmSeries}
            onAdjust={onAdjustSeries}
          />
          <section className="flex flex-col gap-1">
            <div className="text-xs font-semibold text-ink-faint">Ordering variable</div>
            <div className="text-sm text-ink">{orderingVariable ?? "—"}</div>
          </section>
        </>
      )}

      <section className="flex flex-col gap-1.5">
        <div className="text-xs font-semibold text-ink-faint">Representation</div>
        <RepresentationToggle value={representation} onChange={onRepresentationChange} />
      </section>

      <section className="flex flex-col gap-2.5" data-testid="rail-display">
        <div className="text-xs font-semibold text-ink-faint">Display</div>
        <OffsetSlider value={offset} onChange={onOffsetChange} />
        <ScaleToggle value={scaleMode} onChange={onScaleModeChange} />
      </section>

      <section className="flex flex-col gap-1.5" data-testid="rail-export">
        {exportControls}
      </section>
    </aside>
  );
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- SeriesBuilderRail`
Expected: PASS.

- [ ] **Step 5: Run build**

Run: `npm run build`
Expected: FAIL — `SeriesBuilderPage.tsx` doesn't yet pass the new required rail props. (Expected; Task 8 fixes it. Do NOT commit yet — proceed to Task 8 so the build is green, then commit both. If you must commit per-task, note the build is intentionally broken between Task 7 and Task 8.)

Skip the commit here; the rail + page land together at the end of Task 8.

---

### Task 8: Wire SeriesBuilderPage — state, plate, kicker tags, caption, dock

**Files:**
- Modify: `src/pages/SeriesBuilderPage.tsx`
- Test: `test/SeriesBuilderPage.test.tsx`

- [ ] **Step 1: Write the failing tests**

The existing `SeriesBuilderPage.test.tsx` mocks `MultiTracePlot` and records props into `h.lastPlotProps`. Extend that mock's recorded type and add cases. Change the `lastPlotProps` hoist type to include the new props:

```typescript
  lastPlotProps: undefined as undefined | {
    showPeakTicks?: boolean; showPeakLabels?: boolean;
    xType?: "log" | "linear"; workingBandFraction?: number;
  },
```

and the mock component signature/attrs:

```typescript
  MultiTracePlot: (props: {
    showPeakTicks?: boolean; showPeakLabels?: boolean;
    xType?: "log" | "linear"; workingBandFraction?: number;
  }) => {
    h.lastPlotProps = props;
    return (
      <div
        data-testid="mock-multi-trace-plot"
        data-show-peak-ticks={String(props.showPeakTicks)}
        data-show-peak-labels={String(props.showPeakLabels)}
        data-x-type={String(props.xType)}
      />
    );
  },
```

Add cases inside the describe:

```typescript
  it("renders the figure-as-plate container with kicker tags and a caption", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("series-builder-plate")).toBeInTheDocument();
    expect(screen.getByTestId("fig-tags")).toBeInTheDocument();
    expect(screen.getByTestId("fig-caption")).toBeInTheDocument();
  });

  it("forwards the default scale (log) and offset to the plot", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("mock-multi-trace-plot")).toHaveAttribute("data-x-type", "log");
    expect(h.lastPlotProps?.workingBandFraction).toBeCloseTo(0.85, 2);
  });

  it("flips the plot to linear when the scale toggle is switched", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    fireEvent.click(screen.getByTestId("scale-linear"));
    expect(h.lastPlotProps?.xType).toBe("linear");
    // the kicker scale tag reflects the new scale
    expect(screen.getByTestId("fig-tags")).toHaveTextContent("linear q");
  });

  it("changes the plot offset (workingBandFraction) when the slider moves", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    fireEvent.change(screen.getByTestId("offset-slider"), { target: { value: "0.4" } });
    expect(h.lastPlotProps?.workingBandFraction).toBeCloseTo(0.45, 5);
  });

  it("shows the floating offset dock only when the rail is collapsed", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.queryByTestId("offset-dock")).not.toBeInTheDocument();
    fireEvent.click(screen.getByTestId("rail-collapse-toggle"));
    expect(screen.getByTestId("offset-dock")).toBeInTheDocument();
  });
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- SeriesBuilderPage.test`
Expected: FAIL — series-builder-plate / fig-tags / scale-linear not found; workingBandFraction undefined.

- [ ] **Step 3: Implement the page changes**

In `src/pages/SeriesBuilderPage.tsx`:

(a) Add imports at the top:

```typescript
import { OffsetDock } from "../components/OffsetDock";
import { offsetToBandFraction } from "../components/MultiTracePlot";
import type { ScaleMode } from "../components/ScaleToggle";
```

(b) In `SeriesBuilderBody`, after the `const [representation, setRepresentation] = ...` line, add compose state:

```typescript
  // Compose controls (R8 / B-F): the trace-offset slider and the log/linear
  // q-axis scale. Local UI state — a read surface's composition is a local
  // concern (persisting it is a recipe edit, deferred). Offset maps to the
  // plot's working-band fraction; scaleMode maps to xType.
  const [offset, setOffset] = useState(1.2);
  const [scaleMode, setScaleMode] = useState<ScaleMode>("log");
  const workingBandFraction = offsetToBandFraction(offset);
```

(c) Pass the two new props to `MultiTracePlot` (inside the `<MultiTracePlot ... />`):

```typescript
                    xType={scaleMode}
                    workingBandFraction={workingBandFraction}
```

(d) Replace the `SeriesBuilderBody` return's outer plot column markup. Find the block:

```typescript
        <div className="flex-1 min-h-0 flex flex-col p-4 gap-3" data-testid="series-builder-plot">
          {members.length === 0 ? (
```

Wrap the figure content in a plate. Replace that opening `<div ...series-builder-plot>` and the empty/loaded branch wrapper so the loaded branch renders the plate. Use this structure for the whole plot column (keep the empty branch as-is inside the plate-less path):

```typescript
        <div
          className="flex-1 min-h-0 overflow-auto flex flex-col items-center px-8 py-6"
          data-testid="series-builder-plot"
        >
          {members.length === 0 ? (
            <div
              data-testid="series-builder-empty"
              className="flex-1 grid place-items-center text-sm text-ink-faint"
            >
              This series has no members yet.
            </div>
          ) : (
            <div
              data-testid="series-builder-plate"
              className={`w-full ${collapsed ? "max-w-[1336px]" : "max-w-[1180px]"} rounded border border-hair bg-plate p-8 shadow-[0_1px_1px_rgba(60,52,40,.04),0_18px_40px_-20px_rgba(60,52,40,.22)] transition-[max-width] duration-200`}
            >
              <div className="mb-2 flex items-baseline gap-3" data-testid="fig-tags">
                <span className="text-[11px] font-bold uppercase tracking-[0.14em] text-print-accent">
                  Series
                </span>
                <div className="flex gap-1.5">
                  <span className="rounded-full border border-hair px-2 py-px text-[10.5px] text-ink-faint">
                    {members.length} {members.length === 1 ? "sample" : "samples"}
                  </span>
                  <span className="rounded-full border border-hair px-2 py-px text-[10.5px] text-ink-faint">
                    {scaleMode === "log" ? "log q" : "linear q"}
                  </span>
                  <span className="rounded-full border border-hair px-2 py-px text-[10.5px] text-ink-faint">
                    offset waterfall
                  </span>
                </div>
              </div>
              <h1 className="text-display font-medium text-ink">{s.title || "Untitled series"}</h1>
              {s.description && (
                <p className="mt-2 max-w-[64ch] text-sm text-ink-soft" data-testid="fig-sub">
                  {s.description}
                </p>
              )}

              <div className="mt-4 flex items-center gap-3" data-testid="series-builder-controls">
                <GroupingModeToggle mode={groupingMode} onChange={setGroupingMode} />
                <AnnotationToggles />
              </div>
              <div className="mt-3 flex min-h-0 flex-row gap-2" style={{ height: "60vh" }}>
                <div ref={plotColRef} className="min-w-0 flex-1">
                  <MultiTracePlot
                    members={members}
                    traces={traces}
                    xDomain={xDomain}
                    onXDomain={setXDomain}
                    groupingMode={groupingMode}
                    sampleIdFor={sampleIdFor}
                    showPeakTicks={showPeakTicks}
                    showPeakLabels={showPeakLabels}
                    xType={scaleMode}
                    workingBandFraction={workingBandFraction}
                  />
                </div>
                <div className="w-[280px] shrink-0" data-testid="series-builder-gutter">
                  <MemberMetaGutter
                    members={members}
                    panelHeight={panelHeight}
                    mode="review"
                    displayLabelByMemberId={displayLabelByMemberId}
                  />
                </div>
              </div>

              <div
                data-testid="fig-caption"
                className="mt-3 flex gap-2 border-t border-hair pt-3 text-xs leading-relaxed text-ink-soft"
              >
                <span className="font-bold text-ink">Fig.</span>
                <span>
                  {members.length} 1D integration{members.length === 1 ? "" : "s"}, vertically
                  offset by {offset.toFixed(2)}× the band height
                  {s.ordering_variable ? <>, ordered by {s.ordering_variable}</> : null}. Peak ticks
                  coloured by indexed phase.
                </span>
              </div>
            </div>
          )}
        </div>
```

Note: this removes the old `gap-3`/`p-4` wrapper. The `ActiveBandProvider` and the outer `<div ...series-builder-body>` flex row stay unchanged.

(e) Update the `<SeriesBuilderRail .../>` invocation to pass the new props, and render the dock after it. Replace the `<SeriesBuilderRail ... />` element with:

```typescript
        <SeriesBuilderRail
          collapsed={collapsed}
          onToggleCollapsed={() => setCollapsed((c) => !c)}
          representation={representation}
          onRepresentationChange={setRepresentation}
          orderingVariable={s.ordering_variable}
          offset={offset}
          onOffsetChange={setOffset}
          scaleMode={scaleMode}
          onScaleModeChange={setScaleMode}
          sampleCount={s.samples.length || members.length}
          onConfirmSeries={() => { /* read surface: confirm is a no-op visual affordance (recipe edit deferred) */ }}
          onAdjustSeries={() => startSeriesDraftFromRail()}
          exportControls={
            <FigureExportControls
              spec={exportSpec}
              filenameStem={exportFilenameStem}
              ariaContext="series figure"
              disabled={exportDisabled}
            />
          }
          {...(editing
            ? { editControls: <SeriesRecipeEditor seriesId={s.id} members={members} /> }
            : {})}
        />
        <OffsetDock
          show={collapsed && representation === "waterfall"}
          value={offset}
          onChange={setOffset}
        />
```

(f) The rail's `onAdjustSeries` wants to enter edit mode. `SeriesBuilderBody` doesn't currently have the draft-start action. Add it: the top-level page (`SeriesBuilderPage`) already reads `startSeriesDraftFromSeries`. Pass an `onAdjust` down. Simplest: in `SeriesBuilderBody`, read the action from Zustand directly. Add near the other `useAppState` reads in `SeriesBuilderBody`:

```typescript
  const startSeriesDraftFromRail = useAppState((st) => () => st.startSeriesDraftFromSeries(s));
```

Wait — that selector form recreates a function each render. Instead add a plain selector + local callback. Replace the line above with:

```typescript
  const startDraft = useAppState((st) => st.startSeriesDraftFromSeries);
```

and define the handler inline in the JSX: `onAdjustSeries={() => startDraft(s)}`. Update the `<SeriesBuilderRail>` prop accordingly (replace `onAdjustSeries={() => startSeriesDraftFromRail()}` with `onAdjustSeries={() => startDraft(s)}`).

- [ ] **Step 4: Run the page tests**

Run: `npm test -- SeriesBuilderPage.test`
Expected: PASS.

- [ ] **Step 5: Run the edit-mode page test (regression)**

Run: `npm test -- SeriesBuilderPage.edit`
Expected: PASS (edit mode still injects the recipe editor; autogroup hidden).

- [ ] **Step 6: Run build**

Run: `npm run build`
Expected: PASS.

- [ ] **Step 7: Run the full unit suite**

Run: `npm test`
Expected: PASS (capture to a file if slow; grep for failures).

- [ ] **Step 8: Commit (rail + page together)**

```bash
git add src/components/SeriesBuilderRail.tsx test/SeriesBuilderRail.test.tsx \
        src/pages/SeriesBuilderPage.tsx test/SeriesBuilderPage.test.tsx
git commit -m "R8: wire compose controls + figure-as-plate into the builder (#231)

Offset slider + log/linear scale toggle (B-F), floating offset dock in
full-bleed (B-G), kicker tag-row + fig-sub + auto caption (B-H), autogroup
read-rail card (B-I), figure-as-plate container (B-J), edit-input bg-plate +
recessed rail margin (B-K/B-L).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>"
```

---

### Task 9: Live verification (screenshots)

**Files:** none (verification only).

- [ ] **Step 1: Start the dev server against the shared backend**

```bash
cd .claude/worktrees/r8-series-builder/packages/HimalayaUI/frontend
VITE_API_PORT=8091 npm run dev -- --host 127.0.0.1 --port 5195 > /tmp/r8-vite.log 2>&1 &
```

(If 5195 busy, use 5200.)

- [ ] **Step 2: Screenshot `/series/1` and compare to the mockup**

Use Playwright MCP `browser_navigate` to `http://127.0.0.1:5195/series/1`, `browser_take_screenshot` → `/tmp/r8-builder.png`. Confirm: plate container (white plate, hair border, shadow, centered), kicker "Series" + tags, serif title, caption block, autogroup card in the rail, offset slider, log/linear toggle.

- [ ] **Step 3: Exercise the offset slider + scale toggle**

Screenshot before; drag the offset slider (`browser_drag` or `browser_evaluate` set value + dispatch input); screenshot after (`/tmp/r8-offset-after.png`) — the stack spacing must visibly change. Click `scale-linear`; screenshot (`/tmp/r8-linear.png`) — the q-axis must switch to linear and the kicker tag read "linear q".

- [ ] **Step 4: Verify the floating dock in full-bleed**

Click `rail-collapse-toggle`; screenshot (`/tmp/r8-fullbleed.png`) — the rail folds, the plate widens to max-w-1336, and the floating offset dock appears bottom-right.

- [ ] **Step 5: Tear down**

```bash
lsof -ti:5195 | xargs -r kill
```

Leave :8091 running. Never bind 5173/8080.

---

## Self-Review notes

- **Spec coverage:** B-F (Tasks 3,4,8), B-G (Tasks 6,8), B-H (Task 8 plate kicker/sub/caption), B-I (Tasks 5,7), B-J (Task 8 plate), B-K/B-L (Task 7 `bg-plate` edit input + `bg-paper-sunk` rail). Out-of-scope heatmap/track left untouched.
- **Type consistency:** `ScaleMode` defined once in `ScaleToggle.tsx`, imported by rail + page. `offsetToBandFraction` defined+exported in `MultiTracePlot.tsx`, imported by page + tested. `workingBandFraction` prop name identical across `applyNormalization` (existing), `MemberMarksProps`, `MultiTracePlotProps`, page.
- **Risk:** the offset→band-fraction mapping is a UX approximation (the render core has no literal inter-trace gap model); chosen because it's the faithful, low-risk lever the existing normalization already exposes. Export pipeline is intentionally NOT changed (out of scope; export already honors band layout, not workingBandFraction).
