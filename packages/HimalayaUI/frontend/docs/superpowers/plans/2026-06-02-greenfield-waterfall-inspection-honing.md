# WaterfallChart Inspection Honing — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Fresh implementer per task + two-stage review. Steps use `- [ ]`.
> Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Spec: `docs/superpowers/specs/2026-06-02-greenfield-waterfall-inspection-honing.md`

**Goal:** Add a peak-q cursor (full-height guide + q readout, snapped to detected peaks) and a live log/linear q-axis toggle to `WaterfallChart`.

**Architecture:** Centralized pointer handling at the chart level (TracePlot untouched — per-row interaction would trap scroll via PlotFrame's unconditional wheel `preventDefault`). A pure snap helper drives a controllable `hoveredQ`; the guide reuses the shared `makeAxis(qDomain, [4, plotW-4], scale)` so it aligns with beads and axis ticks. Scale and hoveredQ both follow the file's existing controlled-or-internal idiom (`hoveredKey`/`internalHot`).

**Tech stack:** React + TS strict (`exactOptionalPropertyTypes`), Vitest/RTL, `SegmentedControl` primitive, Tailwind token classes (`bg-accent`, `bg-plate`, `border-hair`, `text-meta`, `font-mono`). All changes are design-guard non-exempt → placement-only.

**Working dir:** `packages/HimalayaUI/frontend`. Run `npx tsc --noEmit -p tsconfig.build.json` (not the default root tsc). Commit only the named files (never `git add -A`).

---

### Task 1: Pure snap helper (`cursor.ts`)

**Files:**
- Create: `src/print/waterfall/cursor.ts`
- Test: `test/print-waterfall/cursor.test.ts`

- [ ] **Step 1: Write the failing test** — `test/print-waterfall/cursor.test.ts`

```ts
import { describe, it, expect } from "vitest";
import { cursorCandidates, snapToPeakQ } from "../../src/print/waterfall/cursor";
import type { WaterfallRow } from "../../src/print/waterfall/waterfallModel";

const EMPTY = { q: [], I: [], sigma: [] };
function row(key: string, anchors: { id: number; q: number }[]): WaterfallRow {
  return {
    key, label: key, trace: EMPTY, phase: "Ia3d", state: "indexed",
    anchors: anchors.map((a) => ({ id: a.id, q: a.q, intensity: 100, phase: "Ia3d" })),
    bandHeight: 1, yOffset: 0,
  };
}
// Identity projection: q→px is q*1000 (so q=0.06 → 60px), tol in px.
const toPx = (q: number) => q * 1000;

describe("cursor", () => {
  it("flattens anchors across all rows into {id,q} candidates", () => {
    const rows = [row("a", [{ id: 1, q: 0.06 }]), row("b", [{ id: 2, q: 0.08 }, { id: 3, q: 0.10 }])];
    expect(cursorCandidates(rows)).toEqual([
      { id: 1, q: 0.06 }, { id: 2, q: 0.08 }, { id: 3, q: 0.10 },
    ]);
  });

  it("snaps to the nearest peak q within tolerance", () => {
    const rows = [row("a", [{ id: 1, q: 0.06 }, { id: 2, q: 0.09 }])];
    // pointer at px 62 (q≈0.062); nearest is q=0.06 at px 60, dist 2 ≤ 10 → snaps to 0.06
    expect(snapToPeakQ(62, rows, toPx, 10)).toBe(0.06);
  });

  it("returns null when no peak is within tolerance", () => {
    const rows = [row("a", [{ id: 1, q: 0.06 }])];
    // pointer at px 100 (q≈0.1); nearest peak at px 60, dist 40 > 10 → null
    expect(snapToPeakQ(100, rows, toPx, 10)).toBeNull();
  });

  it("returns null for empty rows", () => {
    expect(snapToPeakQ(60, [], toPx, 10)).toBeNull();
  });

  it("ignores optimistic peaks (id < 0)", () => {
    const rows = [row("a", [{ id: -1, q: 0.06 }])];
    expect(snapToPeakQ(60, rows, toPx, 10)).toBeNull();
  });
});
```

- [ ] **Step 2: Run → fail** — `npm test -- print-waterfall/cursor` → FAIL (module not found).

- [ ] **Step 3: Implement** — `src/print/waterfall/cursor.ts`

```ts
/**
 * cursor — pure snap math for the WaterfallChart peak-q cursor. A pointer x
 * (container-relative px) snaps to the nearest detected peak's q across all
 * rows, within a px tolerance, reusing the plot engine's hitTestPeaks (which
 * ignores optimistic id<0 peaks and resolves ties to the later peak).
 */
import { hitTestPeaks } from "../plot/interaction";
import type { WaterfallRow } from "./waterfallModel";

/** All anchors across all rows as hit-test candidates. id+q suffice; equal q → equal px. */
export function cursorCandidates(rows: WaterfallRow[]): { id: number; q: number }[] {
  return rows.flatMap((r) => r.anchors.map((a) => ({ id: a.id, q: a.q })));
}

/** Snap a container-relative px to the nearest peak q within tolPx; null if none. */
export function snapToPeakQ(
  px: number,
  rows: WaterfallRow[],
  toPx: (q: number) => number,
  tolPx: number,
): number | null {
  const hit = hitTestPeaks(cursorCandidates(rows), px, toPx, tolPx);
  return hit ? hit.q : null;
}
```

- [ ] **Step 4: Run → pass** — `npm test -- print-waterfall/cursor` → PASS (5).
- [ ] **Step 5: tsc** — `npx tsc --noEmit -p tsconfig.build.json` → 0.
- [ ] **Step 6: Commit** — `git add src/print/waterfall/cursor.ts test/print-waterfall/cursor.test.ts && git commit`.

---

### Task 2: Log/linear q-axis toggle

**Files:**
- Modify: `src/print/waterfall/WaterfallChart.tsx`
- Test: `test/print-waterfall/WaterfallChart.test.tsx` (append)

Mirror the existing `hoveredKey`/`internalHot` controlled-or-internal idiom for scale.

- [ ] **Step 1: Write failing tests** — append to `test/print-waterfall/WaterfallChart.test.tsx` (these reference `TRANSITION`/`FULL`-style fixtures already imported at the top of that file; reuse whatever fixture the existing tests use — check the imports and use the same `ROWS`/fixture variable).

```tsx
import { fireEvent } from "@testing-library/react"; // add if not already imported

it("defaults the q-axis to log and exposes it via data-xtype", () => {
  const { getByTestId } = render(<WaterfallChart rows={ROWS} />);
  expect(getByTestId("waterfall").getAttribute("data-xtype")).toBe("log");
});

it("renders a log/linear scale toggle and flips internal scale on click", () => {
  const { getByTestId } = render(<WaterfallChart rows={ROWS} />);
  expect(getByTestId("wf-scale")).toBeTruthy();
  // SegmentedControl segments carry data-value; click the linear one.
  const linear = getByTestId("wf-scale").querySelector('[data-value="linear"]') as HTMLElement;
  fireEvent.click(linear);
  expect(getByTestId("waterfall").getAttribute("data-xtype")).toBe("linear");
});

it("fires onXTypeChange and respects a controlled xType", () => {
  const onXTypeChange = vi.fn();
  const { getByTestId, rerender } = render(
    <WaterfallChart rows={ROWS} xType="log" onXTypeChange={onXTypeChange} />,
  );
  const linear = getByTestId("wf-scale").querySelector('[data-value="linear"]') as HTMLElement;
  fireEvent.click(linear);
  expect(onXTypeChange).toHaveBeenCalledWith("linear");
  // controlled: stays log until the parent updates the prop
  expect(getByTestId("waterfall").getAttribute("data-xtype")).toBe("log");
  rerender(<WaterfallChart rows={ROWS} xType="linear" onXTypeChange={onXTypeChange} />);
  expect(getByTestId("waterfall").getAttribute("data-xtype")).toBe("linear");
});
```

(If `ROWS` isn't the existing fixture name, use the one the file already imports — do not invent a new fixture. Verify `vi`/`fireEvent` imports exist; add only what's missing.)

- [ ] **Step 2: Run → fail** — `npm test -- print-waterfall/WaterfallChart` → new tests FAIL.

- [ ] **Step 3: Implement** in `WaterfallChart.tsx`:
  1. Import: `import { SegmentedControl } from "../ui";`
  2. Add to `WaterfallChartProps`: `onXTypeChange?: (next: "log" | "linear") => void;` (keep `xType?` as the controlled override).
  3. In the component: `const [internalXType, setInternalXType] = useState<"log" | "linear">("log");` then `const scale = xType ?? internalXType;`
  4. Replace every later use of `xType` with `scale` (per-row `TracePlot xType={scale}`, bottom `makeAxis(qDomain, [4, plotW - 4], scale)`).
  5. Add a change handler:
     ```tsx
     const setScale = (next: "log" | "linear"): void => {
       if (xType === undefined) setInternalXType(next);
       onXTypeChange?.(next);
     };
     ```
  6. Add `data-xtype={scale}` to the root `<div ... data-testid="waterfall">`.
  7. Add a thin toolbar row ABOVE the stack `.relative` div (placement-only):
     ```tsx
     <div className="flex justify-end mb-1">
       <SegmentedControl<"log" | "linear">
         options={[{ value: "log", label: "log" }, { value: "linear", label: "linear" }]}
         value={scale}
         onChange={setScale}
         size="sm"
         aria-label="q axis scale"
         testId="wf-scale"
       />
     </div>
     ```

- [ ] **Step 4: Run → pass** — `npm test -- print-waterfall/WaterfallChart` → all green (existing 8 + new 3).
- [ ] **Step 5: Guards** — `npm run lint:design` clean; `npx tsc --noEmit -p tsconfig.build.json` → 0.
- [ ] **Step 6: Commit** — `git add src/print/waterfall/WaterfallChart.tsx test/print-waterfall/WaterfallChart.test.tsx && git commit`.

---

### Task 3: Peak-q cursor (guide + readout)

**Files:**
- Modify: `src/print/waterfall/WaterfallChart.tsx`
- Test: `test/print-waterfall/WaterfallChart.test.tsx` (append)

Controllable `hoveredQ` (so rendering is testable without JSDOM rect math) + a centralized pointer handler that snaps and broadcasts.

- [ ] **Step 1: Write failing tests** — append:

```tsx
it("renders no q-guide when hoveredQ is unset", () => {
  const { queryByTestId } = render(<WaterfallChart rows={ROWS} />);
  expect(queryByTestId("wf-qguide")).toBeNull();
  expect(queryByTestId("wf-qreadout")).toBeNull();
});

it("renders the q-guide + readout at a controlled hoveredQ", () => {
  // pick a q that exists among ROWS anchors:
  const q = ROWS[0]!.anchors[0]!.q;
  const { getByTestId } = render(<WaterfallChart rows={ROWS} hoveredQ={q} />);
  const guide = getByTestId("wf-qguide");
  expect(guide.getAttribute("data-q")).toBe(String(q));
  expect(getByTestId("wf-qreadout").textContent).toBe(q.toFixed(3));
});

it("fires onHoverQ on pointer move near a peak (mocked rect)", () => {
  const onHoverQ = vi.fn();
  const { getByTestId } = render(<WaterfallChart rows={ROWS} onHoverQ={onHoverQ} />);
  const stack = getByTestId("wf-stack");
  // plotW falls back to maxWidth-LABEL_W in JSDOM; mock the rect so px math resolves.
  stack.getBoundingClientRect = () => ({ left: 0, top: 0, width: 1080, height: 420, right: 1080, bottom: 420, x: 0, y: 0, toJSON: () => ({}) }) as DOMRect;
  // Move to the px of the first anchor under the shared projection. We don't assert the
  // exact value (projection detail); we assert SOME q was broadcast for an on-peak move.
  const q = ROWS[0]!.anchors[0]!.q;
  // Approx px for a log/linear axis over qDomain → use a move far from any peak first:
  fireEvent.pointerMove(stack, { clientX: 5 });   // left edge, unlikely on a peak
  fireEvent.pointerLeave(stack);
  expect(onHoverQ).toHaveBeenLastCalledWith(undefined);
});
```

(The third test asserts the leave-path broadcast deterministically; the on-peak px is projection-dependent, so it only checks that pointer leave clears. Keep it — it pins the wiring. If `ROWS[0].anchors` is empty for the fixture used, switch to a fixture/row that has anchors, e.g. the indexed member.)

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement** in `WaterfallChart.tsx`:
  1. Imports: `import { snapToPeakQ } from "./cursor";` and `import { PEAK_HIT_PX } from "../plot/interaction";`
  2. Add to props: `hoveredQ?: number;` (controlled override). `onHoverQ` is ALREADY in the props interface — destructure it now (it is currently dead).
  3. State + effective value: `const [internalHoveredQ, setInternalHoveredQ] = useState<number | undefined>(undefined);` then `const cursorQ = hoveredQ ?? internalHoveredQ;`
  4. Shared projection (build once from the same args as the bottom axis):
     `const sharedX = makeAxis(qDomain, [4, plotW - 4], scale);`
     Use `sharedX` for BOTH the bottom `<Axis axis={sharedX} ...>` (replace the inline `makeAxis(...)` there) and the cursor.
  5. Setter:
     ```tsx
     const setCursorQ = (q: number | undefined): void => {
       if (hoveredQ === undefined) setInternalHoveredQ(q);
       onHoverQ?.(q);
     };
     ```
  6. Pointer handlers on the stack `.relative` div (give it `data-testid="wf-stack"`):
     ```tsx
     onPointerMove={(e) => {
       const rect = e.currentTarget.getBoundingClientRect();
       const px = e.clientX - rect.left;
       const snapped = snapToPeakQ(px, rows, (q) => sharedX.to(q), PEAK_HIT_PX);
       setCursorQ(snapped ?? undefined);
     }}
     onPointerLeave={() => setCursorQ(undefined)}
     ```
  7. Cursor overlay — render INSIDE the `.relative` stack div, AFTER the rows map, placement-only token classes, only when `cursorQ != null`:
     ```tsx
     {cursorQ != null ? (
       <>
         <div
           data-role="wf-qguide"
           data-testid="wf-qguide"
           data-q={cursorQ}
           className="absolute top-0 w-px bg-accent pointer-events-none"
           style={{ left: sharedX.to(cursorQ), height: TOTAL_H }}
         />
         <div
           data-role="wf-qreadout"
           data-testid="wf-qreadout"
           className="absolute top-0 -translate-x-1/2 bg-plate border border-hair rounded-sm px-1 text-meta font-mono text-ink pointer-events-none"
           style={{ left: sharedX.to(cursorQ) }}
         >
           {cursorQ.toFixed(3)}
         </div>
       </>
     ) : null}
     ```

- [ ] **Step 4: Run → pass** — `npm test -- print-waterfall` → all green.
- [ ] **Step 5: Guards** — `npm run lint:design` clean; `npx tsc --noEmit -p tsconfig.build.json` → 0.
- [ ] **Step 6: Commit** — `git add src/print/waterfall/WaterfallChart.tsx test/print-waterfall/WaterfallChart.test.tsx && git commit`.

---

### Task 4: Stories — exercise the toggle + show the cursor

**Files:**
- Modify: `src/print/waterfall/WaterfallChart.stories.tsx`

- [ ] **Step 1: Update stories** (no failing test — verified by build-storybook + manual review):
  - Ensure `Wide` and `Narrow` do NOT pass `xType` (so the live toggle defaults to log and the user can flip it). `LinearQ` keeps `xType="linear"` to demonstrate the controlled path.
  - Add a `QCursor` story that passes a controlled `hoveredQ` equal to a real anchor q from the fixture (e.g. the first indexed member's first anchor), so the guide + readout are visible statically in Storybook:
    ```tsx
    export const QCursor: Story = {
      render: () => {
        const q = TRANSITION[0]?.anchors[0]?.q ?? 0.06;
        return <div className={frame} style={{ width: 920 }}><WaterfallChart rows={TRANSITION} maxWidth={920} hoveredQ={q} /></div>;
      },
    };
    ```
    (Use whatever fixture the file already imports; if `TRANSITION[0]` has no anchors, pick an indexed member that does. Use conditional spread if `exactOptionalPropertyTypes` complains about an optional prop.)
  - Keep `MixedStates` and `ControlledHover` as-is.

- [ ] **Step 2: tsc** — `npx tsc --noEmit -p tsconfig.build.json` → 0.
- [ ] **Step 3: Commit** — `git add src/print/waterfall/WaterfallChart.stories.tsx && git commit`.

---

### Final verification (after all tasks)

- [ ] `npm test -- print-waterfall` → all green.
- [ ] `npm run lint:design` → clean (proves placement-only).
- [ ] `npx tsc --noEmit -p tsconfig.build.json` → 0.
- [ ] `npm run build` → exit 0.
- [ ] `npm run build-storybook` → exit 0.
- [ ] Manual (user): `npm run storybook` → Wide story, flip log/linear, hover a peak to see the guide + q readout; QCursor story shows the static guide.
