# R3-E — Series Finish Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Close issue #258 — apply Print tokens + finish/polish across the series-builder UI: rewrite `GroupingModeToggle` to Print vocabulary, recede the rail "COMPOSE" header, terracotta-ify the "Track reflections" checkbox, add a heatmap outer-hair keyline + rotated ordering-variable axis label, pin the cross-trace-tracking phase-color choice in a docstring, scope a leaky descendant selector, and apply the Plate Lift inset highlight + anchor the scoping Discard link.

**Architecture:** All work is frontend-only (React/TSX + Tailwind utility classes mapped to Print `@theme` tokens in `styles.css`). Visual/token findings are verified by `data-*`-attribute assertions (never Tailwind class-string tests — see `frontend/test/AGENTS.md`); the heatmap keyline is verified by inspecting the mark-factory's captured `Plot.rect` args (the established `MemberHeatmapLayer.test.tsx` pattern); the rotated axis label is a DOM overlay carrying a `data-testid` + the ordering-variable text, verified via RTL. The terracotta checkbox check-state is additionally verified visually by the orchestrator in the live harness.

**Tech Stack:** React 18, TypeScript, Tailwind (Print `@theme` tokens), Observable Plot (`@observablehq/plot`), Vitest + React Testing Library, JSDOM.

---

## Token / reference cheat-sheet (read before starting)

- Print active-segment vocabulary (canonical sibling `ScaleToggle.tsx:27/:36`): `bg-ink text-paper` active, `text-ink-faint` rest. Hover for a ghost toggle: `hover:text-ink hover:bg-paper-sunk` (DESIGN.md §5 Ghost).
- Terracotta accent utilities exist: `text-print-accent`, `border-print-accent`, `accent-print-accent` (native form-control accent; see `OffsetSlider.tsx:31`). The CSS var is `--color-print-accent` (= `--color-accent`, hue 38). Source: `styles.css:32,60`.
- Heatmap hairline keyline colour: `var(--color-hair)` (`styles.css:58`). This is what the mockup `series-builder.html:791-792` strokes the row frame in.
- Plate Lift inset highlight (DESIGN.md §4): `0 1px 0 rgba(255,255,255,.6) inset` — prepend to the existing `shadow-[...]` arbitrary value on each plate (the `.card` class at `styles.css:181` already carries it; the two series plates use inline `shadow-[...]` instead and currently omit it).
- Kicker letter-spacing for the receding rail header (R3-Y03): `tracking-[0.14em]` (matches the terracotta plate kicker at `SeriesBuilderPage.tsx:293`).
- `frontend/test/AGENTS.md`: assert `data-*` attributes, never Tailwind class strings; use RTL queries (`getByTestId`, `getByText().closest(...)`) over `document.querySelector`.

## File Structure

| File | Responsibility | Findings |
|---|---|---|
| `src/components/GroupingModeToggle.tsx` (modify) | Rewrite classNames to Print vocab; delete stale dark-era docstring | R3-Y02 |
| `test/GroupingModeToggle.test.tsx` (modify) | Add `data-active` reflection coverage (already exists; assert it survives the rewrite) | R3-Y02 |
| `src/components/SeriesBuilderRail.tsx` (modify) | Recede COMPOSE header; terracotta checkbox; scope `[&_input]` selector | R3-Y03, R3-Y04, R3-Y09 |
| `test/SeriesBuilderRail.test.tsx` (modify) | `data-*` assertions for recede + accent checkbox + scoped editor class | R3-Y03, R3-Y04, R3-Y09 |
| `src/components/MemberHeatmapLayer.tsx` (modify) | Emit per-row outer hair keyline `Plot.rect` | R3-Y06 |
| `test/MemberHeatmapLayer.test.tsx` (modify) | Assert the keyline rect's `stroke` + `fill:none` in captured args | R3-Y06 |
| `src/pages/SeriesBuilderPage.tsx` (modify) | Render rotated ordering-variable axis label in heatmap mode | R3-Y07, R3-Y10 |
| `test/SeriesBuilderPage.test.tsx` (modify) | Axis label present in heatmap mode / absent in waterfall; carries ordering var | R3-Y07 |
| `src/pages/SeriesScopingPage.tsx` (modify) | Plate Lift inset highlight; anchor Discard link padding | R3-Y10, R3-Y11 |
| `src/components/CrossTraceTrackingLayer.tsx` (modify) | Docstring note pinning phase-color (not terracotta) choice | R3-Y08 |

---

## Task 1: R3-Y02 — GroupingModeToggle to Print tokens

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/GroupingModeToggle.tsx:18-21, 59-65`
- Test: `packages/HimalayaUI/frontend/test/GroupingModeToggle.test.tsx`

The toggle ships five dark-era tokens (`bg-bg-subtle`, `text-fg-dim`, `text-fg`, `bg-bg-hover`, `border-border`) DESIGN.md §6 retires. Works today by R0a alias inheritance; rewrite to the canonical `ScaleToggle`/`RepresentationToggle` vocabulary. The component already exposes `data-active` per option — the TDD artifact is a `data-active` assertion (a `data-*` attribute, NOT a class string), confirming the active option is reflected after the rewrite.

- [ ] **Step 1: Add/confirm the failing test** — assert `data-active` reflects the value prop and toggles on change. (If equivalent coverage already exists, tighten it to also assert the inactive options report `data-active="false"`.)

```tsx
// test/GroupingModeToggle.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { GroupingModeToggle } from "../src/components/GroupingModeToggle";

describe("GroupingModeToggle — active reflection (Print tokens, R3-Y02)", () => {
  it("reflects the active option via data-active and routes onChange", () => {
    const onChange = vi.fn();
    const { rerender } = render(
      <GroupingModeToggle value="bySample" onChange={onChange} />,
    );
    const bySample = screen.getByRole("radio", { name: "By sample" });
    const byPhase = screen.getByRole("radio", { name: "By phase" });
    expect(bySample).toHaveAttribute("data-active", "true");
    expect(byPhase).toHaveAttribute("data-active", "false");

    fireEvent.click(byPhase);
    expect(onChange).toHaveBeenCalledWith("byPhase");

    rerender(<GroupingModeToggle value="byPhase" onChange={onChange} />);
    expect(screen.getByRole("radio", { name: "By phase" })).toHaveAttribute(
      "data-active",
      "true",
    );
  });
});
```

- [ ] **Step 2: Run to verify it passes today (guards the contract), then proceed to the token rewrite**

Run: `node_modules/.bin/vitest run test/GroupingModeToggle.test.tsx`
Expected: PASS (this test pins the behaviour we must preserve; the token rewrite is a refactor under it). Note: this is a refactor-under-test, not red→green — the visual token swap has no behavioural delta, so the data-* test guards against regression rather than driving new behaviour. The visual change is verified by the orchestrator's live harness.

- [ ] **Step 3: Rewrite the classNames + delete the stale docstring**

Replace the block-comment at lines 17-21 (the `**Styling — text-link parity...**` paragraph describing `text-fg-dim`/`bg-bg-hover`/`bg-bg-subtle`) — collapse to a one-line note:

```tsx
 * Spec selectors: container `data-testid="grouping-mode"` with
 * `data-mode={"bySample" | "byPhase" | "distinct"}` reflecting the
 * active option (used by E2E tests). Each option carries `data-value` and
 * `data-active` for assertions.
 *
 * Styling mirrors the canonical sibling toggles (`ScaleToggle`,
 * `RepresentationToggle`): `bg-ink text-paper` active, `text-ink-faint` at
 * rest, ghost hover (`hover:text-ink hover:bg-paper-sunk`). Print vocabulary
 * only — no dark-era `bg-bg-*`/`text-fg-*`/`border-border` tokens (DESIGN.md §6).
 */
```

Replace the `className` array at lines 59-65:

```tsx
            className={[
              "px-1.5 py-0.5 rounded text-xs transition-colors",
              "border border-transparent",
              active
                ? "bg-ink text-paper"
                : "text-ink-faint hover:text-ink hover:bg-paper-sunk",
            ].join(" ")}
```

- [ ] **Step 4: Run the test + grep for dark-era survivors in this file**

Run: `node_modules/.bin/vitest run test/GroupingModeToggle.test.tsx`
Expected: PASS.
Run: `grep -nE "bg-bg|text-fg|border-border|bg-bg-subtle|bg-bg-hover" src/components/GroupingModeToggle.tsx`
Expected: no matches.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/GroupingModeToggle.tsx \
        packages/HimalayaUI/frontend/test/GroupingModeToggle.test.tsx
git commit -m "fix(series): GroupingModeToggle uses Print tokens (R3-Y02, #258)"
```

---

## Task 2: R3-Y03 + R3-Y04 + R3-Y09 — Rail header recede, terracotta checkbox, scoped selector

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/SeriesBuilderRail.tsx:78-89, 92, 127-139`
- Test: `packages/HimalayaUI/frontend/test/SeriesBuilderRail.test.tsx`

Three scoped changes in one file:
- **R3-Y03** — the "Compose" header at line 79 uses `text-xs font-semibold uppercase tracking-wide text-ink`; it competes with the figure-plate kicker. Drop `text-ink`→`text-ink-faint`, `tracking-wide`→`tracking-[0.14em]`, remove `font-semibold`.
- **R3-Y04** — the native checkbox at lines 131-137 renders OS-default blue. Add `accent-print-accent` (mirrors `OffsetSlider.tsx:31`) so the checked state is terracotta.
- **R3-Y09** — the editor `<section>` at line 92 uses the wildcard descendant selector `[&_input]:bg-plate`, which also catches `<input type="range">` slider thumbs nested in injected `editControls`. Replace with a scoped class on a wrapper.

The "Compose" header and the section/checkbox have no stable test hook today. The TDD artifacts are `data-*` attributes: add `data-testid="rail-compose-header"` to the header span and assert the checkbox input carries `data-accent="print"` (a data-* marker the component sets alongside the `accent-print-accent` class so the test never asserts the class string). For R3-Y09, assert the editor wrapper exposes `data-rail-edit-inputs` and that the wildcard selector is gone (grep).

- [ ] **Step 1: Write the failing tests**

```tsx
// test/SeriesBuilderRail.test.tsx — append to the existing describe block
import { render, screen } from "@testing-library/react";
// (existing imports / makeProps helper already in the file)

it("recedes the COMPOSE header behind the plate (R3-Y03)", () => {
  render(<SeriesBuilderRail {...makeProps()} />);
  // The header is a recessed label, not a competing title: it carries the
  // receded data marker the component sets when it drops to ink-faint.
  expect(screen.getByTestId("rail-compose-header")).toHaveAttribute(
    "data-recede",
    "true",
  );
});

it("renders the Track-reflections checkbox in the terracotta accent (R3-Y04)", () => {
  render(<SeriesBuilderRail {...makeProps()} />);
  expect(screen.getByTestId("track-toggle-input")).toHaveAttribute(
    "data-accent",
    "print",
  );
});

it("scopes edit-input styling so it does not target slider thumbs (R3-Y09)", () => {
  render(<SeriesBuilderRail {...makeProps({ editControls: <input data-testid="ec" /> })} />);
  expect(screen.getByTestId("rail-edit")).toHaveAttribute("data-rail-edit-inputs", "");
});
```

> If the file's existing `makeProps` helper has a different name/signature, adapt these three tests to it (it is the helper that builds the `SeriesBuilderRailProps`; the file already constructs full props for the 10 passing tests). Do NOT invent a new helper.

- [ ] **Step 2: Run to verify they fail**

Run: `node_modules/.bin/vitest run test/SeriesBuilderRail.test.tsx`
Expected: 3 FAIL — missing `data-recede`, `data-accent`, `data-rail-edit-inputs`.

- [ ] **Step 3: Implement the three changes**

Header span (lines 78-79):

```tsx
      <div className="flex items-center justify-between">
        <span
          data-testid="rail-compose-header"
          data-recede="true"
          className="text-xs uppercase tracking-[0.14em] text-ink-faint"
        >Compose</span>
```

Editor section (line 92) — replace the wildcard with a scoped wrapper class + data marker:

```tsx
        <section
          className="flex flex-col gap-1.5 rail-edit-inputs"
          data-testid="rail-edit"
          data-rail-edit-inputs=""
        >
          {editControls}
        </section>
```

Add the scoped rule to `styles.css` (so only non-range inputs in the editor pick up the plate fill):

```css
  /* Series-builder rail recipe-editor: plate-fill text inputs only — never
     the range-slider thumbs that injected editControls may also nest
     (R3-Y09, #258). */
  .rail-edit-inputs input:not([type="range"]):not([type="checkbox"]) {
    background: var(--color-plate);
  }
```

Checkbox input (lines 131-137):

```tsx
          <input
            type="checkbox"
            data-testid="track-toggle-input"
            data-accent="print"
            checked={trackOn}
            onChange={(e) => onTrackOnChange(e.target.checked)}
            className="rounded border-hair accent-print-accent"
          />
```

- [ ] **Step 4: Run the tests + grep the wildcard is gone**

Run: `node_modules/.bin/vitest run test/SeriesBuilderRail.test.tsx`
Expected: PASS (13 tests).
Run: `grep -n "\[&_input\]" src/components/SeriesBuilderRail.tsx`
Expected: no matches.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/SeriesBuilderRail.tsx \
        packages/HimalayaUI/frontend/src/styles.css \
        packages/HimalayaUI/frontend/test/SeriesBuilderRail.test.tsx
git commit -m "fix(series): recede rail header, terracotta track-checkbox, scope edit-input selector (R3-Y03/Y04/Y09, #258)"
```

---

## Task 3: R3-Y06 — Heatmap outer-hair keyline

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MemberHeatmapLayer.tsx:113-115, 155-165`
- Test: `packages/HimalayaUI/frontend/test/MemberHeatmapLayer.test.tsx`

Heatmap rows currently push one `Plot.rect` (the cells) with `stroke: null`. Add a second `Plot.rect` per row: a single-cell outer frame at the row's inset bounds, `fill: none`, `stroke: var(--color-hair)`, `strokeWidth: 1` — mirroring the mockup's row keyline (`series-builder.html:791-792`). This makes the row read as "intensity inside a frame" instead of "two coloured boxes". The test asserts on the captured `Plot.rect` args (the established pattern — `Plot.rect` is mocked to echo its args).

- [ ] **Step 1: Write the failing test**

```tsx
// test/MemberHeatmapLayer.test.tsx — add inside the buildMemberHeatmapMarks describe
it("emits an outer hair keyline rect framing each row (R3-Y06)", () => {
  const marks = buildMemberHeatmapMarks({
    member: makeMember(),
    trace,
    yBand: [10, 30],
    qDomain: [0.05, 0.9],
    groupingMode: "byPhase",
    allMembers: [makeMember()],
    sampleIdFor: () => 1,
  });
  // Two rect marks: the binned cells + the framing keyline.
  expect(marks.length).toBe(2);
  const calls = (Plot.rect as unknown as { mock: { calls: unknown[][] } }).mock.calls;
  // The keyline is the call with fill:none + a hair stroke.
  const keyline = calls.find(
    (c) => (c[1] as { fill?: unknown; stroke?: unknown }).fill === "none",
  );
  expect(keyline).toBeDefined();
  const opts = keyline![1] as { stroke: string; strokeWidth: number };
  expect(opts.stroke).toBe("var(--color-hair)");
  expect(opts.strokeWidth).toBe(1);
});
```

- [ ] **Step 2: Run to verify it fails**

Run: `node_modules/.bin/vitest run test/MemberHeatmapLayer.test.tsx`
Expected: FAIL — only one rect today (`marks.length` is 1), no `fill:"none"` call.

- [ ] **Step 3: Append the keyline mark after the cells push (after line 165)**

```tsx
  // Outer hair keyline framing the row (R3-Y06, #258). Without it adjacent
  // rows read as "coloured boxes"; the frame makes them "intensity rows in a
  // frame", matching the mockup's row keyline (series-builder.html:791-792).
  // A single full-span cell at the row's inset bounds, no fill, hair stroke.
  marks.push(
    Plot.rect([{ x1: qDomain[0], x2: qDomain[1] }], {
      x1: "x1",
      x2: "x2",
      y1: y1,
      y2: y2,
      fill: "none",
      stroke: "var(--color-hair)",
      strokeWidth: 1,
    }),
  );

  return marks;
```

(Remove the existing bare `return marks;` at line 167 so the keyline pushes before the single return.)

- [ ] **Step 4: Run the test**

Run: `node_modules/.bin/vitest run test/MemberHeatmapLayer.test.tsx`
Expected: PASS (existing tests + the new keyline test). Confirm the existing "emits no marks when trace missing/empty" test still returns `marks.length === 0` (the early-return at line 109 is before the keyline push — verify the keyline lives AFTER the cells, inside the populated path).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/MemberHeatmapLayer.tsx \
        packages/HimalayaUI/frontend/test/MemberHeatmapLayer.test.tsx
git commit -m "fix(series): heatmap rows get an outer hair keyline frame (R3-Y06, #258)"
```

---

## Task 4: R3-Y07 — Rotated ordering-variable axis label (heatmap mode)

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx:356-381`
- Test: `packages/HimalayaUI/frontend/test/SeriesBuilderPage.test.tsx`

The heatmap lacks the mockup's rotated ordering-variable axis label in the left margin (`series-builder.html:817-822` `.axis-title`), which is what makes the representation read as "a migration map ordered by X" instead of "stacked rows". Render a DOM overlay (a rotated `<div>` positioned in the plot column's left margin) ONLY in heatmap mode, carrying the active series's `ordering_variable`. A DOM overlay (not an injected Plot mark) keeps it RTL-testable and avoids threading label state into `MultiTracePlot`'s Plot config.

The plot column is `<div ref={plotColRef} className="min-w-0 flex-1">` at lines 357-372. Wrap it so the overlay is positioned relative to the plot column.

- [ ] **Step 1: Write the failing test**

```tsx
// test/SeriesBuilderPage.test.tsx — add to the read+states describe
// (reuse the file's existing render helper that mounts SeriesBuilderPage with
//  a populated series; switch representation to heatmap via the rail's
//  RepresentationToggle, or seed representation if the helper supports it.)

it("shows a rotated ordering-variable axis label in heatmap mode (R3-Y07)", async () => {
  renderBuilder(/* populated series whose ordering_variable === "dose" */);
  // switch to heatmap (RepresentationToggle exposes data-value="heatmap")
  fireEvent.click(screen.getByRole("button", { name: /heatmap/i }));
  const axis = await screen.findByTestId("heatmap-axis-title");
  expect(axis).toHaveTextContent("dose");
});

it("omits the heatmap axis label in waterfall mode (R3-Y07)", () => {
  renderBuilder(/* same populated series, default waterfall */);
  expect(screen.queryByTestId("heatmap-axis-title")).toBeNull();
});
```

> Use the file's existing render helper + fixture (the 18 passing tests already mount a populated series and toggle representation). If `ordering_variable` is null in the default fixture, set it to `"dose"` in the test's series override. Match the existing helper's name/signature — do not invent `renderBuilder` if the file uses a different one.

- [ ] **Step 2: Run to verify it fails**

Run: `node_modules/.bin/vitest run test/SeriesBuilderPage.test.tsx`
Expected: FAIL — `heatmap-axis-title` not found.

- [ ] **Step 3: Wrap the plot column with a relative container + conditional overlay**

Replace the plot-column `<div>` (lines 357-372) with:

```tsx
                    <div ref={plotColRef} className="relative min-w-0 flex-1">
                      {representation === "heatmap" && s.ordering_variable && (
                        // Rotated ordering-variable axis title in the left
                        // margin (R3-Y07, #258). Mockup series-builder.html
                        // :817-822 `.axis-title` — what makes the heatmap read
                        // as a migration map "ordered by X", not stacked rows.
                        <div
                          data-testid="heatmap-axis-title"
                          aria-hidden="true"
                          className="pointer-events-none absolute left-0 top-1/2 z-10 origin-left -translate-y-1/2 -rotate-90 whitespace-nowrap text-[10.5px] uppercase tracking-[0.14em] text-ink-faint"
                        >
                          {s.ordering_variable} &rarr;
                        </div>
                      )}
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
                        representation={representation}
                        showCrossTraceTracking={trackOn}
                      />
                    </div>
```

> Note: `s` is the series object already in scope (used at `:332` `s.title`, `:392` `s.ordering_variable`). `representation` is the local state at `:187`. No new state/props.

- [ ] **Step 4: Run the test**

Run: `node_modules/.bin/vitest run test/SeriesBuilderPage.test.tsx`
Expected: PASS (existing 18 + 2 new). Confirm no regression in the waterfall-mode tests.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx \
        packages/HimalayaUI/frontend/test/SeriesBuilderPage.test.tsx
git commit -m "feat(series): rotated ordering-variable axis label in heatmap mode (R3-Y07, #258)"
```

---

## Task 5: R3-Y10 + R3-Y11 — Plate Lift inset highlight + anchored Discard link

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesScopingPage.tsx:228, 236`
- Modify: `packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx:277`

Two pages' plates use inline `shadow-[...]` and omit the DESIGN.md §4 Plate Lift `0 1px 0 rgba(255,255,255,.6) inset` top-rule. Prepend it. Also anchor the scoping Discard link with the mockup's `padding: 7px 4px` (`py-[7px] px-[4px]`) so it doesn't float.

These are pure visual token applications with no behavioural delta. No new test drives them (a shadow value is not a `data-*` contract and class-string assertions are forbidden by `frontend/test/AGENTS.md`); they are verified by the orchestrator's live harness + the existing scoping/builder tests not regressing. Run the affected suites as the verification step.

- [ ] **Step 1: Run the affected suites green first (baseline guard)**

Run: `node_modules/.bin/vitest run test/scoping.test.tsx test/SeriesBuilderPage.test.tsx`
Expected: PASS (establishes the no-regression floor before the edit).

- [ ] **Step 2: Apply the scoping plate inset + Discard padding**

`SeriesScopingPage.tsx:236` — prepend the inset to the shadow arbitrary value:

```tsx
        className="w-full max-w-[760px] rounded-md border border-hair bg-plate px-8 py-7 shadow-[0_1px_0_rgba(255,255,255,.6)_inset,0_1px_1px_rgba(60,52,40,.04),0_18px_40px_-22px_rgba(60,52,40,.22)]"
```

`SeriesScopingPage.tsx:228` — anchor the Discard link (mockup padding `7px 4px`):

```tsx
          className="px-[4px] py-[7px] text-xs font-semibold text-ink-faint hover:text-ink"
```

- [ ] **Step 3: Apply the builder plate inset**

`SeriesBuilderPage.tsx:277` — prepend the inset to the shadow arbitrary value:

```tsx
            className={`w-full ${collapsed ? "max-w-[1336px]" : "max-w-[1180px]"} rounded border border-hair bg-plate p-8 shadow-[0_1px_0_rgba(255,255,255,.6)_inset,0_1px_1px_rgba(60,52,40,.04),0_18px_40px_-20px_rgba(60,52,40,.22)] transition-[max-width] duration-200`}
```

- [ ] **Step 4: Run the affected suites again (no regression)**

Run: `node_modules/.bin/vitest run test/scoping.test.tsx test/SeriesBuilderPage.test.tsx`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/pages/SeriesScopingPage.tsx \
        packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx
git commit -m "fix(series): Plate Lift inset highlight + anchored scoping Discard link (R3-Y10/Y11, #258)"
```

---

## Task 6: R3-Y08 — CrossTraceTrackingLayer phase-color docstring

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/CrossTraceTrackingLayer.tsx:119` (the `const stroke = phaseColor(phase);` line, inside the per-phase loop)

Pin the "phase color, not terracotta" choice in a docstring so a future reviewer doesn't "fix" the migration line to terracotta thinking it's the grease-pencil interaction mark. The migration line is a phase relation, not an interaction. Docstring-only; no test.

- [ ] **Step 1: Add the comment above `const stroke = phaseColor(phase);` (line 119)**

```tsx
    // Phase colour, NOT terracotta. The migration polyline is a *phase
    // relation* (the same reflection tracked across members of one phase),
    // not an interaction mark — so it carries the phase hue, like every other
    // phase-coloured mark on the surface. DESIGN.md §2 Phase-Carries-the-
    // Surface Rule reserves terracotta (the grease pencil) for interaction
    // (reject ✕, primary action, q-link cross-highlight). Do not "fix" this
    // to the accent (cf. R3-F01: terracotta over-use is the recurring miss).
    const stroke = phaseColor(phase);
```

- [ ] **Step 2: Verify the file still type-checks (no behavioural change)**

Run: `node_modules/.bin/vitest run test/CrossTraceTrackingLayer.test.tsx`
Expected: PASS (unchanged — comment-only).

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/CrossTraceTrackingLayer.tsx
git commit -m "docs(series): pin phase-color (not terracotta) choice in CrossTraceTrackingLayer (R3-Y08, #258)"
```

---

## Task 7: Verify — full frontend gate

**Files:** none (verification only)

- [ ] **Step 1: Full Vitest suite**

Run: `npm test > /tmp/r3e-vitest.out 2>&1; grep -E "Test Files|Tests |FAIL" /tmp/r3e-vitest.out`
Expected: all test files pass, 0 failures. Floor: no fewer passing tests than baseline (62 across the six target files) plus the new tests added in Tasks 1–4.

- [ ] **Step 2: Build gate (tsc --noEmit + vite build)**

Run: `npm run build 2>&1 | tail -20`
Expected: `tsc --noEmit` clean, vite build succeeds (no type errors from the new overlay / mark args).

- [ ] **Step 3: Project-wide dark-era grep on touched files (the R3-F05 pattern fix in miniature)**

Run: `grep -nE "bg-bg|text-fg|border-border|bg-bg-subtle|bg-bg-hover" src/components/GroupingModeToggle.tsx src/components/SeriesBuilderRail.tsx`
Expected: no matches (confirms the token sweep is complete in the files this issue scopes).

- [ ] **Step 4: Report to team-lead for PR stage** (do not open the PR until directed by the runbook's PR stage).

---

## Self-Review

**Spec coverage (issue #258 Done-when):**
- GroupingModeToggle uses Print tokens; stale docstring deleted → Task 1. ✓
- Rail "COMPOSE" header recedes → Task 2 (R3-Y03). ✓
- "Track reflections" checkbox check-state is terracotta → Task 2 (R3-Y04); visual verified by orchestrator. ✓
- Heatmap rows render an outer hair keyline → Task 3 (R3-Y06). ✓
- Heatmap mode shows a rotated ordering-variable axis label → Task 4 (R3-Y07). ✓
- CrossTraceTrackingLayer docstring cites §2 → Task 6 (R3-Y08). ✓
- Plate Lift inset highlight on scoping + builder plates → Task 5 (R3-Y10). ✓
- Scoping Discard link anchored per mockup padding → Task 5 (R3-Y11). ✓
- R3-Y09 scoped descendant selector (in scope per issue body) → Task 2. ✓
- Out of scope: R3-Y05 folio chips — not touched. ✓

**Placeholder scan:** every code step shows the actual class strings / mark args / comment text. Render-helper names in Tasks 2/4 are flagged as "adapt to the existing helper" because the test files own those helpers; the plan does not invent new ones.

**Type consistency:** `data-recede`, `data-accent`, `data-rail-edit-inputs`, `heatmap-axis-title`, `rail-compose-header` are the new `data-*` hooks; each is introduced in the component and asserted in the matching test with the same name. The keyline mark uses the same `Plot.rect(data, opts)` shape the existing cells mark and the test mock use. `s.ordering_variable` + `representation` reuse identifiers already in `SeriesBuilderPage` scope.

**Visual "terracotta checkbox check-state" verification (Done-when):** the TDD artifact is the `data-accent="print"` attribute assertion (Task 2) — a `data-*` test per `frontend/test/AGENTS.md`, never a class string. The actual rendered terracotta check-state is verified separately by the orchestrator in the live Playwright harness (the checkbox uses `accent-print-accent`, the native form-control accent token confirmed present in `OffsetSlider.tsx:31` + `styles.css`).
