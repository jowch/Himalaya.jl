# WaterfallChart Inspection Honing — Design Spec

**Date:** 2026-06-02
**Status:** Approved scope, pending spec review
**Surface:** `src/print/waterfall/WaterfallChart.tsx` (greenfield "The Print" Series renderer)

## Why

A SAXS waterfall is a **peak-position comparison instrument**. A due-diligence research pass (domain literature + our own code + the mockups) confirmed the two highest-value reading tasks are:

1. **Read a peak's exact q** — "which peaks did we observe and where."
2. **Judge whether the same peak is present/aligned across the stack** — the visual precursor to "how does scattering change across the experimental treatment."

Today the chart serves *neither*: the rows pass no `interaction` to their `TracePlot`s, so the existing hover q-readout never fires, and `WaterfallChart`'s `onHoverQ` prop is declared but dead (`WaterfallChart.tsx:20`). This spec adds a **peak-q cursor** that delivers both readings in one interaction, plus a **user-facing log/linear q-axis toggle**.

The richer **migration track** (a connector following one reflection across the stack with ghost rings for predicted-but-absent orders) was researched and **deliberately deferred**: even though it is buildable client-side from `confirmed_index.phase + lattice_d` (the legacy `buildAnchorMap` proves this — it is *not* blocked on backend schema as the prior spec claimed), the q-proximity matching is fragile and the resulting connector is hard to read. The two features here get ~90% of the value at a fraction of the complexity and risk.

## Scope

**In:**
- A **peak-q cursor**: hovering near a detected peak draws a full-height vertical guide across the entire stack at that peak's q, and shows the q value in a small readout chip. Centralized at the chart level; reuses the shared x-axis projection so the guide aligns with both the beads and the axis ticks.
- A **log/linear q-axis toggle** (`SegmentedControl`) the user can flip live; default log, mirroring the Focus hero plot.

**Out (non-goals):**
- Migration track / cross-trace reflection connectors / ghost rings (deferred — see above).
- Per-bead q-halo highlighting of the same-q peaks in other rows (the guide line already shows alignment; lighting beads requires touching `TracePlot`).
- A free (non-snapped) floating crosshair anywhere on the canvas (v1 snaps to real detected peaks only — an honest guide marks a measured peak, not an arbitrary cursor position).
- Any change to `TracePlot` / `PlotFrame` (leaving per-row interaction off avoids PlotFrame's unconditional wheel `preventDefault`, which would trap page scroll over the chart).

## Feature 1 — Peak-q cursor (guide + readout)

### Behavior
- On pointer move over the stack, snap the pointer's x to the nearest detected peak (across all rows' anchors) within a horizontal tolerance (`PEAK_HIT_PX` = 10px). When a peak is within tolerance, set the cursor q to that peak's q; otherwise clear it.
- While a cursor q is set:
  - A **1px vertical guide** spans the full stack height (`TOTAL_H`) at `x = sharedX.to(cursorQ)`, in the chart's interaction-accent color. `pointer-events: none` so it never blocks hover.
  - A small **readout chip** shows `cursorQ.toFixed(3)` (mono, paper-on-hairline), anchored at the top of the guide.
  - The guide extends a short tick into the bottom axis strip so the reading ties to the q-axis.
- On pointer leave, clear the cursor (guide + chip disappear).

### Mechanism
- **Centralized, not per-row.** `WaterfallChart` attaches `onPointerMove` / `onPointerLeave` to the stack container. (Per-row `TracePlot` interaction stays off — see non-goals.)
- **Shared projection.** Build once: `const sharedX = makeAxis(qDomain, [4, plotW - 4], xType)` — identical to the bottom axis (`WaterfallChart.tsx:129`), so a row's bead at q, the axis tick at q, and the guide at q all land on the same container-relative pixel. (A row `TracePlot` with `axes={false}` uses 4px margins, making its internal q→px map equal to `[4, plotW-4]` in container coords — that is why the existing axis already uses this range.)
- **Snap helper (pure, unit-tested).** New `src/print/waterfall/cursor.ts`:
  ```ts
  import { hitTestPeaks } from "../plot/interaction";
  import type { WaterfallRow } from "./waterfallModel";

  /** All anchors across all rows as hit-test candidates (id+q suffice; same q → same px). */
  export function cursorCandidates(rows: WaterfallRow[]): { id: number; q: number }[] {
    return rows.flatMap((r) => r.anchors.map((a) => ({ id: a.id, q: a.q })));
  }

  /** Snap a container-relative px to the nearest peak q within tol; null if none. */
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
  (`hitTestPeaks` already ignores `id < 0` and resolves ties to the later peak — reused as-is.)
- **State, controlled-or-internal** (mirrors the existing `hoveredKey`/`internalHot` idiom in this same file): add `hoveredQ?: number` controlled prop; internal `internalHoveredQ` state; effective `cursorQ = hoveredQ ?? internalHoveredQ`. The pointer handler calls `setCursorQ`, which sets internal state only when `hoveredQ === undefined` and always fires the (now wired) `onHoverQ?.(q)` / `onHoverQ?.(undefined)`. Making it controllable lets tests assert guide rendering without fighting JSDOM's zeroed `getBoundingClientRect`.

### Rendering (placement-only — `WaterfallChart` is NOT design-guard-exempt)
- The guide line and readout chip are **positioned HTML `<div>`s with token classes only** (no SVG color literals, no `text-[…]`/`rounded-[…]`):
  - Guide line: `absolute top-0` 1px-wide div, `height: TOTAL_H`, `left: sharedX.to(cursorQ)`, token interaction-accent background class, `pointer-events-none`. `data-role="wf-qguide"`, `data-q={cursorQ}`.
  - Readout chip: `absolute` div near the guide top, `bg-plate border border-hair rounded-sm text-meta font-mono text-ink px-1`, `pointer-events-none`. `data-role="wf-qreadout"`, text = `cursorQ.toFixed(3)`.
- **Refactor-on-contact:** if the interaction-accent background lacks a token utility class (e.g. `bg-accent`), add it in the exempt token layer (`src/index.css` / the styles layer) rather than inlining a raw color — keep `WaterfallChart` placement-only. The implementer verifies the exact class names against the live token set before use.

## Feature 2 — Log/linear q-axis toggle

### Behavior
- A `SegmentedControl` with options **`log` / `linear`** sits in a thin toolbar row above the stack (right-aligned, layout-only flex). Flipping it re-renders the axis, the row traces, and the cursor's shared projection in the chosen scale. Default **log**.

### Mechanism
- Mirror the controlled-or-internal idiom again: keep `xType?: "log" | "linear"` as an optional **controlled** override; add `internalXType` state initialized to `"log"`; effective `scale = xType ?? internalXType`. Add `onXTypeChange?: (next) => void`. The `SegmentedControl.onChange` sets internal state when uncontrolled and always fires `onXTypeChange?`.
- Use `scale` everywhere `xType` is used today (the per-row `TracePlot xType`, the bottom `makeAxis`, and the new `sharedX`).
- Expose the active scale on the root for tests: `data-xtype={scale}`.
- `SegmentedControl<"log" | "linear">` props: `options=[{value:"log",label:"log"},{value:"linear",label:"linear"}]`, `value={scale}`, `onChange`, `aria-label="q axis scale"`, `size="sm"`, `testId="wf-scale"`.

## Files

- **Modify** `src/print/waterfall/WaterfallChart.tsx` — wire `onHoverQ`; add `hoveredQ` + `onXTypeChange` props; internal `internalHoveredQ` + `internalXType` state; stack pointer handlers; guide + readout overlay; toolbar with `SegmentedControl`; `data-xtype` on root.
- **Create** `src/print/waterfall/cursor.ts` — `cursorCandidates`, `snapToPeakQ` (pure).
- **Modify** `src/print/waterfall/index.ts` — export the cursor helpers if useful downstream (optional).
- **Create** `test/print-waterfall/cursor.test.ts` — unit-test `snapToPeakQ` (in tolerance → snaps to nearest q; out of tolerance → null; empty rows → null; ignores `id < 0`).
- **Modify** `test/print-waterfall/WaterfallChart.test.tsx` — add: controlled `hoveredQ` renders `wf-qguide` + `wf-qreadout` with the right `data-q` / text; omitting `hoveredQ` renders neither; `SegmentedControl` present, clicking `linear` sets `data-xtype="linear"` and fires `onXTypeChange("linear")`; default `data-xtype="log"`.
- **Modify** `src/print/waterfall/WaterfallChart.stories.tsx` — drop the fixed `xType` from at least one story so the live toggle is exercised; add a story with a controlled `hoveredQ` so the guide is visible statically in Storybook; keep Wide/Narrow/MixedStates.

## Testing approach

- **Pure first:** `snapToPeakQ` carries the hit-test logic and is fully unit-tested without the DOM.
- **Controlled-prop rendering:** the guide/readout are asserted via the controlled `hoveredQ` path (set the prop, query `data-role`), sidestepping JSDOM's zeroed `getBoundingClientRect`. One optional integration test may mock `getBoundingClientRect` to exercise the pointer→snap path end-to-end; if it proves flaky in JSDOM, rely on the pure helper test + controlled-prop test (documented, not silently dropped).
- **Toggle:** assert presence, default `data-xtype="log"`, click→`linear`, callback fires.
- Gates per task: `npm test -- print-waterfall` green, `npm run lint:design` clean (proves placement-only), `npx tsc --noEmit -p tsconfig.build.json` clean. After the batch: `npm run build` and `npm run build-storybook` exit 0.

## Design-guard compliance

`WaterfallChart` and its helpers live outside the exempt dirs, so: compose primitives (`SegmentedControl`) + layout/token Tailwind classes only; no inline appearance literals. The guide/readout use token classes; any missing token utility is added in the exempt layer as a reviewed refactor-on-contact. `lint:design` staying green is the proof.
