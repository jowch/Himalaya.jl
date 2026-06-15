# Greenfield WaterfallChart — design

> Part of the greenfield "The Print" rebuild (branch `worktree-greenfield-ui-rebuild`, NOT merged).
> Renderer track, after the comb + residual batch. Source mockup:
> `docs/redesign-mockups/2026-05-29-series-plot.html` (the `#waterfall` hero).

## What this is

The **Series** surface's hero figure: a **stacked-offset waterfall** of a series'
member traces, ordered low→high along the recipe variable. Each member's measured
SAXS trace gets its own vertical slot; phase-coloured lines, with light
peak-anchor beads at the member's indexed Bragg peaks. The job of the surface is
*"what changes across the variable?"* — the eye follows a reflection marching in q
and growing/shrinking across the stack.

## Scope

**In:**
- `WaterfallChart` — the stacked-offset multi-member renderer.
- A small refactor of the existing `TracePlot` / `TraceLine` to serve as the
  per-member plot unit (see "TracePlot refactor").

**Deferred, with reasons (NOT in this batch):**
- **Migration track** (hover an anchor → dashed connector of the same reflection
  across the stack + ghost rings for predicted-but-absent). Blocked on data the
  series member snapshot doesn't carry: `MemberSnapshotPeak` has no
  `ratio_position`/phase, and `confirmed_index.peak_ids` is observed-only with no
  slot for a predicted-but-absent reflection — so the cross-stack connector would
  misalign when a reflection drops out, and the "ghost ring" has no data source.
  That's a snapshot-shape (backend) decision, not a renderer concern. The anchors
  this batch *does* build are the interaction handles a future track hangs off.
- **Heatmap toggle.** YAGNI for the typical 6–8-exposure titration (stacked traces
  are strictly more legible at that density); its continuous intensity-colour LUT
  competes with the "colour = meaning / phase colour" rule. A separate renderer
  when a dense series actually needs it.
- **Per-peak phase colouring under coexistence.** Same per-peak phase-identity gap
  as the track — anchors colour by the member's dominant phase for now.
- **Export clean-figure**, rail panels (Members/Reading/Export), the `#phasestrip`
  companion (already a built primitive — see below). These are tier-2 plate / other
  renderers, planned per-batch.

## Key findings that shape the design

1. **TracePlot is single-trace in practice.** Every greenfield `TracePlot` usage
   passes a 1-element `traces={[model]}`; the only consumers are its own stories.
   The multi-trace array is unused generality from the legacy `MultiTracePlot`.
2. **Composing N single-trace plots replaces the normalization library.** The
   legacy `computeYBands`/`computeReference`/`applyNormalization` machinery
   (`src/lib/comparison/`) exists because `MultiTracePlot` drew every member in
   ONE coordinate space and had to mathematically partition + clip. N separate
   per-member plot boxes give that partitioning **structurally** — each box
   auto-scales its single trace. The waterfall needs none of that math.
3. **The vertical phase-strip companion already exists.** `PhaseStrip` (a `ui/`
   primitive) supports `orientation="vertical"`, documented as "the Series
   waterfall companion." The waterfall aligns to it; it does not rebuild it.
4. **Real data is the truth test.** Build/story/test against
   `src/print/fixtures/realSeriesMembers.ts` (the real Sample-9
   `Ia3d → Im3m+Ia3d → Im3m` transition; measured traces in
   `src/print/fixtures/traces/<exposureId>.json`), NOT the mockup's idealized
   synthetic curves.

## TracePlot refactor

`TracePlot` and `TraceLine` (`src/print/plot/`) are design-guard-exempt and only
consumed by their own stories today, so this is a safe, contained change.

- **Single trace.** `TracePlotProps.traces: TraceModel[]` → `trace: TraceModel`.
  Remove the `traces.flatMap(...)` over peaks/extent; operate on the one trace.
  Update `TracePlot.stories.tsx` (two stories: `traces={[heroModel]}` →
  `trace={heroModel}`, `traces={[annotatedModel]}` → `trace={annotatedModel}`).
- **Phase-coloured line.** Add `color?: string` to `TraceLineProps` (default keeps
  the current `var(--color-ink-soft)`). `TracePlot` passes
  `trace.phase ? phaseColor(trace.phase) : "var(--color-ink-faint)"` so the line
  carries phase meaning (previously the line was always neutral; only peaks were
  coloured). This is a latent-gap fix, not waterfall-specific.
- **Top headroom.** Pad the y-domain at the top (a small fraction) so a tall peak
  doesn't slam the plot ceiling — approximates the legacy 70% working-band
  headroom. Applied in `makeProjection` inputs (extend `yExtent` top by a fixed
  factor, e.g. `×1.08`), gated so it doesn't change axis tick labels unexpectedly;
  the simplest form is a `yHeadroom?: number` prop (default `0`, waterfall passes
  a small value). Headroom must not break the existing peak-position / q-readout
  tests.
- Everything else unchanged: `axes`, `xDomain`, `xType`/`yType`, `interaction`,
  `peaks` (PlotPeaks layer + hover + q-readout), `show`, `highlightPeakIds`,
  `overlay`, `width`/`height`, `paperColor`.

## WaterfallChart

### Files

- `src/print/waterfall/waterfallModel.ts` — pure view-model + builder.
- `src/print/waterfall/WaterfallChart.tsx` — the renderer (composition).
- `src/print/waterfall/waterfall.fixtures.ts` — plain module (NOT a stories file)
  deriving `WaterfallRow[]` from `realSeriesMembers` + the trace JSON.
- `src/print/waterfall/WaterfallChart.stories.tsx` — stories.
- `src/print/waterfall/index.ts` — barrel.
- `test/print-waterfall/waterfallModel.test.ts`,
  `test/print-waterfall/WaterfallChart.test.tsx` — unit tests.

### View-model (`waterfallModel.ts`)

```ts
import type { SeriesMember, Trace } from "../../api";

/** One indexed-peak bead position on a member's trace. */
export interface WaterfallAnchor {
  /** The effective-peak id (from peak_ids) — maps 1:1 to PlotPeak.id so the
   *  peaks layer's hover / q-readout key cleanly; ids are unique per member. */
  id: number;
  q: number;
  /** Member dominant phase (drives bead colour); null = unindexed (neutral). */
  phase: string | null;
}

/** One member row, fully resolved for rendering (low→high variable order). */
export interface WaterfallRow {
  /** Stable key — the member id (or exposure id) as a string. */
  key: string;
  /** Right-margin label: label_override ?? the recipe variable value. */
  label: string;
  /** Measured trace (q[]/I[]); empty trace renders no line (still reserves a slot). */
  trace: Trace;
  /** Dominant confirmed phase → line + anchor colour; null = neutral. */
  phase: string | null;
  /** Durable 3-state assignment; "form_factor"/"null" suppress anchors. */
  state?: "indexed" | "form_factor" | "null";
  /** Indexed-peak beads — confirmed_index.peak_ids ∩ effective_peaks (by id). */
  anchors: WaterfallAnchor[];
  /** Relative slot height (durable band_height; default 1). */
  bandHeight: number;
  /** Vertical nudge in slot-fractions (durable y_offset; default 0). */
  yOffset: number;
}

/** Map API members + a trace lookup → rows in display order (low→high). */
export function toWaterfallRows(
  members: SeriesMember[],
  tracesById: Record<number, Trace>,
): WaterfallRow[];

/** Positive q-extent across every row's trace, padded ×0.97 / ×1.03. */
export function waterfallQDomain(rows: WaterfallRow[]): [number, number];
```

`toWaterfallRows` rules:
- Sort by `display_order` ascending (caller renders bottom-up; see layout).
- `phase`: the member's dominant confirmed phase. Derive from
  `snapshot.confirmed_phases` (first/dominant) or `snapshot.confirmed_index.phase`;
  `null` when `assignment_state` is `null`/`form_factor` or there is no index.
- `state`: from `snapshot.assignment_state` (missing → `"indexed"`).
- `anchors`: for each id in `snapshot.confirmed_index.peak_ids`, look up the
  matching `snapshot.effective_peaks` entry by `id`; emit `{ id, q, phase }` with
  the member's dominant phase. Empty for `form_factor`/`null` members or no index.
- `label`: `label_override ?? <variable value>` — for the fixture, the captured
  per-member label (e.g. `"C06-5"`); a series-level variable string is the tier-2
  plate's concern, so the row label uses `label_override` and falls back to the
  exposure id.
- `bandHeight`: `member.band_height ?? 1`; `yOffset`: `member.y_offset ?? 0`.
- `trace`: `tracesById[member.exposure_id] ?? { q: [], I: [] }`.

### Renderer (`WaterfallChart.tsx`)

Props:

```ts
interface WaterfallChartProps {
  rows: WaterfallRow[];          // low→high (display order); rendered bottom-up
  xType?: "log" | "linear";      // default "log"
  hoveredKey?: string;           // controlled row emphasis (hot)
  onHoverRow?: (key?: string) => void;
  onHoverQ?: (q?: number) => void;
  maxWidth?: number;             // fit-to-width ceiling (default 1080, mockup plate)
  minWidth?: number;             // floor below which it stops shrinking (default 560)
  className?: string;            // PLACEMENT ONLY
}
```

Layout & rendering:
- **Fit-to-width:** the figure fills its container up to `maxWidth`, down to
  `minWidth`. Use a single responsive SVG-per-row stack inside one width-measured
  container (or a viewBox that scales) — no horizontal scroll (unlike the comb).
- **Bottom-up stack:** `display_order` 0 sits at the **bottom** (matches the mockup
  hero and the bottom-up vertical `PhaseStrip`). Reverse iteration for vertical
  placement; total vertical space split by `bandHeight` ratios, each row nudged by
  `yOffset`.
- **Per-member box = one refactored `TracePlot`:** `trace`, `axes={false}`, shared
  `xDomain` = `waterfallQDomain(rows)`, `xType`, `color` via phase, `peaks` = the
  row's `anchors` mapped to the PlotPeaks shape (so beads + anchor-hover + the
  q-readout come from the existing peaks layer — no bespoke anchor drawing), a
  small `yHeadroom`. `band={false}` (no ±σ in a dense stack).
- **Shared bottom axis:** one `Axis` (`makeAxis(qDomain, [..], xType)`,
  `orientation="bottom"`, label `"q (Å⁻¹)"`) beneath the stack — the only axis;
  member boxes are axis-less.
- **Sample label** per row at the right margin (`data-role="wf-label"`, mono,
  `text-ink-soft`; emphasized when hot).
- **Row hover (hot/dim):** the wrapper tracks `hoveredKey` (controlled) /
  internal hover. The hot row's box renders at full opacity + heavier weight; the
  others dim (reduced opacity). Implemented as a `data-hot` / `data-dim` attribute
  on each row group with token-driven opacity — NOT inline appearance literals.
  `onHoverRow` fires on row enter/leave; `onHoverQ` bubbles from the member
  TracePlot's peak hover.
- **Empty / form-factor rows:** a `form_factor`/`null`/empty-trace member still
  reserves its slot and shows its label, with no anchors and (for empty trace) no
  line — honest about "a member with no Bragg peaks."

### Design-guard

`WaterfallChart` aims to be **placement-only** (composition of the exempt plot
layer + Tailwind layout/token classes), needing **no** new design-guard exemption.
Appearance (line colour, bead glyphs, axis strokes) lives inside the already-exempt
`print/plot/`. If a genuinely unavoidable raw SVG appearance literal surfaces in
`WaterfallChart` itself, add `print/waterfall/` as the 6th `isExcluded` prefix in
`scripts/check-design.mjs` (+ a `test/check-design.test.ts` block, mirroring the
`print/comb/` exemption) — but only then, and as a reviewed step.

## Alignment with `PhaseStrip`

The vertical `PhaseStrip` companion is already built and is the tier-2 Series
plate's responsibility to place beside the waterfall. This renderer does not embed
it. The contract: both stack in the same low→high order; the plate sizes the strip
to the waterfall's plot height. (The strip currently lays segments top→bottom; the
plate passes segments in the matching order. Cross-component pixel alignment is the
plate's concern, out of scope here.)

## Testing

Unit (`test/print-waterfall/`), asserting via `data-role`/`data-*`/text (never
class strings; JSDOM computes no layout):

`waterfallModel.test.ts`
- `toWaterfallRows` preserves display order and derives `phase` from
  `confirmed_phases`/`confirmed_index`.
- anchors = `confirmed_index.peak_ids` resolved to their `effective_peaks` q;
  count matches; `form_factor`/`null` members yield zero anchors.
- a member with no matching trace yields an empty trace (no throw).
- `waterfallQDomain` returns the padded positive q-extent across all traces.

`WaterfallChart.test.tsx`
- renders one row group per member (`data-role="wf-row"`), bottom-up: the
  display-order-0 row sits at the largest `y` (lowest on screen).
- each row renders a phase-coloured trace line (the refactored TracePlot line
  carries the phase colour, not ink-soft) and a bead per anchor.
- one shared bottom axis (`data-role="axis-bottom"`); member boxes are axis-less.
- hovering a row marks it `data-hot` and dims siblings (`data-dim`); fires
  `onHoverRow` with the key, and `undefined` on leave.
- hovering an anchor fires `onHoverQ` with the bead's q (via the peaks layer).
- `xType="linear"` is honoured (bead x-positions differ from log for the same q).

TracePlot refactor regression: existing `test/print-plot/*` stay green after
`traces`→`trace`; add a test that `TraceLine` honours a `color` prop and defaults
to ink-soft.

Gates per the established cadence: `npm test -- print-waterfall` and
`npm test -- print-plot` green; `npm run lint:design` clean (proves placement-only);
`npx tsc --noEmit -p tsconfig.build.json` clean. After the batch: `npm run build`
and `npm run build-storybook` exit 0. Manual fidelity: `npm run storybook` →
`waterfall/WaterfallChart` vs the `#waterfall` region of the series-plot mockup.
