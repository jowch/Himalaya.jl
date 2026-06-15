# Greenfield Trace-Plot Engine — Design

**Date:** 2026-05-31
**Branch:** `worktree-greenfield-ui-rebuild`
**Status:** design converged; claims verified against live source (2026-05-31).
Plan #1 (engine + hero/mini proof) written:
`docs/superpowers/plans/2026-05-31-greenfield-trace-plot-engine.md`.

## Goal

Replace the Observable-Plot–based plotting layer with an **owned, declarative,
d3-backed trace-plot engine** rendered as React/SVG, aligned with the greenfield
design system (closed-look, DESIGN.md plot invariants). One engine renders every
trace plot in the app at different feature levels.

## Scope

**In scope — all are "traces", one `<TracePlot>` engine, feature-gated:**

- **Focus hero** — a single measured trace + σ band + rich peak glyphs + axes +
  wheel-zoom + q-link cursor + peak add/remove/exclude.
- **Multi / waterfall** — N stacked traces in bands (the Series overlay).
- **Mini** — folio mini-waterfall + scoping sparkline (axes and interaction off,
  cheap enough for a masonry of many cards).

**Not a `<TracePlot>`, but in the Observable-Plot *removal* scope:**

- `MillerPlot` — **dropped, not ported** (user, 2026-05-31): the ratio-vs-q
  regression scatter has been superseded by the residual view. It is **already
  orphaned dead code** — the IndicesCard mount was removed in R4 (verified: no
  JSX mount or import anywhere; only the definition + a stale test remain). So
  this is a plain delete of `MillerPlot.tsx` + `test/MillerPlot.test.tsx`, not
  the removal of a live inset.
- **Figure export** — **four** files import Observable Plot:
  `lib/figure-export/marks/traceExportMarks.ts`, `marks/multiTraceExportMarks.ts`,
  `renderer.ts`, **and `types.ts`** (uses Plot's `Markish` types). All four must
  migrate to serialize the owned `<TracePlot>` SVG; `types.ts` is easy to miss
  and would block the final dependency drop. Must land before Plot can be
  removed from the dependency tree.

**Confirmed standalone, already SVG, NOT part of this engine:**

- `CombPanel` — verified raw SVG (`<svg data-testid="comb-svg">`, no Plot import),
  owns **both** the comb and the **residual** views via its SegmentedControl. The
  residual plot stays SVG. Rides the shared spine (phaseColor, q-link tolerance,
  predicted-q) but is not a `<TracePlot>`. **Future:** extract both the comb and
  residual panels into `src/print/plot/` alongside the trace engine and reskin —
  a separate plan, not this one.

## Background — why own it

`PlotSurface` (`src/components/PlotSurface.tsx`, read in full) is original-design
lineage. Its root problem is **inverted projection ownership**: Observable Plot
owns q→pixel, and our interactive overlay reverse-engineers the projection back
out of Plot's rendered DOM (`plotEl.scale("x")` / `applyQHelper`). This readback
happens in **three** places — `makeCtx` (`:203–218`), the `handleClick` x-scale
read (`:223–226`), and the `overlayContent` IIFE (`:366–389`) — all routed
through `applyQHelper`/`invertQHelper` (`lib/plot/invertQ.ts:38–57`). That single
inversion forces everything that makes it brittle:

1. **A ~420-line imperative bridge** — `Plot.plot()` → `replaceChildren` → manual
   `addEventListener` (click/wheel/dblclick/brush) → manual teardown, in a
   20-dependency `useCallback`, plus a `_resizeKey`/`ready` re-render dance to
   repaint the overlay after Plot builds. This fights React.
2. **Split appearance** — fonts/colours/ticks injected via Plot's `style`/`x`/`y`
   config, then *supplemented* by our own `labelDodge.ts` (142 lines written to
   fix Plot's tick labels). Closed-look cannot own the axes here.
3. **A leaked abstraction** — escape hatches (`onPointerClick` returning a bool to
   preempt, `onPointerMove`, `onBrush`) because MultiTracePlot's per-band hit-test
   didn't fit the "one peak by id" contract; `MillerPlot` excluded outright. It is
   not a universal spine — it is the Focus/Compare q-plot host with bolt-ons.

Observable Plot buys us only: log/linear **scales**, **axis-tick generation**, and
**`linearRegressionY`** (Miller only) — all of which are `d3-scale`/`d3-array`,
which Plot itself depends on (`d3 ^7.9`). So the migration is **subtractive**: keep
the d3 math Plot sits on, drop Plot's rendering host.

## Architecture

### Own the projection

A pure factory, no DOM:

```
projection(domain: [number, number], range: [number, number],
           type: "log" | "linear") -> {
  toX(q: number): number;        // d3-scale
  invertX(px: number): number;   // d3-scale .invert()
  toY(v: number): number;
  ticks(count?: number): number[]; // d3-scale .ticks()
}
```

Backed by `d3-scale` (`scaleLog` / `scaleLinear`). Once *we* own q→px, the entire
plot — data line, σ band, peak glyphs, axes, cursor — is **one declarative SVG
tree** driven by `(samples, peaks, projection)`. No Plot instance, no scale
readback, no imperative event bind (React `onWheel`/`onClick` with `invertX`), no
`_resizeKey`/`ready` dance (a `ResizeObserver`→width state, or pure CSS + viewBox).

### Dependencies

Make **direct** (all already present transitively via Plot — zero install weight):
`d3-scale`, `d3-shape` (line/area `d` strings), `d3-array`, `d3-format`. `visx` is
the documented fallback if axis/grid JSX becomes tedious (it wraps the same d3);
default is direct d3-modules for total closed-look control over ~5 mark types.

### Keep (port from PlotSurface as pure functions — the trusted *behaviours*)

- Log-aware wheel-zoom-about-cursor math (`PlotSurface.tsx:258–284`; log branch
  clamps cursor q to `≥ 1e-6`, min-span guard `1e-4`, `factor = exp(deltaY·0.001)`).
- `hitTestPeaks` (lines 104–123) — pure, keep verbatim (skips `id < 0` and
  non-finite px; `<=` so a later equidistant peak wins).
- dblclick reset (`:290–292`, `onReset?.() ?? onXDomain(null)`).
- The hit radius is **pixel-based**: `PEAK_HIT_PX = 10` (`:26`), plus a `4`px
  brush dead-zone (`:328`). There is **no** relative-q "q-link tolerance"
  constant — the "#180 q-link feel" comment annotates the pixel radius. Keep the
  pixel model verbatim.
- The `overlay(ctx)` API *shape* — but now **we** hand consumers the projection;
  no readback.

### Drop (the machinery)

`Plot.plot()` data layer; the `renderPlot`/`replaceChildren`/`addEventListener`
bridge; scale-readback into Plot internals; the two-SVG (Plot + absolute overlay)
stack → collapse to one SVG; the `_resizeKey`/`ready` dance.

### The shared spine it rides

(see `project_greenfield_ui_phase2_sourcing` memory)

- `phaseColor()` — colour = phase, always; never hue for state.
- `effective_peaks` + exclusion contract — the read-only renderers honour the
  same curated peak set the hero shows (fixes today's sparkline-shows-raw drift).
- Canonical peak vocabulary = greenfield `print/ui/PeakGlyph` (retire the legacy
  `components/ui/PeakGlyph` surrounding-ring version).
- DESIGN.md plot invariants: Second-Channel (glyph silhouette carries provenance,
  not hue), Fixed-Scale, Mono-Means-Measurement (mono tick labels), log-x default.

### Data model

`TraceModel` = `{ samples: {q,I}[]; peaks: Peak[]; phase: string|null; state?: AssignmentState }`.
Feeders stay plural and live in `lib`:

- **measured** — file-backed `Trace` (`{q,I,sigma}`) for the hero + sparkline.
- **modelled** — `buildMiniWaterfall` renders a synthetic Gaussian model
  positioned at `effective_peaks` q-values, with per-peak amplitudes derived from
  measured peak `intensity` (fallback `0.5`); it never reads the raw `Trace.I`
  integration curve. Cheap per card; for the folio waterfall.

The renderer takes a `TraceModel` (or N of them) + a `projection`; the feeder
choice is a data concern, not a renderer fork.

## File layout

The plot layer is its own subsystem with its own internal primitives — a sibling
of `print/ui/`, not folded into it:

```
src/print/plot/
  projection.ts        # d3-scale wrapper: toX/invertX/toY/ticks (pure, unit-tested)
  interaction.ts       # ported pure fns: zoom-about-cursor, hitTestPeaks, q-link tol
  TracePlot.tsx        # the engine / public component (feature-gated)
  PlotFrame.tsx        # SVG root + ResizeObserver width + gesture surface
  Axis.tsx             # x/y axis as JSX over scale.ticks() (Mono-Means-Measurement)
  marks/
    TraceLine.tsx      # line + σ-band paths (d3-shape)
    PlotPeaks.tsx      # peak marks — composes print/ui/PeakGlyph
  index.ts             # barrel
  # FUTURE (separate plan): CombPanel.tsx, ResidualPanel.tsx extracted here
```

Stories colocate (`src/print/plot/*.stories.tsx`, auto-discovered by the
`../src/print/**` glob); unit tests live under `test/print-plot/` per the house
convention (tests are NOT colocated — e.g. `test/print-ui/`). **Design guard:**
add `print/plot/` to the `isExcluded(relPath)` directory exemption in
`scripts/check-design.mjs` (`:134–137`) — the same full exemption `print/ui/`
gets, because the plot *authors* SVG appearance + phase colour. (The
detector/heatmap/figure-export precedent uses the weaker colour-only allowlist,
which would still ban the engine's axis fonts — `isExcluded` is the right lever
for a `print/` appearance subsystem.) Consumer-facing `className` on
`<TracePlot>` stays placement-only.

## Component surface (sketch, to be firmed in the plan)

```
<TracePlot
  traces={TraceModel[]}          // 1 = hero/sparkline, N = waterfall
  projection={Projection}        // or {domain,type} and TracePlot builds it
  layout="single" | "bands"
  axes={false | {x?,y?}}         // off for mini
  interaction={false | {...}}    // off for mini; zoom/q-link/peak-edit for hero
  overlay={(ctx) => ReactNode}   // comb/rings/cursor; ctx carries OUR projection
  size?  className?              // placement only
/>
```

Mini = `axes={false} interaction={false}`; hero = all on. Same component.

## Migration (strangler-fig)

1. Build `<TracePlot>` + `projection` in `src/print/`, honed in Storybook against
   the captured real fixtures (`src/print/fixtures/`).
2. Prove on the **Focus hero first** (hardest: scales + axes + zoom + q-link +
   peak edit). If it carries the hero, mini/waterfall fall out.
3. Migrate consumers one at a time; keep `PlotSurface` alive until the last is
   moved, then retire it.
4. **Observable Plot leaves the tree entirely.** Verified removal surface:
   `PlotSurface`, `TraceViewer`, `MultiTracePlot` (+ `MemberTraceLayer`,
   `MemberHeatmapLayer`, `CrossTraceTrackingLayer`), `peakMark.ts` (both copies —
   `print/ui/peakMark.ts` is Storybook-only), `MillerPlot` (already orphaned —
   plain delete + its test), and the **four** figure-export Plot files
   (`marks/traceExportMarks.ts`, `marks/multiTraceExportMarks.ts`, `renderer.ts`,
   `types.ts`). After the trace engine, the export migration, and the Miller
   removal land, drop `@observablehq/plot` from `package.json` — and add
   `d3-scale`/`d3-shape` as direct deps first (they are transitive via Plot today,
   so removing Plot without promoting them breaks resolution).

## Testing

- **Unit:** `projection` (toX/invertX round-trip, log + linear, ticks), path `d`
  generation, peak hit-test, zoom-about-cursor math. Pure functions, fast.
- **Storybook:** stories per feature level against real fixtures (hero on
  `traces/65.json` + snapshot peaks; waterfall on `transitionSeries`; sparkline on
  a measured trace; empty/form-factor states).
- **Behaviour parity:** q-link feel, zoom clamp, dblclick reset match the
  incumbent (the behaviours are ported verbatim, so this is a guard not a rewrite).
- **Design guard / impeccable:** appearance lives in the engine; consumer
  `className` is placement-only; `npx impeccable detect` clean.

## Resolved (spec review, 2026-05-31)

1. **MillerPlot** — dropped, not ported (superseded by the residual view). ✓
2. **CombPanel + residual** — standalone, already SVG, confirmed; not folded in. ✓
3. **`d3-scale` direct vs `visx`** — **direct d3-modules**; visx evaluated and
   rejected (not needed for ~5 marks). ✓

## Open (defer to the plan)

- **Decimation threshold:** at what point count do we down-sample the line path
  for the hero's drag-zoom (visual resolution ≈ container width in px)?
- **Export sequencing:** migrate figure export in the same plan as the hero, or as
  a fast-follow once `<TracePlot>` SVG is stable?
