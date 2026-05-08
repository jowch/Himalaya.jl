# Figure Export — Copy / Save buttons for plots

**Date:** 2026-05-08
**Status:** Draft

## Context

Two of HimalayaUI's plot surfaces — TraceViewer (Index page) and MultiTracePlot (Compare page) — are the figures users want to drop into progress reports and working presentations. Today there is no way to extract them: the only path is OS-level screenshot, which captures the surrounding chrome and bakes in the dark theme.

This spec proposes a Copy / Save affordance per plot that produces a "presentation ready" version of the figure suitable for pasting into a slide deck or attaching to a report.

**Out of scope:** publication-quality output (PDF, embedded LaTeX fonts, vector typography). For final figures, users are expected to render outside Himalaya.

## Goals

- One-click PNG → clipboard, paste-ready, from each supported plot.
- One-click download with format choice (PNG / SVG).
- Output that reads on its own: light background, larger fonts than the on-screen plot, embedded title and legend.
- No new runtime dependencies — reuse Observable Plot, which is already shipped.

## Non-goals

- Publication-quality export (PDF, embedded glyphs, LaTeX fonts).
- Visual regression / pixel-identical reproducibility against the on-screen render.
- Server-side rendering (Makie / Plots.jl).
- Export for the MillerPlot inset and the DetectorImage canvas. Both are deferred — neither is the typical "figure for a report" candidate.
- R² annotation in the trace export legend (deferred).

## Scope

In:
- TraceViewer (Index page)
- MultiTracePlot (Compare page)

Out (deferred):
- MillerPlot inset
- DetectorImage canvas

## Approach: headless re-render via Observable Plot

Don't snapshot the on-screen plot. When the user clicks Copy or Download, build a fresh `Plot.plot()` call with export-tuned parameters — forced light theme, larger fonts, fixed dimensions, inlined title and legend, wider strokes — and use that as the export source. PNG is produced by drawing the SVG into an offscreen canvas at 2× DPI. SVG is produced by serializing the same node.

Rationale:
- The "presentation ready" requirements (light theme override, larger fonts, fixed size, embedded title and legend) all favour an explicit re-render over a clone-and-mutate path. Clone paths require post-processing the cloned tree to scale strokes and fonts, which is brittle.
- TraceViewer's interactive overlay (peak triangles, predicted-q ticks, cursor) is a sibling SVG positioned on top of the Plot SVG, not a child of it. A naive clone of the Plot SVG loses the overlay; merging the two SVGs by hand is error-prone (margin and coordinate-space reconciliation).
- Server-side rendering would add a network round-trip per Copy click, kill the snappy paste-ready UX, and require keeping two render libraries in sync.

**Tradeoff (accepted):** the export-only mark builders duplicate logic from on-screen renderers (peak triangles, predicted-q tick geometries). When peak-source types or tick geometry change, both paths need updates. This is the price of leaving TraceViewer's load-bearing interactive overlay (peak hit-testing, hover, cursor follow) untouched. The export does not need pixel-identical parity with the screen — only visual parallelism.

## Per-plot figure design

The two plots have different purposes; their exports embed that difference rather than mapping a single generic shape onto both.

### TraceViewer export — "proposed index assignment for this exposure"

Self-explanatory contents:

- **Title:** `experiment · sample · exposure`
- **Trace** (q vs I, log/linear honouring the user's current `xType`, current `xDomain` / `yDomain`)
- **Peak markers** (auto / manual / excluded — distinguishable in the legend)
- **Predicted-q ticks** for the active group's indices, phase-coloured, with each phase labelled inline. Hovered-but-not-active candidates are excluded — hover state is transient and would not be present by the time the thunk runs.
- **Legend row at bottom:** peak-source key + the phases of the active group's indices

R² for active phases — deferred.

### MultiTracePlot export — "trace overlay across these members"

Self-explanatory contents:

- **Title:** comparison name (and parent experiment)
- **Stacked traces** in current band layout, current `xDomain`
- **Per-member labels** rendered inline at each band's vertical position, inside the SVG. Today these labels live in `MemberMetaGutter`, a sibling component to the plot; the export folds them into the figure so it reads alone.
- **Honours current toggles** (`showPeakTicks`, `showPeakLabels`) — these already exist on the page and are the user's "what to include" controls; respecting them means no separate export config
- **Legend:** minimal — only when grouping mode is `"sample"` (then sample → colour); in `"member"` mode each label is its own legend entry

## Architecture

```
src/lib/figure-export/
  types.ts                    — ExportSpec, LegendSpec
  renderer.ts                 — buildExportSvg, buildExportPng
  clipboard.ts                — copyPngToClipboard
  download.ts                 — downloadBlob
  filename.ts                 — slugifyForFilename, buildFilename
  presets.ts                  — EXPORT_DIMS, fonts, strokes, LIGHT_PALETTE
  marks/
    traceExportMarks.ts       — emits trace + peaks + predicted-q ticks for export
    multiTraceExportMarks.ts  — emits stacked traces + peak ticks/labels for export
  adapters/
    traceAdapter.ts           — TraceViewer state → ExportSpec
    multiTraceAdapter.ts      — MultiTracePlot state → ExportSpec

src/components/
  FigureExportControls.tsx    — Copy + split-Download buttons
```

### `ExportSpec` shape

```ts
interface ExportSpec {
  title: { primary: string; secondary?: string };
  width: number;
  height: number;
  plot: PlotConfig;       // marks + x/y configs for Plot.plot()
  legend?: LegendSpec;    // structured legend rows the renderer paints below
}
```

The adapter pre-bakes the figure design — title text, dims, marks, palette — and the renderer is dumb: lay out the spec, produce SVG or PNG.

### Adapters

Each adapter is a pure function from live plot state to `ExportSpec`. They speak the same `ExportSpec` to the renderer but their internals look completely different — the trace adapter threads peaks and active-group indices into mark builders that emit phase-coloured ticks with labels; the multi-trace adapter threads members into stacked-band marks with inline member labels.

```ts
function buildTraceExportSpec(args: {
  trace: Trace;
  peaks: Peak[];
  activeGroupIndices: IndexEntry[];
  experimentName: string;
  sampleName: string;
  exposureLabel: string;
  xDomain: [number, number] | null;
  yDomain: [number, number] | null;
  xType: "log" | "linear";
  qUnits: string;
}): ExportSpec
```

### Renderer

```ts
buildExportSvg(spec: ExportSpec): SVGElement
  // 1. Plot.plot(spec.plot) → inner plot SVG
  // 2. Wrap in outer SVG with explicit xmlns: white-bg <rect>, title <text>
  //    nodes, legend rows beneath
  // 3. Marks were built with literal palette colours at adapter time, so no
  //    var(--color-*) resolution is needed at render time

buildExportPng(spec: ExportSpec, scale = 2): Promise<Blob>
  // Serialize the SVG, draw onto an OffscreenCanvas at width*scale × height*scale,
  // canvas.convertToBlob({ type: "image/png" })
```

### Colour resolution

The export-only mark factories accept literal colours from `LIGHT_PALETTE` (in `presets.ts`) directly — they don't read CSS variables. The rendered SVG never contains `var(--color-*)`, so the renderer doesn't need a string-replace post-pass.

### Button component

```tsx
<FigureExportControls
  spec={() => buildTraceExportSpec({ /* ... */ })}
  filenameStem="himalaya-trace-{experiment}-{sample}-{exposure}"
/>
```

`spec` is a thunk evaluated at click time so it captures fresh state — `xDomain`, current peaks, active group, etc. — rather than stale state captured at component mount.

The Download button is a split:
- Main click → download PNG (the most common case)
- Chevron → focus-trapped popover with "Download as PNG" / "Download as SVG", pattern matching `ForksPopover`

### Placement

- **PlotCard title strip** (Index page): append after the q-range reset button, separated by a thin divider. The strip already holds fit-features + x-scale + q-range; placing export controls *next to* the controls that change what gets exported keeps the adjacency right.
- **ComparePage review header** (`compare-review-header`): append at the right end of the toggles row, after `ForksPopover`.

### Data flow at click time

```
user clicks Copy
  → FigureExportControls invokes spec()
  → adapter pulls live state via the closure
  → renderer.buildExportPng(spec) → Blob
  → clipboard.copyPngToClipboard(blob)
  → toast.success("Copied figure to clipboard")
```

No state lifted to Zustand or TanStack Query — the export is a pure derivation of state already present in the parent component when the user clicks.

## Browser compatibility

### Clipboard

`navigator.clipboard.write([new ClipboardItem({ "image/png": blob })])`. Three failure modes handled:

- **Insecure origin** — `navigator.clipboard` undefined. Pre-flight check; if missing, Copy is disabled with a tooltip "Clipboard requires HTTPS."
- **`ClipboardItem` unavailable** — same path; disable Copy, leave Download.
- **Runtime denial** (user-interaction policy, browser block) — caught and surfaced as `toast.error("Couldn't copy. Try Download instead.")`.

### SVG → PNG canvas pipeline

- Use `font-family: ui-sans-serif, system-ui, sans-serif` for export. Avoids `@font-face url(...)` references → no canvas-taint risk. Won't match the app's Plus Jakarta Sans exactly; accepted as the "progress report" tradeoff.
- Outer wrapper SVG must declare `xmlns="http://www.w3.org/2000/svg"`. Plot's inner SVG already does; the outer wrapper sets it explicitly.
- Standard recipe: serialize SVG → blob URL → `Image` → `OffscreenCanvas.drawImage` → `convertToBlob({ type: "image/png" })`.

## Disabled states

- PlotCard: Copy + Download disabled if `!traceQ.data || !peaksQ.data` (matches the existing `canFit` gate pattern).
- ComparePage: disabled if `members.length === 0 || traces.size === 0`.

## Filenames

- TraceViewer: `himalaya-trace-{experiment}-{sample}-{exposure}-{YYYY-MM-DD}.{ext}`
- Compare: `himalaya-comparison-{experiment}-{name}-{YYYY-MM-DD}.{ext}`

`slugifyForFilename(s)`: lowercase, non-alphanumeric runs → `-`, collapse repeated dashes, trim. Empty / all-special-char input falls back to a sentinel.

## Testing

Three layers, matching the existing project conventions.

### Vitest (unit)

- `traceAdapter.test.ts` — given fixture state, returns expected `ExportSpec` shape (title parts, dims, presence of tick marks per index, legend rows).
- `multiTraceAdapter.test.ts` — same for Compare.
- `renderer.test.ts` — `buildExportSvg(spec)` produces an SVG with no `var(--…)` literals (string regex), title `<text>` node present at expected position, legend group at expected y, white background `<rect>`.
- `filename.test.ts` — slugify covers spaces, slashes, parens, unicode, empty input.
- `clipboard.test.ts` — mock `navigator.clipboard.write` via `Object.defineProperty` (Clipboard is read-only on `navigator`); assert the expected `ClipboardItem` was passed.

### JSDOM constraint

`OffscreenCanvas.convertToBlob()` and `Image.decode()` aren't real in JSDOM. The PNG path is smoke-testable only (assert the function exists and returns a Promise). The SVG path is fully testable in JSDOM.

### Playwright E2E (mocked)

Two specs in `e2e/`:
- "Copy figure → clipboard contains PNG" — `context.grantPermissions(['clipboard-read', 'clipboard-write'])`, click Copy, read clipboard via `evaluate`, assert a PNG data URL.
- "Download → PNG / SVG" — `page.waitForEvent('download')`, assert the suggested filename matches the pattern and the file is non-empty.

### Visual regression

Out of scope. The "good enough for progress reports" framing means structural snapshot tests cover drift; pixel-perfect reproducibility is not a goal.

## Toast messages

- Copy success: "Copied figure to clipboard"
- Copy failure (no clipboard support): "Clipboard not available — try Download instead"
- Copy failure (denied / error): "Couldn't copy figure: <reason>"
- Download: silent (the file lands; the OS confirms)

## Risks

- **Mark-builder duplication drift.** The export marks for trace + peaks + predicted-q ticks live in a separate file from the on-screen overlay logic. New peak-source types or tick-geometry changes need updates in both places. Mitigation: name the export-mark file alongside the on-screen file in CLAUDE.md so the dependency is discoverable; structural snapshot test asserts the basic mark inventory.
- **Font fidelity.** Using `ui-sans-serif` instead of Plus Jakarta Sans means exported figures don't typographically match the app. Accepted: the "progress report" framing lowers the bar, and embedding the app's webfont would require base64-encoding into every export — not worth the bytes.
- **Browser compat on Firefox.** Older Firefox stable does not support `image/png` `ClipboardItem`. Mitigation: pre-flight check + graceful fallback to Download.
