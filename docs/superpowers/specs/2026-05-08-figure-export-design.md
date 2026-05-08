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

Peak-source key reads "auto, manual, and excluded auto peaks" — `excluded` is a boolean on `Peak` (an auto peak the user explicitly hid), not a third source value.

R² for active phases — deferred.

### MultiTracePlot export — "trace overlay across these members"

Self-explanatory contents:

- **Title:** comparison title — the primary line is `Comparison.title` (the on-disk field name; the adapter arg is `comparisonTitle` to match). When the parent experiment is known (per-experiment scope, `/experiments/:eid/compare/:id`), a secondary line shows the experiment name. In **global scope** (`/compare/all`), `eid` is `undefined` and `Comparison` has no `experiment_id` field — the secondary line is omitted entirely (not "all experiments"; the absence is intentional).
- **Stacked traces** in current band layout, current `xDomain`. Members with `exposure_id === null` (placeholder rows that haven't been bound to an exposure, or whose exposure was deleted) are **filtered out** before plotting — they can't index into the `traces` map, and the gutter renders them as `"(deleted exposure)"` in `MemberMetaRow.defaultLabel`. The disabled-state gate is augmented to also disable export when *every* member is null-exposure (band collapse → empty figure).
- **Per-member labels** rendered inline at each band's vertical position, inside the SVG. Today these labels live in `MemberMetaGutter`, a sibling component to the plot; the export folds them into the figure so it reads alone. The label string is `defaultLabel(member)` from `MemberMetaRow.tsx:61-65` — either `member.label_override` if set, or `"Exposure #${member.exposure_id}"` (post-filter, exposure_id is non-null). Sample names, R-factors, and other metadata are out of scope for v1, matching the on-screen gutter's primary line.
- **Honours current toggles** (`showPeakTicks`, `showPeakLabels`) — these already exist on the page and are the user's "what to include" controls; respecting them means no separate export config.
- **Per-member colours** are resolved via `colorFor(member, ctx, mode)` from `src/lib/comparison/coloring.ts` — the same function the on-screen plot uses. The export adapter pre-computes a `Map<number /* member.id */, string /* css colour */>` once and threads it into the export-marks builder, so the mark builder doesn't repeat colour resolution per band.
- **Legend** depends on `groupingMode`:
  - `"bySample"` — legend is a sample → colour key (each unique `sampleIdFor(m)` gets one row).
  - `"byPhase"` — legend is a phase → colour key (each unique `member.snapshot.confirmed_index.phase` gets one row), matching `colorFor`'s palette for that mode.
  - `"distinct"` — no legend; per-member labels inline at each band already encode the colour.

## Architecture

```
src/lib/figure-export/
  types.ts                    — ExportSpec, LegendSpec
  renderer.ts                 — buildExportSvg, buildExportPng
  clipboard.ts                — copyPngToClipboard
  download.ts                 — downloadBlob
  filename.ts                 — slugifyForFilename, buildFilename
  presets.ts                  — EXPORT_DIMS, fonts, strokes, LIGHT_PALETTE
                                 (TRACE_DIMS = 800×600, COMPARE_DIMS = 1000×400;
                                  numeric values pinned here to keep adapters terse
                                  and to make figure-size changes a one-line edit)
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
import * as Plot from "@observablehq/plot";

interface ExportSpec {
  title: { primary: string; secondary?: string };
  width: number;
  height: number;
  plot: Plot.PlotOptions; // marks + x/y configs for Plot.plot() (only a subset
                          // of PlotOptions fields are populated by the adapter:
                          // typically marks, x, y, width, height, margin*, style)
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
  qUnits?: string;             // matches TraceViewerProps; falls back to "A-1"
}): ExportSpec

import type { GroupingMode } from "../../comparison/coloring";

function buildMultiTraceExportSpec(args: {
  members: ComparisonMember[];           // already in display_order; null-exposure members
                                         // are filtered by the adapter before plotting
  traces: Map<number, Trace>;            // keyed by exposure_id, shape { q, I, sigma }
  comparisonTitle: string;               // = Comparison.title from api.ts
  experimentName?: string;               // optional: undefined in /compare/all global scope
                                         // (Comparison has no experiment_id). When defined,
                                         // ComparePage sources it via useExperiment(eid).
  xDomain: [number, number] | null;
  showPeakTicks: boolean;
  showPeakLabels: boolean;
  groupingMode: GroupingMode;            // "bySample" | "byPhase" | "distinct" — drives
                                         // colorFor() and legend shape
  sampleIdFor: (m: ComparisonMember) => number | null; // already used by MultiTracePlot
  // highlightedMemberId is intentionally NOT in the input — like the trace
  // adapter's hoveredIndex, transient highlight state is excluded from the export.
}): ExportSpec
```

### Renderer

```ts
buildExportSvg(spec: ExportSpec): SVGSVGElement
  // 1. Plot.plot(spec.plot) → inner SVGSVGElement.
  //    Adapters MUST NOT set `title`, `caption`, or `figure: true` on
  //    `spec.plot` — those make Plot wrap the SVG in a <figure> and return
  //    HTMLElement instead of SVGSVGElement, which breaks the rest of the
  //    pipeline. Title and legend rendering belong to the renderer, not
  //    Plot's wrapper.
  // 2. Wrap in outer SVGSVGElement with explicit xmlns: white-bg <rect>,
  //    title <text> nodes, legend rows beneath
  // 3. Marks were built with literal palette colours at adapter time, so no
  //    var(--color-*) resolution is needed at render time

async buildExportPng(spec: ExportSpec, scale = 2): Promise<Blob>
  // 1. const svg = buildExportSvg(spec)
  // 2. const url = URL.createObjectURL(new Blob([new XMLSerializer().serializeToString(svg)],
  //                                            { type: "image/svg+xml" }))
  // 3. try {
  //      const img = new Image();
  //      img.src = url;
  //      await img.decode();
  //      const off = new OffscreenCanvas(spec.width * scale, spec.height * scale);
  //      off.getContext("2d")!.drawImage(img, 0, 0, off.width, off.height);
  //      return await off.convertToBlob({ type: "image/png" });
  //    } finally { URL.revokeObjectURL(url); }
  //
  // `width`/`height` in ExportSpec are logical (CSS) pixels. scale=2 produces a
  // 2× retina-quality raster suitable for print and HiDPI displays.
```

`OffscreenCanvas.convertToBlob` is genuinely async (returns `Promise<Blob>`); no callback wrapping needed. The `finally` is load-bearing — without `URL.revokeObjectURL` every export leaks an object URL. Browsers without `OffscreenCanvas` support are out of scope (Safari 16.4+, Chromium evergreen, Firefox 105+ all have it); we don't propose an `HTMLCanvasElement.toBlob` fallback because the `OffscreenCanvas` floor matches HimalayaUI's existing browser-support policy.

### Colour resolution

The export-only mark factories accept literal colours from `LIGHT_PALETTE` (in `presets.ts`) directly — they don't read CSS variables. The rendered SVG never contains `var(--color-*)`, so the renderer doesn't need a string-replace post-pass.

### Button component

```tsx
<FigureExportControls
  spec={() => buildTraceExportSpec({ /* ... */ })}
  filenameStem="himalaya-trace-jc23-sample4-exp7"   // already-resolved by parent
  ariaContext="trace plot"                           // used in aria-labels
/>
```

**`filenameStem` contract**: the parent passes a fully-resolved, already-slugified stem. The component does NOT do template substitution. Internally it appends `-{YYYY-MM-DD}.{ext}` (date is **local time**, not UTC — local matches user expectation when exporting near midnight) and the format extension. The `slugifyForFilename` helper in `filename.ts` is for the parent to call when building the stem from name fields.

**`spec` is a thunk** evaluated at click time so it captures fresh state — `xDomain`, current peaks, active group, etc. — rather than stale state captured at component mount.

**Internal state**: a `pending` boolean disables both buttons during an in-flight render to prevent double-click re-entry. The `pending` flip-back lives in a `try { … } finally { setPending(false) }` so a thrown render (malformed SVG → `Image.decode()` rejects, `convertToBlob` throws, `clipboard.write` errors) does not brick the button until refresh.

**Render-cost assumption.** "Sub-second render" is an assumption to verify during implementation, not a measured number. If profiling shows >500ms (likely candidates: a TraceViewer export with many predicted-q ticks across multiple high-rank phases, plus a 2× DPI raster on a slower laptop), replace the `pending` flag with an `AbortController` so the user can cancel an in-flight render. Calling this out explicitly here so it doesn't get missed.

**Buttons:**
- **Copy** — single icon button. `aria-label="Copy {ariaContext} to clipboard"`.
- **Download (split)** — primary icon button + sibling chevron button.
  - Primary: `aria-label="Download {ariaContext} as PNG"`. Click → download PNG directly.
  - Chevron: `aria-label="Other download formats"`, `aria-haspopup="menu"`, `aria-expanded`. Click → toggle a small popover `<ul role="menu">` with "Download as PNG" / "Download as SVG" rows.

**Popover behaviour** — match `ForksPopover` pattern (lightweight, click-to-toggle, NOT focus-trapped because it isn't a modal). Add Esc-to-close and click-outside dismiss for accessibility — these were noted as future enhancements in `ForksPopover` and are appropriate here. Do not call `useFocusTrap`.

**Note on render cost:** The parent creates a new `spec` function reference on every render. This does not affect correctness — the thunk is only evaluated on click — but it defeats `React.memo` on `FigureExportControls` if applied. In practice the component is small (two buttons) so the extra re-render is negligible. If profiling later shows a concern, wrap the thunk in `useCallback` with stable deps.

### Placement

- **PlotCard title strip** (Index page): append after the q-range reset button, separated by a thin divider. The strip already holds fit-features + x-scale + q-range; placing export controls *next to* the controls that change what gets exported keeps the adjacency right.
- **ComparePage review header** (`compare-review-header`): append at the right end of the toggles row, after `ForksPopover`.

### Data flow at click time

```
user clicks Copy
  → FigureExportControls sets `pending = true` (buttons disabled)
  → try {
      → invokes spec()
      → adapter pulls live state via the closure
      → renderer.buildExportPng(spec) → Blob
      → clipboard.copyPngToClipboard(blob)
      → showToast("Copied figure to clipboard", "success")
    } catch (err) {
      → showToast("Couldn't copy figure — try Download instead.", "error")
      → console.warn(err)
    } finally {
      → setPending(false)
    }
```

No state lifted to Zustand or TanStack Query — the export is a pure derivation of state already present in the parent component when the user clicks.

## Browser compatibility

### Clipboard

`navigator.clipboard.write([new ClipboardItem({ "image/png": blob })])`. Three failure modes handled:

- **Insecure origin** — `navigator.clipboard` undefined. Pre-flight check; if missing, Copy is disabled with a tooltip "Clipboard requires HTTPS."
- **`ClipboardItem` unavailable** — same path; disable Copy, leave Download.
- **Runtime denial** (user-interaction policy, browser block) — caught and surfaced as `showToast("Couldn't copy. Try Download instead.", "error")`.

### SVG → PNG canvas pipeline

- Use `font-family: ui-sans-serif, system-ui, sans-serif` for export. Avoids `@font-face url(...)` references → no canvas-taint risk. Won't match the app's Plus Jakarta Sans exactly; accepted as the "progress report" tradeoff.
- Outer wrapper SVG must declare `xmlns="http://www.w3.org/2000/svg"`. Plot's inner SVG already does; the outer wrapper sets it explicitly.
- Standard recipe: serialize SVG → blob URL → `Image` → `OffscreenCanvas.drawImage` → `convertToBlob({ type: "image/png" })`.
- Browsers without `OffscreenCanvas` (Safari < 16.4) are out of scope; we deliberately do **not** propose an `HTMLCanvasElement.toBlob()` fallback. `toBlob` is callback-based, not synchronous, so a fallback would require wrapping in `new Promise(resolve => canvas.toBlob(resolve))` — and the support floor (Safari 16.4+, Chromium evergreen, Firefox 105+) matches HimalayaUI's existing browser-support policy, so the extra code path isn't worth it.

## Disabled states

- **PlotCard**: Copy + Download disabled if `!traceQ.data || !peaksQ.data`. This is intentionally **stricter** than the existing `canFit` gate (which only checks `traceQ.data`) because the export emits peak markers and predicted-q ticks as part of the figure design — without `peaksQ.data`, the figure misses load-bearing content. The stricter gate is correct, not an oversight.
- **ComparePage**: disabled if `members.length === 0 || traces.size === 0 || members.every(m => m.exposure_id == null)`. The third clause covers the all-placeholder case where every row is unbound — without it the export would render a band-collapse with no traces. Here `traces` is the `Map<number, Trace>` returned by `useMemberTraces(exposureIds)`, where `Trace = { q, I, sigma }`.

## Filenames

- TraceViewer parent builds: `himalaya-trace-{slug(experiment)}-{slug(sample)}-{slug(exposure)}` — passed to `FigureExportControls` as `filenameStem`. The component appends `-{YYYY-MM-DD}.{ext}`.
- Compare parent builds: `himalaya-comparison-{slug(experimentName ?? "all")}-{slug(comparison.title)}` — same. (`name` here resolves from `Comparison.title` — the on-disk field. The fallback `"all"` covers `/compare/all` global scope where there's no parent experiment.)

The date is **local time** (in the user's timezone, not UTC — local matches user expectation when exporting near midnight). The locale is **pinned to `en-CA`** which always produces `YYYY-MM-DD`:

```ts
new Intl.DateTimeFormat("en-CA", { year: "numeric", month: "2-digit", day: "2-digit" }).format(d)
```

Default-locale formatting is unsafe — `en-US` produces `05/08/2026`, and `/` is an invalid filename character on every major OS. `en-CA` is locale-pinned for filename safety only; it does not localize anything user-facing. (Equivalent alternative: `formatToParts()` and assemble manually. Either is fine; pick one and stick to it.)

`slugifyForFilename(s)`: lowercase, replace non-alphanumeric runs with `-`, collapse repeated dashes, trim leading/trailing dashes. Empty input or input that slugifies to an empty string returns the sentinel `"figure"`.

**Applied per segment**, not to the concatenated filename — each template variable (`slug(experiment)`, `slug(sample)`, etc.) is slugified individually before interpolation. The static separator dashes between segments remain unambiguous delimiters, so a slugified segment containing dashes (e.g., an experiment named "SAXS-run-042") cannot be confused with the boundary between segments.

## Testing

Three layers, matching the existing project conventions.

### Vitest (unit)

All test files live in `frontend/test/figure-export/` (matching the project's `frontend/test/` convention seen in `test/coloring.test.ts`, `test/TraceViewer.test.tsx`, etc.).

- `traceAdapter.test.ts` — given fixture state, returns expected `ExportSpec` shape (title parts, dims, presence of tick marks per index, legend rows).
- `multiTraceAdapter.test.ts` — same for Compare, with one case per `GroupingMode` value (`"bySample"`, `"byPhase"`, `"distinct"`) to verify legend variation. Also a case where every member has `exposure_id === null` (asserts the disabled-state condition fires) and a case where one member is null-exposure (asserts the filter drops it from the figure).
- `renderer.test.ts` — `buildExportSvg(spec)` produces an SVG with no `var(--…)` literals (string regex), title `<text>` node present at expected position, legend group at expected y, white background `<rect>`. Also: round-trip `XMLSerializer.serializeToString(svg)` and re-parse to confirm the output is well-formed XML. Also: assert `Plot.plot(spec.plot)` returns `SVGSVGElement` (not `HTMLElement`) — guards against an adapter accidentally setting `title` / `caption` / `figure: true`.
- `filename.test.ts` — slugify covers spaces, slashes, parens, unicode, empty input. Empty/all-special input returns `"figure"`. Date suffix uses `en-CA` locale and produces `YYYY-MM-DD` regardless of system locale (test by stubbing `Date`; verify with `vi.setSystemTime` and check the formatted output literally equals e.g. `"2026-05-08"`).
- `clipboard.test.ts` — mock `navigator.clipboard.write` via `vi.stubGlobal("navigator", { ...navigator, clipboard: { write: vi.fn() } })`; assert the expected `ClipboardItem` was passed. (Cleaner than `Object.defineProperty` on the read-only `clipboard` getter; restores correctly between tests.)

### JSDOM constraint

`OffscreenCanvas.convertToBlob()` and `Image.decode()` aren't real in JSDOM. The PNG path is smoke-testable only (assert the function exists and returns a Promise). The SVG path is fully testable in JSDOM.

### Playwright E2E (mocked)

Two specs in `e2e/`:
- "Copy figure → clipboard contains PNG" — Chromium-only (WebKit/Firefox restrict `clipboard.read({ image/png })` even with permissions). Set `permissions: ['clipboard-read', 'clipboard-write']` in `playwright.config.ts` `use:` block (project-level — `context.grantPermissions(...)` inside the test fires too late for some Chromium versions). Click Copy, read clipboard via `page.evaluate(() => navigator.clipboard.read())`, assert at least one item with `image/png` MIME. The test skips on non-Chromium projects via `test.skip(browserName !== "chromium", ...)`.
- "Download → PNG / SVG" — cross-browser. `page.waitForEvent('download')`, assert the suggested filename matches the pattern and the file is non-empty. This is the contingency test: even on browsers where the clipboard test is skipped, download is fully covered.

### Visual regression

Out of scope. The "good enough for progress reports" framing means structural snapshot tests cover drift; pixel-perfect reproducibility is not a goal.

## Toast messages

The toast API is `showToast(msg: string, kind: ToastKind)` from `src/lib/toast.ts` (where `ToastKind = "info" | "success" | "warning" | "error"`).

- Copy success: `showToast("Copied figure to clipboard", "success")`
- Copy failure (denied / runtime error): `showToast("Couldn't copy figure — try Download instead.", "error")` — the raw error object is logged to `console.warn` for debugging but not surfaced to the user, who cannot act on browser-internal error messages like "ClipboardItem's type image/png is not supported".

There is **no** "no clipboard support" toast: the pre-flight check disables the Copy button when `navigator.clipboard` or `ClipboardItem` is missing, so the click handler can't fire. The disabled state surfaces a tooltip (`title="Clipboard requires HTTPS"` or similar) — that's the user-facing signal.
- Download: silent (the file lands; the OS confirms)

## Risks

- **Mark-builder duplication drift.** The export marks for trace + peaks + predicted-q ticks live in a separate file from the on-screen overlay logic. New peak-source types or tick-geometry changes need updates in both places. Mitigation: name the export-mark file alongside the on-screen file in CLAUDE.md so the dependency is discoverable; structural snapshot test asserts the basic mark inventory.
- **Font fidelity.** Using `ui-sans-serif` instead of Plus Jakarta Sans means exported figures don't typographically match the app. Accepted: the "progress report" framing lowers the bar, and embedding the app's webfont would require base64-encoding into every export — not worth the bytes.
- **Browser compat on Firefox.** Older Firefox stable does not support `image/png` `ClipboardItem`. Mitigation: pre-flight check + graceful fallback to Download.
