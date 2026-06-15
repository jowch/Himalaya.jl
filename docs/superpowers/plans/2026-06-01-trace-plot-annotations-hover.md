# Trace-Plot Toggleable Annotations + Hover Interaction — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]`.
> Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`

**Goal:** Make the greenfield trace-plot engine's annotation layers individually toggleable, colour peaks per-assigned-index, and add the mockup's hover interaction (hot peak + q-line + a q-value readout), with an index/phase highlight that fades others to neutral — all engine-level and proven in Storybook (no Focus/Series wiring; that's Plan #2).

**Architecture:** Pure additions to `src/print/plot/`. A `show` prop gates layers; `PlotPeak` carries per-peak `color`+`label`; a new labels mark renders mono Miller/q text with a horizontal dodge; the engine owns hover via a `PlotFrame` pointer-move callback feeding `hitTestPeaks`, driving a hot peak + an axis-anchored q chip; a `highlightPeakIds` set fades non-members to neutral gray. Hover has a keyboard/focus equivalent; motion respects `prefers-reduced-motion`.

**Tech Stack:** React 18 + TS strict (`exactOptionalPropertyTypes`), d3-scale/shape, Vitest+RTL, Storybook CSF3. `src/print/plot/**` design-guard-exempt.

---

## Context

The Plan-#1 engine + the just-shipped "The Print" visual refit render the line, σ band, decade grid, and peak glyphs. Three gaps remain, all confirmed against the mockups/docs:
1. **No per-layer visibility control.** Consumers can't show/hide peaks, labels, or the band independently.
2. **Peak colour is per-trace, not per-peak.** Real assigned-index colouring is per-peak (`claimColorByPeakId`, `TraceViewer.tsx:301`): a peak takes the colour of whichever index claims it; unclaimed auto → `--color-ink-faint`.
3. **No hover interaction and no peak labels.** The mockup defines `.hot`/`.preview`/`.dim` (90ms), a `.pk-qline`, and mono `.pk-label`s (Miller/`?`), but the engine has none of it; `PlotFrame` only does wheel/click/dblclick.

### Locked decisions
- **q-readout:** axis-anchored q chip on hover (mono, contrast-safe ink); NOT the label-swap.
- **dim:** peak-hover does NOT dim others (matches focus mockup). Dimming belongs to `highlightPeakIds` (index/phase highlight) and fades non-members to **neutral gray** per `himalayaui-design.md §1.3`, not paler-hue.
- **a11y (impeccable C + audit P1):** the hover highlight has a **focus/keyboard equivalent** (focusing a peak = hovering it); all new transitions are wrapped in `@media (prefers-reduced-motion: reduce)` no-op.
- **labels source:** the engine does NOT compute Miller (needs the index). `PlotPeak.label` is a pre-resolved string passed by the consumer; Storybook synthesizes labels from fixtures.

### Verbatim values from the mockup (`docs/redesign-mockups/2026-05-29-focus-plot.html`)
| Element | Value |
|---|---|
| `.pk-qline` | width 1; opacity 0 → **0.75** hot → **0.85** preview; `transition: opacity 90ms ease-out` |
| `.pk-label` | mono **11px / weight 700**, fill = peak colour, `text-anchor: middle`, y = glyph apex − **11px** |
| `.pk.dim` glyph+label | opacity **0.32** (we override to neutral-gray fade for index-highlight per §1.3) |
| hot glyph | grow ×1.5 + terracotta **ring** (already in `peakGlyph`, via `hot`) — never a recolour |
| label dodge | pixel-space horizontal dodge, `LABEL_OFFSET_PX = 12`, leader line when shifted (`lib/plot/labelDodge.ts`) |
| Miller map | `{2:110, 3:111, 4:200, 6:211, 8:220, 9:221}` else `√N` (consumer-side, not engine) |

---

## Task 1 — `show` layer-visibility prop

**Files:** Modify `src/print/plot/TracePlot.tsx`, `marks/TraceLine.tsx`; Test `test/print-plot/TracePlot.test.tsx`

Add an optional `show` prop controlling which annotation layers render. Grid + axes stay governed by the existing `axes` prop.

- [ ] **Test first** (TracePlot.test.tsx): with `show={{ band: false }}` no `[data-role="trace-line"] path[opacity]` (band) renders but the line does; with `show={{ peaks: false }}` no `[data-role="plot-peaks"]` content; default (no `show`) renders peaks + band as today.
- [ ] **Implement:**
  - Add to `TracePlotProps`: `show?: { peaks?: boolean; labels?: boolean; band?: boolean }`.
  - Resolve with defaults at top of `TracePlot`: `const layers = { peaks: true, labels: false, band: true, ...show };` (labels default **off** — opt-in).
  - Pass `band={layers.band}` to `<TraceLine>` (TraceLine already has a `band` prop — wire it through instead of its internal default).
  - Gate the `<PlotPeaks>` map on `layers.peaks`.
  - (Labels layer added in Task 3, gated on `layers.labels`.)
- [ ] Run `npm test -- TracePlot` → PASS. Commit `feat(plot): show={} layer-visibility prop`.

---

## Task 2 — Per-peak colour + label on `PlotPeak`

**Files:** Modify `src/print/plot/marks/PlotPeaks.tsx`, `TracePlot.tsx`; Test `test/print-plot/PlotPeaks.test.tsx`

- [ ] **Test first:** a peak with `color: "var(--color-success)"` renders its glyph stroke/fill in that colour (overriding the layer `color`); a peak without `color` uses the passed `color` prop (back-compat).
- [ ] **Implement:**
  - Extend `PlotPeak`: add `color?: string` (resolved per-peak colour) and `label?: string` (pre-resolved label text).
  - In `PlotPeaks`, the per-glyph colour becomes `const c = p.color ?? color;` and feeds `peakGlyph({ color: c, ... })` and the hot q-line stroke. The `color` prop stays the layer fallback.
  - `TracePlot` still passes the trace's phase colour as the layer fallback; no change needed there beyond keeping the existing `color={t.phase ? phaseColor(t.phase) : UNINDEXED_COLOR}`.
- [ ] Run `npm test -- PlotPeaks` → PASS. Commit `feat(plot): per-peak color + label fields on PlotPeak`.

---

## Task 3 — Peak labels layer with horizontal dodge

**Files:** Create `src/print/plot/marks/PlotLabels.tsx`; Modify `TracePlot.tsx`; Create `test/print-plot/PlotLabels.test.tsx`

A new mark renders each peak's `label` (skip empty/excluded) as mono text above the glyph, dodged horizontally so dense labels don't overlap, with a leader line when a label is shifted off its peak.

- [ ] **Test first** (PlotLabels.test.tsx): given 2 peaks with labels and a projection, renders 2 `[data-role="peak-label"]` `<text>` nodes with the label text; a peak with `label: ""` or `excluded` renders none; tightly-spaced peaks produce distinct (non-equal) label x-positions (dodge applied).
- [ ] **Implement** `PlotLabels({ peaks, projection, color, baselineI, offsetPx = 12, labelWidthPx = 30 })`:
  - Compute each labelled peak's `(px, py)` like `PlotPeaks` (skip `id < 0`, `excluded`, empty/undefined `label`).
  - Dodge: port the simplified two-pass pixel dodge from `src/lib/plot/labelDodge.ts` (sort by px; left-to-right enforce ≥ `labelWidthPx` gap; recenter on the mean so the cluster stays symmetric; snap unmoved labels to natural px). Keep it pure and unit-testable (a local `dodgeX(positions, labelWidthPx)` helper).
  - Render per label: a `<text data-role="peak-label">` at `(labelX, py - offsetPx)`, `textAnchor="middle"`, style `{ fontFamily: "var(--font-mono)", fontSize: 11, fontWeight: 700, fill: p.color ?? color }`; when `labelX !== px`, a thin leader `<line>` from `(px, py - 6)` to `(labelX, py - offsetPx + 2)` in `--color-hair`.
  - Wrap in `<g data-role="plot-labels">`.
  - In `TracePlot`, render `<PlotLabels>` per trace, gated on `layers.labels`, after `<PlotPeaks>`.
- [ ] Run `npm test -- PlotLabels` and `npm test -- TracePlot` → PASS. Commit `feat(plot): toggleable peak-label layer with horizontal dodge`.

---

## Task 4 — Engine-owned hover: hot peak + q-line + axis q-chip (+ keyboard)

**Files:** Modify `src/print/plot/PlotFrame.tsx`, `TracePlot.tsx`; add a small `QReadout` render; `src/styles.css` (reduced-motion guard if a class is used); Test `test/print-plot/TracePlot.test.tsx`

- [ ] **PlotFrame:** add optional `onPointerMovePx?: (px, py) => void` and `onPointerLeave?: () => void`; bind `onPointerMove`/`onPointerLeave` on the container (container-relative px, same rect math as the wheel handler). Test: moving the pointer fires `onPointerMovePx` with container-relative coords.
- [ ] **TracePlot hover state:** add `const [hoverId, setHoverId] = useState<number | null>(null)`. On pointer move (only when `interaction !== false`): subtract margins, hit-test with `hitTestPeaks(peaks, plotPx, q => projection.x.to(q), tol)`; `setHoverId(hit?.id ?? null)`. On leave: `setHoverId(null)`. The hovered peak is rendered **hot** (pass `hot: true` into its `peakGlyph`/q-line) — merge with any externally-supplied `hot`.
- [ ] **q-chip:** when `hoverId` is set, render an axis-anchored readout: a `<g data-role="q-readout">` at `x = projection.x.to(hoveredPeak.q)`, `y = plotHeight`, containing a small rounded `<rect>` (fill `--color-plate`, stroke `--color-hair-strong`) + mono `<text>` (`var(--font-mono)`, 10.5px, fill `--color-ink` — contrast-safe, NOT ink-faint) showing the q value formatted via `formatAxis(q)` (or `q.toFixed(3)`). Position centered under the q-line foot.
- [ ] **Keyboard/focus equivalent (impeccable):** make each peak hit-target focusable. Simplest: in `PlotPeaks`, wrap each peak `<g>` with `tabIndex={0}` + `role="button"` + `aria-label={`peak at q ${formatAxis(p.q)}`}` and `onFocus`/`onBlur` that call back up to set the same hover id (add an `onPeakFocus?: (id|null)=>void` prop to `PlotPeaks`, wired from `TracePlot` to the same `setHoverId`). Focusing a peak shows the same hot + q-chip. Test: firing `focus` on a peak `<g>` shows `[data-role="q-readout"]`.
- [ ] **Reduced-motion:** the hover/q-line opacity transitions (90ms) must be disabled under `prefers-reduced-motion`. If transitions are expressed via a CSS class in `styles.css`, add the `@media (prefers-reduced-motion: reduce) { transition: none }` guard; if inline, gate the transition on a `usePrefersReducedMotion()` check (or simply omit transitions and rely on React state swaps — acceptable). Document the choice in a comment.
- [ ] Run `npm test -- TracePlot` and `npm test -- PlotFrame` → PASS. Commit `feat(plot): engine-owned hover — hot peak, q-line, axis q-readout, focus equiv`.

---

## Task 5 — `highlightPeakIds`: index/phase highlight fades others to neutral

**Files:** Modify `src/print/plot/TracePlot.tsx`, `marks/PlotPeaks.tsx`, `marks/PlotLabels.tsx`; Test `test/print-plot/PlotPeaks.test.tsx`

- [ ] **Test first:** with `highlightPeakIds={new Set([1])}` on a 3-peak trace, peak 1 renders in its phase colour while peaks 2 & 3 render with the neutral-gray fade (assert via a `data-dimmed="true"` attribute on the non-members' glyph `<g>`); with no `highlightPeakIds` (or empty), nothing is dimmed.
- [ ] **Implement:**
  - Add `highlightPeakIds?: ReadonlySet<number>` to `TracePlotProps`; thread to `PlotPeaks`/`PlotLabels`.
  - In `PlotPeaks`/`PlotLabels`, when the set is non-empty and a peak's id is NOT in it, render the glyph/label in neutral gray (`--color-ink-faint`) — fade to neutral per §1.3, not opacity-only — and set `data-dimmed="true"`. Members render normally (their per-peak colour). Implement by overriding the resolved colour `c` to `var(--color-ink-faint)` for non-members (the silhouette/second-channel still distinguishes provenance).
  - Hover (`hoverId`) and `highlightPeakIds` compose: a hovered peak is always hot regardless of dim.
- [ ] Run `npm test -- PlotPeaks` → PASS. Commit `feat(plot): highlightPeakIds index/phase highlight (fade others to neutral)`.

---

## Task 6 — Storybook stories + build gate

**Files:** Modify `src/print/plot/TracePlot.stories.tsx`

- [ ] Add stories demonstrating the new capabilities, using the real fixtures + synthesized labels/colours:
  - **`Annotated`** — labels on (synthesize Miller-ish labels + per-peak phase colours on `heroModel`'s peaks), band on, with Storybook `argTypes` toggles for `show.peaks` / `show.labels` / `show.band` so the layers flip live.
  - **`Hover`** — note in the story description that pointer/focus over a peak shows the hot state + axis q-chip (interactive; nothing to assert in CSF, it's a manual-fidelity story).
  - **`PhaseHighlight`** — pass `highlightPeakIds` for one phase's peaks to show members lit + others faded to neutral.
  - Keep `Hero`/`Mini`/`TransitionOverlay`/`HeroPlate`.
- [ ] **Verify:** `npm test -- print-plot` all green; `npm run build` exit 0; `npm run build-storybook` exit 0.
- [ ] Commit `feat(plot): Storybook stories for annotations, hover, phase-highlight`.

---

## Verification (end-to-end)

1. `cd packages/HimalayaUI/frontend`
2. `npm test -- print-plot` — all unit tests green (projection, Axis, TraceLine, PlotPeaks, PlotLabels, PlotFrame, TracePlot).
3. `npm run build` — lint:design + tsc + vite all clean.
4. `npm run build-storybook` — succeeds.
5. **Manual fidelity:** `npm run storybook` → `plot/TracePlot`:
   - `Annotated` — flip the `show.*` arg toggles; labels/peaks/band appear/disappear independently.
   - `Hover` — pointer over a peak: it grows + terracotta ring, q-line shows, q-chip appears at the axis foot; **Tab** to a peak gives the same highlight (keyboard equiv); others do NOT dim.
   - `PhaseHighlight` — highlighted phase's peaks stay coloured; the rest fade to neutral gray.

## Out of scope (deferred → Plan #2 / consumer wiring)
- Computing Miller labels / claim-colours from the live assignment model (`IndexEntry` → `claimColorByPeakId`), and wiring real Focus/Series controls.
- `.preview` candidate ghost-comb, detector-ring q-link, series migration-track (separate surfaces).
- Tap-target enlargement for coarse pointers beyond the existing `PEAK_HIT_PX` (the audit's `@media (pointer:coarse)` pass).
