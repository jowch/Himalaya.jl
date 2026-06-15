# Greenfield Comb + Residual Charts — Design

> Part of the greenfield "The Print" rebuild (branch `worktree-greenfield-ui-rebuild`, **not merged** — the whole rebuild lands as one piece). Renderer track, following the detector renderer (Batch 2a). Date: 2026-06-02.

## Purpose

Two standalone, **pure SVG renderers** for the reflection-comb / indexing diagnostic on the Focus page:

- **`CombChart`** — the reflection comb. For each assigned phase, draws the predicted Bragg reflection q-positions ("teeth") on a shared log-q rail, so the user can see how the predicted series lands against the observed peaks and whether reflections are missing.
- **`ResidualChart`** — the goodness-of-fit diagnostic. **Structurally tracks `CombChart`** (same gutter, same rows in the same vertical positions, same q-axis), but inside each phase row it draws a local zero baseline and that phase's fractional residuals `(q_obs − q_pred)/q_pred`, with the phase's R² in the gutter.

Both are pure: props in, hover callbacks out. No Zustand, no TanStack Query, no URL building, no API calls. The future `CombsPanel` composite owns data fetching, the comb/residual SegmentedControl toggle, Zustand `hoveredQ` wiring, and legend placement.

## Why two components, not one with a `view` prop

The comb (multi-row teeth, hkl labels, leftover rings) and the residual (zero-line scatter, tolerance band, R²) share almost nothing *internally* — only the log-q x-projection and the rowed scroll/gutter scaffold. A single `view`-prop component would be one large branchy file doing two jobs (violates "one responsibility per file"). Two siblings are each independently testable and iterable in Storybook, and the genuinely-shared geometry is factored into small helpers. They are *visually* parallel by design (see "tracked layout" below), which is a property of sharing the scaffold — not a reason to fuse them.

## Location & design-guard

New design-guard-exempt directory **`src/print/comb/`** (appearance literals allowed, like the renderers before it). Add `relPath.startsWith("print/comb/")` as the 5th prefix in `isExcluded()` in `scripts/check-design.mjs` (currently: `components/ui/`, `print/ui/`, `print/plot/`, `print/detector/` — at `scripts/check-design.mjs:137-143`).

Files:

```
src/print/comb/
  combModel.ts        # view-model types (CombSeries, CombTooth) + pure derivations
  CombScaffold.tsx    # shared rowed layout: pinned left gutter + horizontally-scrollable log-q pane + bottom q-axis
  CombChart.tsx       # teeth renderer (consumes CombScaffold)
  ResidualChart.tsx   # residual renderer (consumes CombScaffold)
  CombChart.stories.tsx
  ResidualChart.stories.tsx
  index.ts            # barrel
test/print-comb/
  combModel.test.ts
  CombChart.test.tsx
  ResidualChart.test.tsx
```

`CombScaffold` may instead start inline and be extracted once the second consumer (ResidualChart) confirms the shared shape — but since we *know* both consumers up front, the plan builds it as a shared unit from the start.

## Data model (pure view-model, decoupled from the API)

The charts consume a render-ready shape. The composite maps the backend `IndexEntry` (`src/api.ts:276-293`: `phase`, `predicted_q: number[]`, `peaks: IndexPeakRef[]`, `r_squared`, `lattice_d`) into it later — that mapping is **out of scope** here.

```ts
export interface CombTooth {
  q: number;            // predicted reflection q (Å⁻¹)
  label: string;        // hkl or √N label, e.g. "√6" or "211"
  observed: boolean;    // true → a claimed observed peak exists (solid tooth); false → predicted-absent (caret)
  residual?: number;    // fractional residual (q_obs − q_pred)/q_pred for the claimed peak; used by ResidualChart
}

export interface CombSeries {
  phase: string;        // "Pn3m"
  color: string;        // phaseColor(phase) — passed in by the composite (token ref, allowed)
  latticeLabel?: string; // comb gutter line 2, e.g. "a = 197 Å"  (from lattice_d)
  rSquared?: number;    // residual gutter line 2, e.g. 0.998      (from r_squared)
  teeth: CombTooth[];
}
```

Row source (from the user's "tracked" decision):

- **assigned** rows — `CombSeries[]`, the committed phase assignments (solid baseline).
- **hovered** preview row — one optional `CombSeries`, the candidate currently hovered in the candidate list (dashed baseline, the mockup's "hypothetical" treatment).
- **leftover** — `number[]` of unindexed observed peak q-values. **Comb only**: rendered as a final rings row. Residual omits it (unindexed peaks have no predicted q → no residual).

Component props (both charts share the first block):

```ts
interface CombScaffoldRow { gutterTitle: string; gutterSub?: string; color: string; preview?: boolean; }
// CombChart props
{
  assigned: CombSeries[];
  hovered?: CombSeries;
  leftover: number[];
  hoveredQ?: number;          // incoming q-link; lights the matching tooth/ring across rows
  onHoverQ?: (q?: number) => void;
  maxWidth?: number;          // fixed max; min is enforced by the scaffold's per-decade floor
  className?: string;         // placement only
}
// ResidualChart props: same as above minus `leftover` (no leftover row)
```

## Tracked layout (shared scaffold)

`CombScaffold` renders a flex row:

1. **Pinned left gutter** (fixed width, HTML — *not* scrolled): one block per row, vertically aligned to the row's baseline. Line 1 = phase name (`text-data`); line 2 = `latticeLabel` (comb) or `R² = 0.998` (residual), `text-meta text-ink-faint`. The leftover row's gutter reads "leftover" / unindexed.
2. **Horizontally-scrollable pane**: holds the SVG (rail/rows + teeth/residuals + bottom q-axis).

Because both charts use the same gutter, the same per-row vertical positions, and the same bottom q-axis, toggling comb↔residual in the panel keeps each phase in the same place — a clean per-row cross-fade.

Row geometry: rows stack top-to-bottom in order **assigned → hovered preview → (comb only) leftover**. Each row has a baseline; row height is a constant pixel value (text never scales). The hovered preview row's baseline is dashed.

## Responsive: min-width + horizontal scroll

Text and glyphs are a **constant pixel size** — never scaled down. The SVG pane enforces a **minimum pixels-per-q-decade** so teeth never crowd below legibility. When the container is narrower than the resulting min-width, the **pane scrolls horizontally** while the gutter stays pinned (you always know which row you're reading). `maxWidth` caps the pane; between min and max the pane fills the available width.

Q-axis tick density reuses `axisTicks()` from `src/print/plot/projection.ts:50` — decade-anchored major/mid/minor ticks that already thin labels adaptively, so the axis stays legible and labelled at any width.

## Q-scale: always log (both charts)

Both charts are **fixed to log-q** — no per-panel toggle (clutter), and not coupled to the hero trace's log/linear state. Rationale: the comb is a *ratio* diagram (teeth at √N multiples of q₁ spread evenly in log, crowd at low-q in linear); log is the canonical indexing view and the mockup's choice; it keeps the responsive min-width math stable and the chart self-consistent regardless of the hero. Q-domain derived via `positiveExtent()` (`projection.ts:137`) over all predicted + observed + leftover q-values, with the mockup's ~10% pad.

## Vocabulary, color, hover

Glyph vocabulary follows the **mockup's comb** (`docs/redesign-mockups/2026-05-29-focus-plot.html:880-987`), which is the comb's *intrinsic* vocabulary — distinct from the trace-plot peak glyphs:

- **Observed/claimed reflection** → tooth: phase-colored stem (baseline → apex) + filled circle cap; hkl/√N label above the cap.
- **Predicted-absent** → dashed faint stem (~72% height) + hollow caret (chevron) at the apex.
- **Leftover peak** (comb only) → unindexed ring on the leftover baseline.
- **Residual point** → small filled dot at (q, Δq/q-offset) in the phase color, inside a faint per-row ±~2.2% tolerance band, on a local zero baseline. The per-row residual y-domain is **fixed and symmetric** (e.g. ±3%, band drawn at ±2.2%); a residual beyond the domain clamps to the edge and is marked (e.g. an open/overflow dot) rather than rescaling the row.

Color: from the `color` prop (`phaseColor(phase)`, `src/phases.ts:45`). The renderer is color-agnostic pass-through; the composite supplies it. Unindexed/leftover use a neutral ink token.

**Hover / q-link** (the rule established this session for the detector rings and trace-plot engine — *the accent is never a hover highlight*): an incoming `hoveredQ` that matches a tooth/ring/point makes it "hot" via a **wider stem + larger ringed cap in the row's own phase color** — never a terracotta-accent recolor. This **intentionally diverges from the mockup CSS** (`.tooth.hot { stroke: var(--accent) }`). Hovering a tooth/ring/point emits `onHoverQ(q)`; leaving emits `onHoverQ(undefined)`.

### Legend note (flagged)

The already-built `CombLegend` (`src/print/components/CombLegend.tsx`) legends a **4-item peak-glyph vocab** (triangle/diamond/caret/excluded) — that is really the *trace-plot peak* vocabulary, which differs from the comb's intrinsic tooth/caret/ring. Decision: **CombChart embeds no legend** (the panel places one); the comb's 3-item legend (tooth / absent caret / leftover ring) and the CombLegend naming get reconciled later when the `CombsPanel` composite is built. Not in scope here.

## Reuse vs. own

- **Reuse** from `src/print/plot/projection.ts`: `makeAxis` / `Axis1D` (log-q scale), `axisTicks` (adaptive q-axis ticks), `positiveExtent` (q-domain).
- **Own**: the rowed scaffold, gutter, scroll behavior, teeth/caret/ring glyph geometry, residual zero-line + band + points, hkl labels, hover state. These charts are **not** `TracePlot` and do not extend the trace overlay engine.

## Out of scope (deferred to `CombsPanel` composite / later)

Real API data + `IndexEntry`→`CombSeries` mapping; Zustand `hoveredQ` wiring; the comb/residual SegmentedControl; legend placement + CombLegend reconciliation; pixel-exact alignment to the hero trace; the residual view's interaction with speculative/anchor index building.

## Testing & gates

Vitest structural assertions (via `data-role` / `data-*` / text — never class strings):

- `combModel`: q-domain derivation, tooth `observed` flag from claimed peaks, leftover separation, label/√N formatting.
- `CombChart`: one row per assigned series + a dashed preview row + a leftover ring row; `data-role` `tooth`/`tooth-cap`/`caret`/`leftover-ring`; gutter shows phase name + lattice; hkl label present; min-width floor enforced (pane width ≥ floor at narrow container); hot tooth keeps own color (no `var(--color-accent)`); `onHoverQ` fires q on enter, undefined on leave.
- `ResidualChart`: same rows/positions as CombChart minus leftover; per-row zero baseline; `data-role` `resid-point` at expected y-offset for a known residual; tolerance band present; gutter shows phase name + `R²`; no terracotta on hot.

Stories: a **wide** comb, a **narrow (scrolling)** comb proving the pinned gutter + horizontal scroll, a comb with a hovered preview row, a residual chart (wide + narrow), and an interactive q-link story (comb + residual side by side sharing a `hoveredQ`).

Per component: `npm test -- print-comb/<Name>` green; `npm run lint:design` clean (proves the new exempt dir is correctly registered and nothing leaks appearance outside it); `npx tsc --noEmit -p tsconfig.build.json` clean. After the batch: `npm run build` exit 0 and `npm run build-storybook` exit 0.
