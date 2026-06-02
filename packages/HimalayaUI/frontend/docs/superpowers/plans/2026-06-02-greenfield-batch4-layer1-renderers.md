# Greenfield Batch 4 — Close Layer 1 (the three remaining renderers)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development (fresh implementer per task + spec review + quality review). Steps use `- [ ]` tracking.
> Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Commit ONLY the named files per task. NEVER `git add -A`. Do not touch pre-existing modified files (`src/bones/registry.ts`, `src/bones/contact-sheet.bones.json`) or the untracked detector-renderer plan doc.
> Typecheck with `npx tsc --noEmit -p tsconfig.build.json` (the default root tsc has unrelated errors).
> Branch `worktree-greenfield-ui-rebuild` stays UNMERGED — do NOT push / PR / finish.

**Goal:** Build the last three Layer-1 SVG renderers — `CardFigure` (mini-waterfall), `CustomPreview` (predicted-vs-observed comb), `CleanFigure` (plain export idiom) — completing Layer 1 and unblocking `SeriesCard`/`CustomIndexModal`/`ExportSheet`.

**Architecture:** Each renderer reuses appearance-free math (`makeAxis`/`makeProjection` from `src/print/plot/projection.ts`) for log-q scaling. `CardFigure` *also* reuses the branded `TracePlot` engine per-row (it wants the Print look) and lives beside `WaterfallChart`. `CustomPreview` hand-rolls one comb row and lives in the (already guard-exempt) `comb/` dir, sourcing phase color via `phaseColor()`. `CleanFigure` deliberately sheds the Print look (Arial + literal hex, no branding), so it reuses ONLY the math, hand-rolls its SVG, and lives in a NEW guard-exempt dir `src/print/export/` (the guard's `isExcluded` allowlist gains one prefix).

**Tech Stack:** React + TypeScript, d3-scale (via `makeAxis`), Vitest + Testing Library, Storybook (`@storybook/react-vite`). Tests assert via `data-*`/roles/text — NEVER class strings (`test/AGENTS.md`).

---

## Verified API facts (cite these; already checked against live source)

- `phaseColor(phase: string): string` — `src/phases.ts:45`. Returns an `oklch(...)` string at runtime (invisible to the line-based guard).
- `makeAxis(domain: [number,number], range: [number,number], type: "log"|"linear"): Axis1D` — `projection.ts:18`. `Axis1D.to(v)`, `.invert(px)`, `.ticks(count?)`, `.domain`, `.range`, `.type`.
- `axisTicks(axis: Axis1D, count=6): {value:number, kind:"major"|"mid"|"minor"}[]` — `projection.ts:50` (decade-anchored log ticks).
- `makeProjection({xDomain,yDomain,plotWidth,plotHeight,xType?,yType?}): {x:Axis1D,y:Axis1D}` — `projection.ts:110` (y is inverted: domain max → px 0).
- `positiveExtent(values: number[], fallback?): [number,number]` — `projection.ts:137`.
- `WaterfallRow` — `waterfall/waterfallModel.ts:17`: `{ key:string; label:string; trace:Trace; phase:string|null; anchors:WaterfallAnchor[]; ... }`. `Trace = { q:number[]; I:number[]; sigma:number[] }`. `WaterfallAnchor` (`:6`): `{ id; q:number; intensity; phase:string|null }`.
- Waterfall fixtures — `waterfall/waterfall.fixtures.ts`: `FULL` (5 rows), `TRANSITION` (3-row hero), `MIXED_STATES`.
- `CombSeries` — `comb/combModel.ts:17`: `{ phase:string; color:string; latticeLabel?:string; rSquared?:number; teeth:CombTooth[] }`. `CombTooth` (`:4`): `{ q:number; label:string; observed:boolean; residual?:number }`.
- Comb fixtures — `comb/comb.fixtures.ts`: `PN3M`, `IM3M` (both `CombSeries`), `LEFTOVER: number[]` (`[0.156, 0.205]`).
- `TracePlotProps` — `plot/TracePlot.tsx:37`: `{ trace:TraceModel; xDomain?; xType?; axes?:boolean; interaction?:TracePlotInteraction|false; height?; width?; show?:{peaks?;labels?;band?}; yHeadroom?; className?; ... }`. `TracePlot` + `TraceModel` exported from `plot/index.ts`.
- Design guard — `scripts/check-design.mjs:137` `function isExcluded(relPath)`: returns true when `relPath` (POSIX, relative to `src/`) `startsWith` one of `components/ui/`, `print/ui/`, `print/plot/`, `print/detector/`, `print/comb/`. `print/waterfall/` is NOT listed — `WaterfallChart` passes only by never writing a literal color (uses `phaseColor()` + `var(--color-*)`). Rule #5 flags literal `oklch(`/`rgba(`/quoted-hex anywhere on a line (after stripping `var(--color-*)`/`shadow-[…]`).
- Story idiom — `waterfall/WaterfallChart.stories.tsx`: `import type { Meta, StoryObj } from "@storybook/react-vite";` + loose `const meta: Meta<typeof X> = {…}` (NOT `satisfies` for render-only multi-row stories). `@storybook/test` is NOT installed — use a plain local `noop`, never `fn()`.
- JSDOM canonicalizes OKLCH on `style` write (`oklch(0.570 …)` → `oklch(0.57 …)`). Phase-color assertions must round-trip `phaseColor(phase)` through a reference element's `style` and compare read-backs — never `toBe(phaseColor(...))`.

---

## Task 1: `CardFigure` — frozen mini-waterfall (series-folio card)

**Files:**
- Create: `src/print/waterfall/CardFigure.tsx`
- Create: `src/print/waterfall/CardFigure.stories.tsx`
- Test: `test/print-waterfall/CardFigure.test.tsx` (new dir — mirror `test/print-plot/`)

**Mockup:** `../../../docs/redesign-mockups/series-folio.html` — the `SeriesCard` figure (`#wf-*`, CSS `.wf-base/.wf-fill/.wf-line/.wf-tick`). Frozen 340-wide stack, no axis, no label gutter, ~31px/row, dashed baselines, phase-colored filled curves + short peak ticks.

**Design intent:** A non-interactive miniature of `WaterfallChart` for embedding in a card. Reuse the branded `TracePlot` engine per row (it *should* look like The Print). Lives in `waterfall/` — passes the guard because it writes NO literal color (`phaseColor()` + tokens only).

**Spec / props:**
```tsx
import type { WaterfallRow } from "./waterfallModel";
export interface CardFigureProps {
  rows: WaterfallRow[];
  /** Overall SVG/figure width in px. Default 340 (the card figure width). */
  width?: number;
  /** Placement only. */
  className?: string;
}
```
Layout constants from the mockup: outer padding `padL=14, padR=14, padT=12, padB=12`; per-row band `stepY = 31`; plot width `pw = width - padL - padR`; total height `H = padT + rows.length*stepY + padB`. Shared log-q x-domain across all rows (use `positiveExtent` over every `row.trace.q`, then pad ~8% each side in log space, OR reuse the same domain helper `WaterfallChart` uses — read `WaterfallChart.tsx` and mirror its `qDomain` computation for consistency).

**Implementation approach (reuse, don't re-derive intensity math):** Stack rows bottom→top. Render each row as a FROZEN `TracePlot` mirroring how `WaterfallChart.tsx` builds its per-row plot, but mini and inert:
- `interaction={false}`, `axes={false}`, `show={{ peaks: true, labels: false, band: false }}`
- fixed `height={stepY}`, `width={pw}`, shared `xDomain` + `xType="log"`
- positioned at `y = padT + (rows.length-1-i)*stepY` (bottom→top) via absolute layout or nested `<g transform>` / stacked flex.
Read `WaterfallChart.tsx` for the exact `WaterfallRow.trace` → `TraceModel` conversion and per-row `phase`→color wiring, and copy that wiring (frozen). Do NOT invent a parallel intensity model.

**Required `data-*` for tests (no class assertions):**
- Root element: `data-testid="card-figure"`, `data-row-count={rows.length}`.
- One element per row carrying `data-row-key={row.key}` and `data-phase={row.phase ?? "none"}`.

- [ ] **Step 1: Write the failing test** — `test/print-waterfall/CardFigure.test.tsx`:
  - renders `<CardFigure rows={TRANSITION} />` (import `TRANSITION` from `../../src/print/waterfall/waterfall.fixtures`); asserts `getByTestId("card-figure")` has `data-row-count="3"`.
  - asserts exactly 3 elements with `[data-row-key]` exist, and their `data-phase` values match `TRANSITION.map(r => r.phase ?? "none")`.
  - phase-color round-trip: pick the first row; render a reference `<span style={{color: phaseColor(TRANSITION[0].phase!)}}/>`; assert the row's rendered color (read back from the element that carries the curve color, via `getComputedStyle` or the inline `style`) equals the reference read-back (NOT a literal compare).
- [ ] **Step 2: Run → fail.** `npm test -- CardFigure` → FAIL (module not found).
- [ ] **Step 3: Implement `CardFigure.tsx`** per the approach above.
- [ ] **Step 4: Write `CardFigure.stories.tsx`** — loose `Meta<typeof CardFigure>`; stories `Transition` (rows=`TRANSITION`), `Full` (rows=`FULL`), `Mixed` (rows=`MIXED_STATES`), each wrapped in `<div className="rounded border border-hair bg-plate p-4" style={{ width: 360 }}>`.
- [ ] **Step 5: Green.** `npm test -- CardFigure` PASS · `npm run lint:design` clean · `npx tsc --noEmit -p tsconfig.build.json` clean.
- [ ] **Step 6: Commit** `git add src/print/waterfall/CardFigure.tsx src/print/waterfall/CardFigure.stories.tsx test/print-waterfall/CardFigure.test.tsx` → `feat(print): CardFigure frozen mini-waterfall renderer`.

---

## Task 2: `CustomPreview` — predicted-vs-observed comb (custom-index modal)

**Files:**
- Create: `src/print/comb/CustomPreview.tsx`
- Create: `src/print/comb/CustomPreview.stories.tsx`
- Modify: `src/print/comb/index.ts` (add `export { CustomPreview } from "./CustomPreview";` + `export type { CustomPreviewProps }`)
- Test: `test/print-comb/CustomPreview.test.tsx` (new dir if absent)

**Mockup:** `../../../docs/redesign-mockups/2026-05-29-focus-plot.html` — `#cs-preview` (`viewBox="0 0 520 150"`, aria-label "Custom index reflections against observed peaks") inside the `.custom-sheet` modal. A single comb row: log-q baseline, observed peaks as small markers, predicted reflections as upright capped stems (phase-colored), predicted-absent reflections as faint dashed carets.

**Design intent:** A STATIC single-row comb showing one candidate phase's predicted teeth (at the user's chosen lattice) against the observed peaks. `CombChart`/`CombScaffold` are too heavy (gutter + scroll + interactivity) — hand-roll one row, but reuse `makeAxis` for log-q and `phaseColor()` for color. Lives in `comb/` (guard-exempt).

**Spec / props:**
```tsx
import type { CombSeries } from "./combModel";
export interface CustomPreviewProps {
  /** The candidate phase's predicted comb at the chosen lattice. */
  series: CombSeries;
  /** Observed peak q-values (Å⁻¹) to overlay. */
  observed: number[];
  width?: number;   // default 520
  height?: number;  // default 150
  className?: string;
}
```
Layout: margins `mL=10, mR=10, mT=14, mB=26`; `pw = width - mL - mR`; baseline `baseY = height - mB`; stem top `topY = mT`. Log-q domain over `[...series.teeth.map(t=>t.q), ...observed]` via `positiveExtent`, padded ~8% each side; `x = makeAxis(domain, [mL, mL+pw], "log")`.

Draw, per element (all colors via `phaseColor(series.phase)` or `var(--color-ink-faint)` token — NO literals):
- baseline: horizontal `<line>` at `baseY`, `stroke="var(--color-hair-strong)"`, sw 1.
- observed markers: for each `q` in `observed`, a small `<circle r=2.5>` at `(x.to(q), baseY)`, `fill="var(--color-ink-soft)"`. Carry `data-observed-q={q}`.
- predicted teeth: for each `t` in `series.teeth`:
  - `observed:true` → solid stem `<line>` from `(x.to(t.q), baseY)` to `(x.to(t.q), topY + k)` (`k` a small per-tooth offset is optional), `stroke={phaseColor(series.phase)}` via inline `style`, sw 1.8, round cap; plus a tiny `√n` label (`<text>`, `var(--font-mono)`, `fill="var(--color-ink-faint)"`).
  - `observed:false` → faint dashed caret/stem, `stroke="var(--color-ink-faint)"`, `stroke-dasharray="1.5 1.8"`, sw 1.5.
  - Carry `data-tooth-q={t.q}`, `data-tooth-observed={t.observed}`, `data-tooth-label={t.label}`.
- Root `<svg>`: `viewBox="0 0 {width} {height}"`, `role="img"`, `aria-label="Custom index reflections against observed peaks"`, `data-testid="custom-preview"`, `data-phase={series.phase}`.

**Phase color lives via inline `style`** (mirror `ScoreBar`/`Swatch`): `style={{ stroke: phaseColor(series.phase) }}` on the solid stems — never a literal in source.

- [ ] **Step 1: Failing test** — `test/print-comb/CustomPreview.test.tsx`:
  - renders `<CustomPreview series={PN3M} observed={[0.0712, 0.0872, 0.156]} />` (import `PN3M` from `../../src/print/comb/comb.fixtures`).
  - asserts `getByRole("img", { name: /custom index reflections/i })` exists and has `data-phase="Pn3m"`.
  - asserts one `[data-observed-q]` per observed value (3).
  - asserts one `[data-tooth-q]` per `PN3M.teeth` entry; asserts at least one `[data-tooth-observed="true"]` and (if `PN3M` has an absent tooth) one `[data-tooth-observed="false"]`.
  - phase-color round-trip on a solid stem vs a reference `style={{stroke: phaseColor("Pn3m")}}` element (read-back compare, not literal).
- [ ] **Step 2: Run → fail.** `npm test -- CustomPreview` → FAIL.
- [ ] **Step 3: Implement `CustomPreview.tsx`** + add the `index.ts` export.
- [ ] **Step 4: Stories** — loose `Meta`; `Default` (series=`PN3M`, observed from its observed teeth q's), `Cubic` (series=`IM3M`), `WithAbsent` (a `series` whose teeth include `observed:false` — build inline from `PN3M` with one tooth flipped). Wrap each in a `bg-plate` framed `<div>`.
- [ ] **Step 5: Green.** `npm test -- CustomPreview` · `npm run lint:design` · `npx tsc --noEmit -p tsconfig.build.json`.
- [ ] **Step 6: Commit** `git add src/print/comb/CustomPreview.tsx src/print/comb/CustomPreview.stories.tsx src/print/comb/index.ts test/print-comb/CustomPreview.test.tsx` → `feat(print): CustomPreview predicted-vs-observed comb renderer`.

---

## Task 3: `CleanFigure` — plain export idiom (+ design-guard exempt-dir)

**Files:**
- Modify: `scripts/check-design.mjs` (add `print/export/` to `isExcluded`)
- Create: `src/print/export/CleanFigure.tsx`
- Create: `src/print/export/index.ts` (`export { CleanFigure } from "./CleanFigure";` + `export type { CleanFigureProps }`)
- Create: `src/print/export/CleanFigure.stories.tsx`
- Test: `test/print-export/CleanFigure.test.tsx` (new dir)

**Mockup:** `../../../docs/redesign-mockups/2026-05-29-series-plot.html` — `.clean-card` + `#clean-fig` (`viewBox="0 0 600 340"`). A GraphPad-style figure: white card, centered bold Arial title, plain SVG plot with black L-shaped axes, Arial tick labels + bold axis titles, hex-colored trace lines (no fill), right-margin Arial sample labels, small gray Arial footer. **Deliberately un-branded** — this is what users export, so it must NOT use Print tokens.

**Why a guard edit:** `CleanFigure` writes literal hex (`#111`, `#c8841f`, `#4a4ba8`, `#666`, `#e3e3e3`) and `font-family: Arial`. Rule #5 of `check-design.mjs` flags literal hex anywhere. Existing renderer dirs (`plot/`/`comb/`/`detector/`) are exempted via `isExcluded`; `CleanFigure` gets the same treatment with a new `print/export/` prefix. This mirrors how the legacy `lib/figure-export/**` export palette is allowlisted.

**Spec / props:**
```tsx
import type { WaterfallRow } from "../waterfall/waterfallModel";
export interface CleanFigureProps {
  rows: WaterfallRow[];
  title: string;
  footer: string;
  width?: number;   // default 600
  height?: number;  // default 340
  className?: string;
}
```
Render the full `.clean-card`: a `<div>` (white bg `#ffffff`, `1px solid #e3e3e3`, radius 3px, padding 16/14/10) containing a centered bold Arial `#111` title, the `<svg>`, and a centered 9px Arial `#666` footer. Use inline `style` for all of this (NOT Tailwind) so the export stays self-contained and literal.

SVG: `viewBox="0 0 {width} {height}"`, margins `m={l:52,r:64,t:12,b:40}`; `pw=width-m.l-m.r`, `ph=height-m.t-m.b`; `baseY=m.t+ph`. Log-q x-axis via `makeAxis(domain,[m.l,m.l+pw],"log")` (domain = padded `positiveExtent` over all `row.trace.q`); ticks via `axisTicks(x)`.
- Axes: x-line and y-line, `stroke="#111"`, sw 1.6.
- X ticks: 5px black ticks (sw 1) + Arial size-9 `#111` labels (`q` formatted, e.g. `.toFixed(2)`); x-axis title `q (Å⁻¹)` Arial size-11 weight-700 `#111` centered under axis.
- Y-axis title `Intensity (a.u.)` Arial 11 bold `#111`, rotated −90, centered on the y-axis.
- Traces: per row (bottom→top, offset by `i*(ph/rows.length)` baseline), build a path `M…L…` over the row's `(q,I)` points mapped through `x.to(q)` and a per-row y-normalization (normalize `I` by the GLOBAL max across rows → amp); `stroke = HEX[row.phase] ?? "#777"`, sw 2, `fill: "none"`. Local palette: `const HEX: Record<string,string> = { Pn3m: "#c8841f", Lamellar: "#4a4ba8", Im3m: "#4a4ba8" }` (extend minimally to cover fixture phases; unknown → `#777`).
- Sample labels: right margin (`m.l+pw+6`), Arial 9 `#111`, one per row at its baseline.
- Root `<svg>`: `role="img"`, `aria-label={title}`.

**Required `data-*`:** card root `data-testid="clean-figure"`; title text rendered verbatim; one trace `<path>` per row carrying `data-phase={row.phase ?? "none"}`; footer text rendered verbatim.

- [ ] **Step 1: Guard edit first** — in `scripts/check-design.mjs` `isExcluded`, add after the `print/comb/` line:
  ```js
      relPath.startsWith("print/comb/") ||
      relPath.startsWith("print/export/")   // CleanFigure: deliberately un-branded export idiom (Arial + literal hex)
  ```
  Also update the header comment block (`:134`) listing exempt dirs to include `src/print/export/**`.
- [ ] **Step 2: Failing test** — `test/print-export/CleanFigure.test.tsx`:
  - renders `<CleanFigure rows={TRANSITION} title="LL37 Titration: Pn3m → Lamellar" footer="Lipid 1-2 · SSRL Apr 2026" />`.
  - asserts `getByTestId("clean-figure")` exists; `getByText("LL37 Titration: Pn3m → Lamellar")` and `getByText(/Lipid 1-2/)` present.
  - asserts `getByRole("img", { name: /LL37 Titration/ })` present.
  - asserts exactly `TRANSITION.length` `[data-phase]` trace paths.
  - (No phase-color round-trip needed — colors are literal hex; assert via `data-phase` only.)
- [ ] **Step 3: Run → fail.** `npm test -- CleanFigure` → FAIL.
- [ ] **Step 4: Implement** `src/print/export/CleanFigure.tsx` + `src/print/export/index.ts`.
- [ ] **Step 5: Stories** — loose `Meta`; `Default` (rows=`TRANSITION`, the mockup's title/footer), `FullSeries` (rows=`FULL`). No Print framing wrapper (it's a standalone card).
- [ ] **Step 6: Green.** `npm test -- CleanFigure` PASS · `npm run lint:design` clean (proves the exempt-dir edit works — if it FAILS on `print/export/`, the guard edit is wrong) · `npx tsc --noEmit -p tsconfig.build.json` clean.
- [ ] **Step 7: Commit** `git add scripts/check-design.mjs src/print/export/CleanFigure.tsx src/print/export/index.ts src/print/export/CleanFigure.stories.tsx test/print-export/CleanFigure.test.tsx` → `feat(print): CleanFigure plain export renderer + print/export guard exemption`.

---

## Task 4: Ledger + close-out

**Files:**
- Modify: `docs/greenfield-component-ledger.md`

- [ ] **Step 1:** Flip the three Layer-1 rows (`CardFigure`, `CustomPreview`, `CleanFigure`) ⬜→✅ with their file paths; update the Layer-1 summary line to `9 ✅ · 0 ⬜` and the coverage summary; record the `print/export/` exempt-dir addition in the Out-of-scope/decisions registry (so it isn't re-litigated); update the closing "build frontier" paragraph to point at the now-fully-unblocked Layer-3 tier-2 panels.
- [ ] **Step 2: Commit** `git add docs/greenfield-component-ledger.md` → `docs(print): ledger — Layer 1 complete (CardFigure/CustomPreview/CleanFigure)`.

## Final verification (after all tasks)

- `npm test` — full suite green (capture once; grep the file).
- `npx tsc --noEmit -p tsconfig.build.json` — exit 0.
- `npm run lint:design` — clean (proves `print/export/` exemption + literal-free `CardFigure`/`CustomPreview`).
- `npm run build-storybook` — exit 0.
- Visual fidelity: serve `storybook-static`, screenshot the three new renderer stories, compare against the cited mockup regions.
