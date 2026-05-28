# "The Print" round-3 fidelity findings

**Date:** 2026-05-28 · validates `main` @ `4ae6073` (milestone #3 closed end-to-end — Phase 1 R0a–R10 + Phase 2 #207/#208/#209 all merged).
**Method:** live app + Playwright over all seven surfaces (samples / loupe-empty / loupe-rich / focus-empty / focus-indexed / folio / scoping / builder-waterfall / builder-heatmap / builder-tracking-on), with three parallel per-surface audit agents independently cross-checking shipped (screenshots + source) against `docs/redesign-mockups/*.html`, `DESIGN.md`, and `PRODUCT.md`.
**Inputs:** the round-2 findings doc `docs/2026-05-27-the-print-round2-findings.md` (used as a disposition checklist; closures verified at the screenshot *and* the source level, not just by reading PR descriptions).

---

## One-paragraph headline

The full Print is shipping. **Every workflow-breaking finding from round 2 is closed**, including the destination surface (`COMPARE_PALETTE` paper-tuned, single builder header, heatmap + cross-trace-tracking landed), the entry surface (contact-sheet fan-out fixed, permanent thumb chrome removed, loupe heading is finally Newsreader), and the focus surface (Plate Lift everywhere, breathing card-header, reflections-table panel built and wired to the q-link triple). What round 3 surfaces is **finish drift**: four P1 items, all of the same shape — a rule from `DESIGN.md` is being *misapplied*, not ignored. (i) Samples thumbnails are mouse-only because the thumb chrome was correctly removed but the `<div onClick>` root wasn't promoted to `<button>` (R3-S01). (ii) Focus auto-peak triangles paint in terracotta (R3-F01), turning the grease-pencil accent into a permanent decoration on every indexed sample. (iii) "+ New series" — the *named example* of `button-accent` in DESIGN.md §5 — renders as ink-solid instead (R3-Y01). (iv) `GroupingModeToggle`, sitting inline above the figure plate, still ships five dark-era token names (R3-Y02). Plus ~17 P2/P3 finish items mostly clustered around three patterns: legacy class-name carryovers in newly-touched files, Fixed-Scale Rule one-offs (`text-[25px]`, `text-[10px]`, `text-[9.5px]`) in new components, and mockup chrome under-implementation on the focus tools cluster + loupe sidebar meta. **The Print is now the unconditional shipped identity.** The R-residuals close it; they do not re-open it.

---

## Round-2 disposition — the inventory closed

| Surface | Closed | Partial / Open by design | Open / Regressed |
|---|---|---|---|
| Samples (R2-M11 thumb chrome, R2-M14 fan-out, R2-M15 cull-bar legend) | 3 | — | — |
| Loupe (R2-T9 heading, R2-M12 rep-box border, R2-T10 legacy tokens) | 3 | 1 (R2-M13 — meta-rows still Kind/Status, see R3-S02) | — |
| Focus (R3-N1 header crush, R3-N2 truncate, R3-N3 Plate Lift, R4-N1 dup copy) | 4 | 3 (R4-N2 Speculative · R4-N3 q-range chrome · R5-N1 Notes textarea — all verify-intentional in round 2, **judged unintentional in round 3**, see R3-F02/F03/F04) | — |
| Folio (R6-N1) | — | — | 1 (R3-Y05 — chips still silently zero) |
| Scoping (R7-N1, R7-N2) | 1 (R7-N2 closed-by-design: "the variable IS the name") | 1 (R7-N1 chevron — dev-DB empty path, source still has `text-xs` chevron at 15px parent) | — |
| Builder (R8-N1 dup header, R8-N2 palette retune) | 2 | — | — |

Net: **15 of 16 round-2 findings closed cleanly**; 1 finish-tier (R2-M13 loupe meta) still open and re-filed as R3-S02; 3 "verify intentional" items re-judged as unintentional this round (R3-F02/F03/F04). **No regressions.**

---

## Round-3 net-new findings — the residuals

**Severities:** P0 blocking · P1 major (accessibility violation or anti-pattern regression) · P2 finish · P3 polish.

### P1 — the four

- **R3-S01 (P1, a11y) — Contact-sheet thumbnails are not keyboard-operable.**
  **Where:** `packages/HimalayaUI/frontend/src/components/ContactSheetRow.tsx:75-92` (`ExposureThumb`'s root is `<div onClick onDoubleClick>` with no `role`, `tabIndex`, or `focus-visible`).
  **Why it matters:** the footer legend at `SamplesPage.tsx:199-217` advertises *click — select / ⇧ click — extend / double-click — open the loupe / Esc — clear / X — drop*. Only `Esc` and `X` are reachable without a mouse — a scientist on a tablet keyboard or VoiceOver cannot screen a sample. This is the entry surface; the workflow contract is broken for AT users. R2-M11 (strip overlay buttons) was the right move but it took the thumb's only keyboard-reachable focus targets with it.
  **Fix:** promote `ExposureThumb` root from `<div>` to `<button type="button">` (the onClick already fires on `Enter` for free); add `focus-visible:ring-2 focus-visible:ring-print-accent focus-visible:ring-offset-1`. Wire `onKeyDown` to `Enter`/`Space` for `onOpenLoupe`. ~30-line change.

- **R3-F01 (P1, anti-pattern) — Auto-peak triangles painted in terracotta, breaking accent rationing.**
  **Where:** `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx:410-413` (peak fill) + `packages/HimalayaUI/frontend/src/components/PlotCard.tsx:684` ("auto peak" legend swatch).
  **Why it matters:** DESIGN.md §2 names terracotta as **the grease pencil** — used on a small fraction of any screen for interaction. Every focus-indexed screen shows 5+ permanent terracotta peak triangles; the mockup paints peak markers in the **phase colour** of the index that claims them (`focus-workspace.html:740-753`), with `--accent` reserved for `.pk.hot` (the q-link cross-highlight). The shipped code's comment "Bright/neon for active workflow" is a dark-era hold-over. This is the single most load-bearing Print rule and it's the most-visible miss on the most-visited surface.
  **Fix:** pass `activeGroupIndices` + `claimOf` into peak-fill resolution; terracotta only when `peak.q === hoveredQ`. Unindexed peaks → `--color-ink-faint`. Remove the "auto peak" terracotta swatch from PlotLegend in the same pass.

- **R3-Y01 (P1, anti-pattern) — "+ New series" is `bg-ink` not `bg-print-accent`.**
  **Where:** `packages/HimalayaUI/frontend/src/pages/SeriesFolioPage.tsx:95`.
  **Why it matters:** DESIGN.md §5 lists this exact button as the canonical `button-accent` use ("the single grease-pencil action where the brand mark is the point (e.g. reject, '+ New series')"). The shipped ink fill makes it indistinguishable from a generic primary button and forfeits the one piece of brand colour the folio landing surface is supposed to carry. The mockup is unambiguous; the doc names this button by name.
  **Fix:** swap `bg-ink border-ink` → `bg-print-accent border-print-accent` (or hoist into a shared `<AccentButton>`).

- **R3-Y02 (P1, theming) — `GroupingModeToggle` ships five dark-era tokens; stale docstring.**
  **Where:** `packages/HimalayaUI/frontend/src/components/GroupingModeToggle.tsx:59-65` (classNames) + `:18-21` (block-comment still describes the dark vocabulary).
  **Why it matters:** sits inline above the figure plate; the active "By sample" pill is the most prominent above-the-fold control after the title. DESIGN.md §6 "Don't" explicitly retires `bg-bg-subtle` / `text-fg-dim` / `text-fg` / `bg-bg-hover` / `border-border`. Works today by R0a inheritance; any future cleanup of those aliases will silently break this control. Sibling toggles (`ScaleToggle`, `RepresentationToggle`) follow the canonical `bg-ink text-paper` active vocabulary.
  **Fix:** rewrite the classNames to mirror `ScaleToggle.tsx:27` / `:36`; delete the stale dark-era docstring.

### P2 — finish (organized by surface)

#### Samples + Loupe

- **R3-S02 (P2)** — Loupe sidebar meta-list still ships `Filename / Kind / Frame / Signal` (capitalised, app-internal); mockup ships `frame / integration / collected / signal` (lowercase, instrument-facing). "Kind" leaks SQLite schema nouns ("averaged" / "background_subtracted") at users who think in "2s exposure at 14:23". **Round-2 R2-M13 carryover.**
  **Where:** `LoupeSidebar.tsx:233-247`.
  **Fix:** replace `Filename` + `Kind` with `Integration` + `Collected`, sourced from exposure backend fields (stub to "—" until plumbed, matching the mockup's mock-data approach). Drop capitalisation to mockup style.

- **R3-S03 (P2)** — Verdict toggle button lost its mono `X` keycap. Mockup `v-toggle` reads `Drop X` / `Restore X` with the `X` set in mono — the same keycap idiom every other action on the page uses (CullBar's `Drop X`, footer-legend chips). Shipped reads bare `Drop` / `Restore`, so the keyboard shortcut goes undiscoverable from the right rail.
  **Where:** `LoupeSidebar.tsx:272-279`.
  **Fix:** append `<span className="ml-1 font-mono text-[10px] opacity-60">X</span>` to the button.

- **R3-S04 (P2, Fixed-Scale Rule)** — `text-[25px]` inline on the screened-progress numeral.
  **Where:** `SamplesPage.tsx:119` — `<div className="text-display !text-[25px] text-ink">`.
  **Fix:** add `--text-progress-numeral: 25px` + `.text-progress-numeral` role in `styles.css`, or reuse the existing `text-headline-lg` at 26px (1px delta invisible at this size).

#### Focus

- **R3-F02 (P2)** — "Speculative" section ships above-the-fold despite not being in the mockup. The rail's job is the calm "output" — every section above the fold that isn't in the mockup competes with the Phase call for the eye. Round-2 R4-N2 marked "verify intentional"; round 3 judges *unintentional*.
  **Where:** `PhasePanel.tsx:353-378`.
  **Fix:** move Speculative below the candidates note paragraph, or behind a `<details>` disclosure; only render the "+ Add" CTA on hover/focus or when at least one speculative exists.

- **R3-F03 (P2, anti-reference)** — Tools cluster ships QRange numeric inputs + reset instead of the mockup's `Auto-fit` + `+ Peak`. Reads as "legacy scientific software" — the number-input pair with manual reset is exactly the toolbar-soup pattern PRODUCT.md anti-references. Round-2 R4-N3 marked "verify intentional"; round 3 judges *unintentional*. Also: this cluster is where ~half of the surviving dark-era class names live (R3-F05).
  **Where:** `PlotCard.tsx:474-496`.
  **Fix:** replace QRange with `Auto-fit` ghost (calls existing `fitFeatures`) + `+ Peak` ghost (enters add-peak mode). Move numeric q-range edit into a popover off the segmented control if the affordance must survive.

- **R3-F04 (P2, anti-reference)** — Notes margin is a single textarea; missing threaded notes + count badge + mono q-refs. Mockup notes margin is the warmest piece of the focus surface (author initials in terracotta, threaded entries with timestamps, mono `q ≈ 0.064` references that mirror the q-link triple). Bare textarea is the "bare academic utility" anti-reference. Round-2 R5-N1 marked "verify intentional"; round 3 judges *unintentional*, at least for the visual treatment.
  **Where:** `FocusNotesMargin.tsx:24-51`.
  **Fix:** stage 1 (no schema change) — render `sample.notes` as quiet body text with a `✎` accent-prefixed "Add a note…" line à la mockup `.nm-add`; show count `1` if non-empty. Stage 2 (schema) — threaded notes.

- **R3-F05 (P2, theming)** — Six survivor dark-era token class names on the focus surface (mirrors R2-T10's loupe sweep). Work today by R0a inheritance; relics nonetheless.
  **Where:** `PlotCard.tsx:448-471, 481-484, 514-518, 559, 565, 577-579, 622-624, 682-683, 694`.
  **Fix:** mechanical sed: `text-fg-dim` → `text-ink-faint`, `text-fg-muted` → `text-ink-soft`, `text-fg` → `text-ink`, `bg-bg-hover` → `bg-paper-sunk`, `bg-bg-subtle` → `bg-paper-sunk`, `border-border` → `border-hair-strong`, `bg-bg` → `bg-plate`.

#### Series

- **R3-Y03 (P2)** — Rail "COMPOSE" header is louder than the figure-plate kicker. Rail uses `text-xs font-semibold uppercase tracking-wide text-ink` with a chevron, drawing the eye to the rail head before the plate kicker (terracotta, 11px, tracking 0.14em). Rail should recess behind the plate, not compete (§4 Flat-Except-the-Plate Rule).
  **Where:** `SeriesBuilderRail.tsx:78-89`.
  **Fix:** drop `text-ink` to `text-ink-faint`, match the kicker letter-spacing (`tracking-[0.14em]`), remove `font-semibold`.

- **R3-Y04 (P2, theming)** — "Track reflections" checkbox is a native HTML checkbox; renders OS-default blue in the rail (visible in `08-builder-tracking-on.png`).
  **Where:** `SeriesBuilderRail.tsx:131-137`.
  **Fix:** `accent-print-accent` on the input (mirrors `OffsetSlider`'s slider treatment), or build a custom-styled checkbox keyed to the accent.

- **R3-Y05 (P2)** — Folio `Has transition` / `Cross-experiment` chips silently filter to empty. The "No series match" copy implies there IS a search to clear when there isn't.  **Round-2 R6-N1 carryover.**
  **Where:** `lib/series/folioFilter.ts:42`; `SeriesFolioPage.tsx:168-171`.
  **Fix:** when `cross` is selected without source data, render a one-line explanation ("Cross-experiment matching not yet wired"), OR disable + tooltip the chip until the join exists.

- **R3-Y06 (P2)** — Heatmap rows have no outer hair keyline; reads as "two coloured boxes" instead of "two intensity rows in a frame". Mockup `series-builder.html:791-792` ships cells inside a keyline.
  **Where:** `MemberHeatmapLayer.tsx:113-115, 155-165` (`Plot.rect` emission has no `stroke`).
  **Fix:** emit a per-row `Plot.rect` outer keyline at `--color-hair`, OR reduce `rowInset` to 0.5 and let cells fuse vertically.

- **R3-Y07 (P2)** — Heatmap representation lacks the rotated ordering-variable axis label the mockup ships in the left margin. Even at N=5+ the heatmap will lack the "ordered by X" axis that grounds the migration story.
  **Where:** `MemberHeatmapLayer.tsx` has no axis; `SeriesBuilderPage.tsx:356-381` provides no rotated axis label in heatmap mode.
  **Fix:** in heatmap mode, render a rotated y-axis label using `s.ordering_variable` along the plot's left margin (mockup `series-builder.html:817-822` `.axis-title`).

### P3 — polish

| ID | Surface | Headline | Where |
|---|---|---|---|
| R3-S05 | Samples | "N dropped" kept-count label uses bare `text-print-accent`, overspending accent on a passive fact | `ContactSheetRow.tsx:378-381` |
| R3-S06 | Loupe | Last legacy `text-fg-muted` survivor (the placeholder "No image") | `DetectorImage.tsx:197` |
| R3-S07 | Loupe | Tag-form inputs lack focus indicator (`bg-transparent outline-none`); mockup wants 1px accent ring | `LoupeSidebar.tsx:131-167` + `ContactSheetRow.tsx:413-457` |
| R3-S08 | Samples | `+ tag` is `opacity-0 group-hover:opacity-100`; keyboard users cannot reach it | `ContactSheetRow.tsx:472` |
| R3-F06 | Focus | Reflections column header label "order" diverges from mockup "hkl"; not documented in DESIGN.md §5 | `FocusReflectionsTable.tsx:116` |
| R3-F07 | Focus | Reflections panel header lacks the right-side "exposure" cluster the mockup's `.panel-head` shows | `FocusReflectionsTable.tsx:224-226` |
| R3-F08 | Focus | Inline `text-[10px]` (Fixed-Scale Rule) on FocusDetectorPanel's "Set rep" button | `FocusDetectorPanel.tsx:150-152` |
| R3-F09 | Focus | Inline `text-[9.5px]` (Fixed-Scale Rule) on Reflections sticky `<th>` | `FocusReflectionsTable.tsx:110-126` |
| R3-F10 | Focus | Reflections "unindexed" row dims only phase + dash cells; mockup dims the whole `<tr>` | `FocusReflectionsTable.tsx:169-174` |
| R3-F11 | Focus | Reflections row `onMouseLeave` `hoveredQ` clear may race with row→row sweeps; verify symmetric with ring overlay | `FocusReflectionsTable.tsx:154-160` |
| R3-Y08 | Series | `CrossTraceTrackingLayer` uses `phaseColor()` correctly but lacks a docstring note pinning the choice (future "fix" risk to terracotta) | `CrossTraceTrackingLayer.tsx:119` |
| R3-Y09 | Series | Builder rail descendant selector `[&_input]:bg-plate` catches slider thumbs; scope it | `SeriesBuilderRail.tsx:92` |
| R3-Y10 | Series | Scoping + builder plate shadows omit the `0 1px 0 rgba(255,255,255,.6) inset` highlight DESIGN.md §4 Plate Lift specifies | `SeriesScopingPage.tsx:236`; `SeriesBuilderPage.tsx:277` |
| R3-Y11 | Series | Scoping `Discard` link lacks the mockup's `padding: 7px 4px`; floats unanchored | `SeriesScopingPage.tsx:228` |

---

## Cross-surface patterns

Three patterns surfaced by *more than one agent independently*. These are the load-bearing ones because they predict where the next regression will land.

1. **Terracotta is being misapplied in both directions.**
   R3-F01 *over*-uses the accent (every auto-peak triangle paints terracotta as decoration). R3-Y01 *under*-uses it (the named example in DESIGN.md §5 — "+ New series" — renders as ink). Both findings point at the same root issue: the team is treating terracotta as a "default primary action colour" rather than as the grease-pencil interaction mark. Tightening the rule requires fixing both sides at once, otherwise the next reviewer will fix one half and double-down on the other.

2. **Dark-era class names re-enter the codebase whenever a new component is added or a busy file is touched.**
   R2-T10 cleaned the loupe; R3-F05 finds six survivors in PlotCard (touched by #249); R3-Y02 finds five in GroupingModeToggle (touched by #251); R3-S06 finds one in DetectorImage (touched by #250). Each of the three PRs cleared the residuals in *its own* focal files but added or carried forward classes in adjacent files. **Pattern fix:** a project-wide grep `(bg-bg|text-fg|border-border|bg-bg-subtle|bg-bg-hover)` in CI, or a Stylelint custom rule. The R0a remap is doing the right thing at the value level; the class names are the leak.

3. **Fixed-Scale Rule one-offs recurred in three places in the new code.**
   `text-[25px]` (R3-S04, samples), `text-[10px]` (R3-F08, focus), `text-[9.5px]` (R3-F09, focus). R2-T9's fix added `text-headline-lg` to the scale — which was the right move — but the *anti-pattern* of "I need a one-off size, I'll just inline" did not get caught. **Pattern fix:** mention §3 in PR templates; consider a Stylelint deny-list on `text-\[\d+px\]`.

A fourth, narrower pattern worth naming: **mockup chrome detail is being under-implemented on controls.** Loupe Drop button missing its X keycap (R3-S03), loupe meta-rows shipping app-internal nouns (R3-S02), focus tools cluster shipping numeric inputs instead of Auto-fit/+Peak (R3-F03), focus Notes shipping a textarea (R3-F04). Each is a place where the implementation chose a "good enough" generic control over the mockup's specific affordance. The mockups are the spec; this is the gap between "feature shipped" and "feature shipped *the way the mockup said*".

---

## 5-dimension audit health score

| # | Dimension | Score | Key finding |
|---|-----------|-------|-------------|
| 1 | Accessibility | **2** | Samples thumbs are mouse-only (R3-S01, P1); loupe tag inputs lack focus indicators (R3-S07); native checkbox in builder rail (R3-Y04). Phase palette AA holds; semantic-dot pattern in focus is correct. |
| 2 | Performance | **3** | Round-2 fan-out closed (R2-M14); `IntersectionObserver`-gated thumbs; stable callback wiring on focus q-link; batched per-series queries on builder. No new perf debt. |
| 3 | Theming | **3** | All paper-tuned values land. Loses a point for the legacy class-name pattern across PlotCard / GroupingModeToggle / DetectorImage (R3-F05, R3-Y02, R3-S06) — works by inheritance, not by name. |
| 4 | Responsive | **3** | Folio columns step properly, builder plate widens 1180→1336px, loupe `max-w-[1100px]`. Untested below 900px; no mobile evidence in screenshots. |
| 5 | Anti-patterns | **2** | Two anti-pattern regressions surface as P1 (R3-F01 accent over-use; R3-Y01 accent under-use), one as P2 (R3-F03 QRange = legacy toolbar; R3-F04 textarea = bare academic). Identity holds firmly at the macro level (paper field, dark detector windows, single elevated plate, no card grids); the misses are *inside* the rules. |
| **Total** | | **13/20** | **Acceptable** — finish work needed, no foundational rework |

---

---

## Operator-added findings (2026-05-28, post-agent-audit)

Five additional findings raised by direct operator review of the running app; the per-surface agents structurally could not have surfaced them (audit screenshots can't show first-paint timing, missing detector data variance, or the LUT's hue choice — only what's *visible right now*).

### U-1 (P1, theming) — Detector window LUT is theme-agnostic neutral gray.

**Where:** `packages/HimalayaUI/src/image.jl:39-64` — `load_and_lognormalize` returns `Matrix{Gray{Float32}}` (pure neutral grayscale, no chroma).
**Why it matters:** `DESIGN.md §2` ships `frame-edge` (oklch 0.150 0.010 55) as the warm near-black for the *window*, but the image *inside* the window is theme-agnostic — neutral gray sits cold against the warm paper field. The dark windows have no warmth bridging them to the surrounding paper. Visible on every contact-sheet row and every loupe/focus big frame.
**Fix:** warm the LUT — either (a) a warm-tinted grayscale (multiply normed by a warm tint, oklch 0.x 0.005 60) for a subtle "ink on dark paper" effect, or (b) a true paper-tone heat LUT (warm yellow → terracotta → frame-edge) for a more aggressive Print restatement. (a) is the lower-risk first move.

### U-2 (P2, theming + perf) — Loading skeletons and missing-image placeholder don't fit The Print.

**Where:** boneyard skeleton fixtures across multiple surfaces; `DetectorImage.tsx:193-200` for "No image".
**Why it matters:** skeletons render as black blocks on paper — dark-era leftover, the very thing the milestone retired. Missing-image placeholder uses legacy `text-fg-muted` on the wrong surface (sits on `paper` instead of inside a `frame-edge` window). Empty/loading states are where The Print is at its most exposed because they ship the most chrome and the least data.
**Fix:** skeletons → `paper-sunk` shimmer on non-detector surfaces, `frame-edge` shimmer for detector-image skeletons specifically. Missing-image placeholder → render inside a `frame-edge` window with `frame-tag` caption text, matching the live detector treatment.

### U-3 (P2, consistency) — Contact-sheet thumbnails rotate orientation row-to-row.

**Where:** `packages/HimalayaUI/frontend/src/components/DetectorImage.tsx:184-210` — size-driven rotation logic (`transform: rotate(90deg)` when container is much wider than image).
**Why it matters:** the right call for the loupe big frame (where rotation maximizes the data area on a wide canvas), but it makes contact-sheet thumbs inconsistent — some thumbs portrait, others landscape-then-rotated, depending on container width. The contact sheet's job is to read as a uniform grid of windows so the *content* of each window is what the eye compares.
**Fix:** gate the rotation logic on `size === "full"`; for `size === "thumb"` always render portrait. Thumb aspect-ratio differences resolve via `object-fit: contain` inside the fixed-aspect window.

### U-4 (P2, focus chrome) — No indicator when trace is zoomed off-default.

**Where:** `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx` (zoom state) + `PlotCard.tsx` TitleStrip.
**Why it matters:** the trace plate is the single most-consulted surface on the focus workspace; when the user has zoomed in on a peak region, no chrome reflects "you're not seeing the whole curve". Subtle but corrosive — leads to false confidence that nothing else is happening at higher q. Adjacent to R3-F03 (replace QRange with Auto-fit/+Peak); the zoom-indicator + Auto-fit together turn the tools cluster from "manual q-range editing" into "the chart knows it's zoomed and offers to reset".
**Fix:** small terracotta or `text-meta` indicator in the TitleStrip when `xRange !== fullRange` ("zoomed · reset"); make it a ghost button that calls `fitFeatures`.

### U-5 (P1, perf) — Thumb pipeline recomputes lognormalize on every cold request.

**Where:** `packages/HimalayaUI/src/image.jl:39-64` + `routes_exposures.jl:48-55`.
**Why it matters:** the *plumbing* is correct — every thumb caller passes `size="thumb"` → `?thumb=1` → server resizes to 128px. But the math runs at full TIFF resolution: log1p over ~4M pixels, sort of a 4M-element vector for the percentile clip, clamp+normalize over 4M pixels — *then* resize to 128² = 16K pixels. **A ~250× per-pixel waste on every cold visit, multiplied by ~139 thumbs on a fresh contact-sheet load.** Visible as the staggered fill the user notices on cold reloads. Browser cache (`Cache-Control: immutable` + `?v=<token>`) handles repeat visits, but cold loads pay full cost every time. Server cold-start (Julia/Oxygen.jl warming) is a separate one-time hit; this is per-request.
**Fix (three composable, all in Issue H below):** (1) downscale TIFF → 128px *before* lognormalize for thumb requests — order-of-magnitude per-request speedup, every visit. (2) disk-backed PNG cache keyed on `image_version_token` — collapses any future cold-load to a file read. (3) **mandatory** thumb pre-warm in `init` and `reingest`, thread-parallel (`Threads.@threads`), skip-with-log on missing TIFF — makes the contact sheet fast on the very first visit to a freshly-ingested corpus.

### Also addressed by Issue G (not a separate finding)

`size="full"` requests (loupe big frame, focus big detector) currently send the unresized native-resolution TIFF→PNG (multi-MB for a 2048² detector) to fit a ~600-800px canvas. Capping the full path at ~1024-1536px max-side trims bandwidth with no visible quality loss. Pairs with U-5 since both touch `image.jl` / `routes_exposures.jl`.

---

## Recommended actions — final issue tree

This round-3 batch is more uniform than round-2 — clear per-surface clusters, three perf / theming items that share `image.jl`, no schema-affecting work, no architectural changes. Recommendation: **new milestone "HimalayaUI — The Print finish, round 3"** (milestone #4); 8 issues; ~5-6 small PRs (B + G + H could share PRs since they all touch `image.jl` / detector-window code).

| Issue | Title | Findings rolled in |
|---|---|---|
| **A** | Accent rationing — both directions | R3-F01 (peaks → phase color) + R3-Y01 ("+ New series" → `button-accent`) |
| **B** | Detector window warmth — LUT, skeleton, missing-image | U-1 (paper-tone LUT) + U-2 (Print skeleton + placeholder) + R3-S06 (legacy token cleanup falls out) |
| **C** | Contact-sheet keyboard + chrome cleanup | R3-S01 (a11y P1) + R3-S02/S03/S04/S05/S07/S08 + U-3 (lock thumbs portrait) |
| **D** | Focus chrome alignment | R3-F02/F03/F04 + R3-F06/F07/F09/F10 + U-4 (zoom indicator) |
| **E** | Series finish | R3-Y02 (token sweep P1) + R3-Y03/Y04 + R3-Y06/Y07 (heatmap keyline + axis) + R3-Y08–Y11 polish |
| **F** | Token-name sweep + CI grep guard | R3-F05 + cross-cutting dark-era class names + grep guard against `bg-bg|text-fg|border-border` re-introduction |
| **G** | Cap full-image PNG bytes | `LoupeFrame` + `FocusDetectorPanel` full path resized to ~1024-1536px max-side server-side |
| **H** | Thumbnail render perf | U-5 — downscale-before-lognormalize + disk-backed PNG cache + mandatory `init`/`reingest` pre-warm (thread-parallel) |

**Shelved by operator decision:** "eliminate monospace font" — contradicts DESIGN.md §3 Monospace-Means-Measurement Rule; needs separate design conversation before becoming an issue.

**Different track, not in this milestone:** beam-center detection — backend SAXS physics, belongs with the saxs-physics-reviewer's beat, not Print fidelity.

---

## Verdict

**The Print is the unconditional shipped identity.** Three rounds of validation in, the system holds at the macro level on every surface: paper field, dark detector windows, rationed terracotta accent in principle, Newsreader serif on display headlines, semantic phase palette at AA, single elevated plate, no card grids, no toolbar soup, no AI-slop tells at the page level. The residuals are *inside the rules*, not against them — places where a rule got *misapplied*, *under-applied*, or *carried forward from the previous era as a class-name relic*. None require re-opening the milestone; none change the architecture. The next pass is finish, not finish-of-the-finish-of-the-redesign.

End of round-3 audit.
