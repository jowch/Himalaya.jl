# "The Print" — Round-4 fidelity findings (post-milestone-4 close-out)

> **Status:** authored 2026-05-29 after milestone #4 ("HimalayaUI — The Print finish, round 3") closed all 8 issues (#254–#261) via PRs #262–#272. This is the fourth audit round, mirroring rounds 1–3 (live app + Playwright screenshots + parallel per-surface agents synthesized against `DESIGN.md`), but with the audit scored on the impeccable 5-dimension rubric and a deliberate **code-level verification pass** (live Julia image math reproduced numerically; keyboard paths traced through source). The deeper verification is why this round surfaces things screenshots alone could not — and why the headline score dips below round-3's 13/20 despite real fidelity progress.

## Headline

Milestone #4 landed most of what it promised, and landed it cleanly: the dark-era token system is **fully excised** (zero survivors), accent rationing is largely fixed (terracotta is correctly the single grease-pencil mark on every surface now), the focus zoom indicator, GroupingModeToggle Print tokens, full-image 1536px cap, and mandatory thumbnail pre-warm at ingest all shipped and verify against source. But the round exposes one **systemic gap the prior screenshot-only rounds missed — keyboard accessibility is broken across the app's core interactions** (frame culling, the q-link triple, rail phase-calls are all mouse-only or focus-invisible; WCAG 2.1.1 and 2.4.7 fail on three surfaces). Two milestone-4 items **shipped code that misses its goal** (the heatmap keyline is too low-contrast to frame anything; the R3-F CI guard was never built, so the token excision is unprotected), one perf optimization introduced a **real signal-fidelity regression** (thumbnail rings lose contrast under linear-space averaging), and the Fixed-Scale Rule leaks inline `text-[Npx]` on every surface. Score: **11/20 — Acceptable, significant work needed**, weighted down almost entirely by accessibility.

## Audit health score

| # | Dimension | Score | Key finding |
|---|-----------|-------|-------------|
| 1 | Accessibility | **1/4** | Core tasks (cull, q-link, phase-call) keyboard-unreachable; no focus-visible on rail/reflections — WCAG 2.1.1 + 2.4.7 fail across 3 surfaces |
| 2 | Performance | **3/4** | Excellent infra (lazy-gate, disk cache, mandatory prewarm, 1536 cap); dinged by thumb math regression + unbounded cache + main-thread LUT |
| 3 | Theming | **2/4** | Dark-era tokens fully excised (great), but Fixed-Scale leaks everywhere + several "landed-but-ineffective" fidelity misses |
| 4 | Responsive | **2/4** | Contact sheet (densest surface) has a fixed ~62rem grid, no breakpoints, no table scroll container; focus/series are fine |
| 5 | Anti-Patterns | **3/4** | Strong: no SaaS grid, terracotta rationed, dark-window held, ink-fill segments. Dinged by a self-defeating test + dead phase-chip code reported as shipped |
| **Total** | | **11/20** | **Acceptable (significant work needed)** |

## Anti-patterns verdict

**Pass — this does not look AI-generated.** Across all seven surfaces the design is distinctive and intentional: warm paper field, ink + hairline chrome, Newsreader serif titles, the folio as a wall of distinct mini-figures (not a card grid), detector images as dark windows, terracotta correctly rationed to interaction marks. No gradient text, no glassmorphism, no hero-metric tiles, no identical card grids, no bounce easing. The two anti-pattern *process* tells are non-visual: (1) a regression test (`test_image.jl`) written to **avoid** the failure mode it nominally covers, and (2) dead phase-chip code reported as "shipped" on a surface where it can never render. Both are honesty-of-coverage issues, not visual slop.

## Milestone-4 disposition

| Issue | Item | Verdict | Note |
|---|---|---|---|
| #254 R3-A | Peaks → phase color; "+New series" → terracotta | **Landed clean** (trace + table + folio) / **incomplete** (rings) | Focus big-detector *rings* still gray, not phase color — see P1-D1 |
| #255 R3-B | Paper-tone LUT + skeleton + missing-image placeholder | **Landed** (LUT, placeholder) / **partial** (skeleton) | LUT is client-side, applied equally to both paths; per-detector-image skeleton missing — see P2-CS3 |
| #256 R3-C | Contact-sheet keyboard a11y + portrait-lock | **Portrait-lock clean** / **a11y half-wired** | Keyboard users still cannot select/cull — see P1-CS1 |
| #257 R3-D | Focus chrome + zoom indicator + reflections | **Landed clean** | Zoom indicator verified (`PlotCard.tsx:567`) |
| #258 R3-E | GroupingModeToggle tokens + rail + heatmap | **Mostly clean** / **2 ineffective** | Keyline too low-contrast; Track-reflections checkbox checked-only — see P2-S1/S2 |
| #259 R3-F | Token sweep + CI grep guard | **Sweep clean** / **guard MISSING** | Zero dark-era survivors, but no `.github/` and no lint guard exists — see P1-F1 |
| #260 R3-G | Cap full-image PNG at ~1536px | **Landed clean** | Verified 854KB@1536px live + `test_routes_image.jl` |
| #261 R3-H | Thumb downscale + disk cache + prewarm | **Landed (4/4)** / **math bug + unbounded cache** | Prewarm mandatory in both entry points; downscale order regresses ring fidelity — see P1-IMG1, P2-IMG2 |

## P1 — Major (fix before next release)

### [P1-IMG1] Thumbnail downscale averages in linear-count space, degrading thin-ring fidelity
- **Dimension:** Performance (correctness-via-perf) · **Location:** `packages/HimalayaUI/src/image.jl:107-125` (`load_and_lognormalize_thumb`) vs `:48-66` (`load_and_lognormalize`)
- **Evidence:** The thumb path box-averages raw photon counts (`resize_to_fit` → `imresize`) **before** `log1p`, then percentile-clips on the downscaled raster; the full path takes `log1p` **first**. A 1–2px Bragg ring at count ~40 over background ~3 averages toward ~4 counts under 8× downscale, landing near the p5 low-clip floor, so the ring loses contrast against the field. Confirmed two ways: (a) numeric reproduction of `image.jl`'s math (same ring 1.0 full vs 0.287 thumb in a high-contrast synthetic); (b) **live measurement on the running app** — on a ringy frame (exp 3) the `%`-of-bright-pixels drops 70.6→58.5 thumb-vs-full, while mean brightness is only 5–20% lower.
- **Calibration (important):** the synthetic fixture overstates real-data impact. Live thumbs are **not** "crushed to near-pure-black" — they retain the beam, field, and gross intensity; what degrades is *thin-ring contrast / fine structure*. The "black thumbnails" perception on the contact sheet is partly this and partly the dark-window LUT rendering low-signal frames as dark squares at 128px (intended aesthetic colliding with triage legibility). Hence P1, not P0: the loupe (full path, unaffected) is one double-click away.
- **Impact:** The contact sheet's job is triage ("flip frames, drop the ones with flares or artifacts"). Losing thin-ring distinctness at a glance is a real degradation of that job, introduced by the #261 optimization.
- **Recommendation:** Log-space downscale — compute `lv = log1p(max(counts,0))` at full res, then `resize_to_fit(colorview(Gray, lv), max_px)`, then percentile-clip on the downscaled log raster (averaging in log space preserves ring contrast; the ~250× pixel win on the percentile sort is retained). Cheaper alt: compute p5/p99 on the full-res log raster, apply to the downscaled raster. Bump `IMAGE_PROCESSING_VERSION` to `v4` to invalidate cached PNGs. **Add a thin-ring regression test** — the current `test_image.jl` fixture deliberately uses a broad block and comments that a 1px ring "can legitimately wash out," so it structurally cannot catch this.

### [P1-CS1] Frame culling is keyboard-unreachable (WCAG 2.1.1)
- **Dimension:** Accessibility · **Location:** `ContactSheetRow.tsx:98-103`, `:258-278`
- **Evidence:** R3-C made each thumb a real `<button>` (good), but the keydown handler maps **both Enter and Space to `onOpenLoupe`** and `preventDefault`s the synthesized click, so `onSelect` never fires from the keyboard. The footer legend advertises "X — drop the selected frames" and "⇧click — extend the range," but the `X`/`Esc` batch handler only mounts once `hasSelection` is true, which a keyboard user can never reach. Culling is mouse-exclusive.
- **Impact:** WCAG 2.1.1 (Keyboard) failure on a core task. R3-C's a11y goal is half-delivered.
- **Recommendation:** Bind a key (Space = select / Enter = loupe, or `x` on a focused thumb) to `onSelect`; update the legend with keyboard equivalents.

### [P1-FL1] No visible focus on rail candidate rows (WCAG 2.4.7)
- **Dimension:** Accessibility · **Location:** `PhasePanel.tsx:93-121`
- **Evidence:** Candidate rows are `role="checkbox" tabIndex={0}` with correct Space/Enter handlers, but there is **no `focus-visible` ring/outline** — the inline border only changes on `inCall`/hover. Grep for `focus-visible|outline` in `PhasePanel.tsx`/`IndicesCard.tsx` returns nothing.
- **Impact:** Keyboard users get no focus state on the primary phase-call surface. DESIGN.md §5 designates the terracotta accent as the focus color; it isn't wired here.
- **Recommendation:** `focus-visible:outline focus-visible:outline-2 focus-visible:outline-print-accent focus-visible:outline-offset-2` on candidate rows and phase-call toggles.

### [P1-FL2] Reflections-table q-link cross-highlight is mouse-only
- **Dimension:** Accessibility · **Location:** `FocusReflectionsTable.tsx:142-166`
- **Evidence:** Rows are `<tr>` with `onMouseEnter`/`onMouseLeave` only — no `tabIndex`/`onFocus`/`onKeyDown`/`role`. The q-link triple (row ↔ peak ↔ ring), the headline feature of the focus surface, is unreachable by keyboard. The rings are `aria-hidden` and mouse-only too.
- **Impact:** A whole interaction channel is unavailable to keyboard/AT users.
- **Recommendation:** Make rows focusable (`tabIndex={0}` + `onFocus`/`onBlur` mirroring the mouse handlers driving `hoveredQ`), matching the rail-candidate pattern.

### [P1-D1] Focus big-detector rings are not phase-colored (R3-A incomplete)
- **Dimension:** Theming · **Location:** `DetectorRingOverlay.tsx:110`
- **Evidence:** `stroke={hot ? "var(--color-accent)" : "var(--color-ink-faint)"}` — verified directly. DESIGN.md §5 Detector frame: "rings render in phase colour on the focus big detector where the question is identity." R3-A recolored trace peaks and reflections-table dots but left the third q-link surface (rings) gray.
- **Impact:** The q-link triple is visually inconsistent — two of three surfaces carry phase identity, the rings don't, on the one surface the doc explicitly calls out for phase rings.
- **Recommendation:** Map `peak.q → claiming phaseColor` (same `claimColorByPeakId` logic as `TraceViewer.tsx:412`); keep terracotta for the `hot` ring, `ink-faint` for unclaimed.

### [P1-F1] R3-F CI grep guard does not exist — the token excision is unprotected
- **Dimension:** Tooling/CI · **Location:** repo root — no `.github/` directory exists at all
- **Evidence:** Verified directly: `ls .github` → absent; no `lint:tokens` script in `package.json`; no settings.json hook referencing the tokens. The sweep itself landed clean (0 dark-era survivors), but issue #259's "Done when" required a guard "that fails the build if they re-enter," and none was built.
- **Impact:** Nothing prevents `bg-bg`/`text-fg-dim`/`border-border` from silently re-entering on the next feature branch — exactly the cross-surface pattern round-3 flagged as systemic. The excision's durability depends on the guard that's missing.
- **Recommendation:** Add a `package.json` script wired into the existing pre-PR `npm run build` gate (the repo has no GitHub Actions): `grep -rnE` the token regex over `src/ --include='*.tsx' --include='*.ts' --include='*.css'`, fail non-zero on any hit. Anchor the regex so bare `text-fg`/`bg-bg`/`border-border` match alongside the `-dim/-muted/-subtle/-hover/-soft` suffixes.

### [P1-S1] Heatmap intensity encoded by hue-mix alone — no second channel, no scale
- **Dimension:** Accessibility/Theming · **Location:** `MemberHeatmapLayer.tsx:138-153`, `SeriesBuilderPage.tsx:371-384`
- **Evidence:** Cell fill is `color-mix(in oklab, <phaseHue> <pct>%, plate)` — intensity maps purely to mix-percentage of one hue, with no colorbar/legend and a per-row normalizer (so equal visual saturation encodes *different* absolute intensities across rows).
- **Impact:** Violates the Second-Channel Rule (DESIGN.md §2) — meaning rests on color saturation alone; unreadable under deuteranopia and gives no quantitative reference.
- **Recommendation:** Add an intensity scale legend (mono caption "low → high"), note per-row normalization in the fig-caption, and pair hue with the existing axis ticks; ideally a luminance ramp the eye can rank.

## P2 — Minor (next pass)

- **[P2-IMG2] Thumb disk cache is unbounded** — `image.jl` never prunes; every TIFF rewrite (new mtime → new token) or `IMAGE_PROCESSING_VERSION` bump orphans a PNG forever. `test_image.jl:137` even asserts the stale entry is "left in place (bounded)" — it is *un*bounded. Fix: in prewarm, glob `<id>-*.png` and delete entries not matching the current token. Low urgency for static corpora, matters for live-updating ones.
- **[P2-S1] Heatmap row keyline ships but is near-invisible** — `MemberHeatmapLayer.tsx:172-182` strokes `var(--color-hair)` (0.882) on `plate` (0.992), ~0.11 ΔL at 1px; combined with `rowInset=2` + the inter-member band gap, rows read as two floating boxes (the exact thing R3-Y06 set out to prevent). Fix: stroke `hair-strong` and/or draw one outer frame around the whole stack.
- **[P2-S2] "Track reflections" checkbox is native OS chrome when unchecked** — `SeriesBuilderRail.tsx:150-157` uses `accent-print-accent`, which tints only the *checked* fill; the unchecked box is the browser default (R3-Y04 half-shipped). Fix: `appearance-none` + explicit Print box (plate fill, hair border, 5px) + terracotta check glyph + accent focus ring.
- **[P2-CS3] No per-detector-image skeleton** — R3-B's skeleton is row/page-level only; the individual `DetectorImage` canvas is unpainted while fetch+decode+LUT runs (`DetectorImage.tsx:107-164`), so scrolling the 139-sample corpus shows a wave of blank windows. Fix: `frame-edge` fill + subtle shimmer until first `putImageData`.
- **[P2-CS4] Contact-sheet table has no responsive behavior** — fixed `grid-cols-[15.25rem_minmax(22.5rem,1fr)_4.875rem_10.5rem_9.375rem]` (~62rem min), no breakpoints, and only the inner exposure strip scrolls — the table container doesn't, so narrow viewports clip or force body h-scroll. Relevant to this repo's documented side-by-side dev workflow. Fix: wrap in an h-scroll container; ideally collapse Kept/Tags/Status at a breakpoint.
- **[P2-CS5] 556 detector canvases share one `aria-label="Detector image"`** — frame identity lives only on the parent button; the nested `role="img"` double-announces. Fix: `aria-hidden` the canvas on thumbs (button already names it).
- **[P2-FL3] Pervasive inline `text-[Npx]` (Fixed-Scale Rule)** — worst in `LoupeSidebar.tsx` (16 sites), also `PhasePanel.tsx`, `FocusReflectionsTable.tsx`, `SeriesBuilderPage.tsx`, `SampleStatusChip.tsx`, `ContactSheetRow.tsx`. Many (`text-[11.5px]`/`text-[10.5px]`) are literally token sizes; others (`text-[13px]`/`text-[9.5px]`/`text-[8.5px]`) are off-scale. Fix: route through `text-sm`/`text-caption`/`text-data-strong`; add reviewed scale steps for genuine off-scale needs.
- **[P2-FL4] Notes drawer uses a cool `shadow-xl`** — `FocusWorkspaceLayout.tsx:146`, a non-plate surface, violates Flat-Except-the-Plate with a generic Tailwind drop shadow (only visible `<xl`). Fix: use the warm Plate Lift shadow.
- **[P2-S3] Divergent segmented-control ARIA** — `GroupingModeToggle` uses `radiogroup`/`radio`/`aria-checked`; the visually-identical sibling `ScaleToggle`/`RepresentationToggle` use `group`/`aria-pressed`. Pick one pattern app-wide (radiogroup fits mutually-exclusive segments).

## P3 — Polish

- **[P3-CS6] Phase-chip seam is dead code on the contact sheet** — `SampleStatusChip` renders a chip when `sample.phase` is truthy, but no corpus phase-rollup is wired, so every row reads "Not indexed" and the chip never renders. R3-A's "chips replace 'Not indexed'" is structurally present but unobservable. Track the rollup as the real deliverable, or label the seam deferred so it isn't reported as shipped.
- **[P3-CS7] Selection ring + rep dot both terracotta** — `ContactSheetRow.tsx:109-110, 136-138`; a persistent *state* badge (rep) competes with the rationed interaction accent. Consider an ink/labeled treatment for the rep, reserving terracotta for the transient selection + reject mark.
- **[P3-FL5] Notes textarea `outline-none` with no replacement** — `FocusNotesMargin.tsx:118`; add a focus-visible ring.
- **[P3-FL6] Off-scale spacing one-offs** — `FocusNotesMargin.tsx:77,79,105` (`gap-[15px]`, `px-[19px]`, `py-[22px]`); snap to the 8/16/34 spacing scale.
- **[P3-S4] Folio sort group unlabeled + inactive segments have no rest affordance** — `SeriesFolioPage.tsx:130-148`; wrap in `role="group" aria-label="Sort"`, give inactive segments a hover affordance.
- **[P3-S5] Folio "+New series" hover uses opacity, not a darker terracotta** — `SeriesFolioPage.tsx:95` (`hover:opacity-90` lets paper bleed through). Hover to a darker accent shade.
- **[P3-CS8] React Router v7 future-flag warnings** on every page (`v7_startTransition`, `v7_relativeSplatPath`) — benign dev noise; opt into the future flags to silence.

## Cross-surface systemic patterns

1. **Keyboard accessibility is the app's structural weak point.** Three independent surfaces fail WCAG keyboard requirements on their *core* task: contact-sheet culling (2.1.1), focus rail focus-visible (2.4.7), reflections q-link (no keyboard path). The pattern is consistent — interactions are built mouse-first and the keyboard path is either absent or routes to the wrong handler. This is the single highest-leverage fix area and the reason A11y scores 1/4. It was invisible to rounds 1–3 because those audited screenshots, not interaction code.
2. **The Fixed-Scale Rule leaks on every surface.** Inline `text-[Npx]` appears in contact-sheet, focus/loupe (16+ sites in LoupeSidebar alone), and series components. The type scale is authoritative in DESIGN.md but not enforced — the same class of leak the R3-F guard was meant to catch for *colors* applies to *type sizes* with no guard at all.
3. **"Landed but ineffective / unprotected."** Several milestone-4 items shipped code that misses the goal: the heatmap keyline (too low-contrast to frame), the Track-reflections checkbox (checked-state only), the focus rings (recolor missed the third surface), and most starkly the R3-F CI guard (never built). The milestone closed issues whose "Done when" was only partially met. A "Done when = verified in the live UI/CI, not just merged" gate would catch this.

## Positive findings (maintain + replicate)

- **Dark-era token excision is complete and verified** — zero survivors across `frontend/src/`. The R0a→R3-F arc worked.
- **Accent rationing is now correct in both directions** — the round-3 P1 (over-rationed on focus, under-rationed on folio) is resolved; terracotta is the single grease-pencil mark on every surface. Confirmed `var(--color-accent)` only on hot q-link states, armed `+Peak`, zoom indicator, reject ✕, selection, and the "+New series" CTA. (The one remaining gap is *under*-application on focus rings — P1-D1.)
- **Image perf infrastructure is genuinely strong** — IntersectionObserver lazy-gating, atomic disk cache keyed on a version token resilient to source-delete, mandatory thread-parallel prewarm wired into both `init` and `reingest`, the 1536px full-image cap (verified 854KB live). The #261/#260 engineering is excellent; only the downscale color-space order and cache GC need work.
- **The detector dark-window rule holds** — both focus and loupe use `bg-frame-edge` with no Plate Lift; the warm LUT ramp (`frame-edge → frame-signal`) is correct and applied equally to both image paths.
- **Console is clean** — no runtime errors on any surface; only benign React Router v7 deprecation nudges.
- **Anti-pattern hygiene is high** — folio is a wall of distinct figures (not a card grid), segments are ink-fill, serif/sans/mono roles are honored, no gradient/glass/hero-metric anywhere.

## Recommended next steps

In priority order (the rubric's command mapping; all are `/impeccable` sub-commands or plain fixes):

1. **[P1] `harden` — keyboard accessibility sweep.** The highest-leverage work: wire keyboard culling (CS1), rail focus-visible (FL1), reflections q-link keyboard path (FL2). One coherent a11y pass across contact-sheet + focus.
2. **[P1] thumbnail log-space downscale + thin-ring regression test (IMG1)** + bump `IMAGE_PROCESSING_VERSION` to v4. Backend/`optimize`.
3. **[P1] build the R3-F CI guard (F1)** — close the issue that was marked done. Plus the type-scale leak (FL3) deserves the same guard treatment.
4. **[P1] `colorize` — finish the q-link triple (D1)** phase-color the focus rings; **(S1)** add the heatmap intensity legend + second channel.
5. **[P2] series + contact-sheet finish** — heatmap keyline contrast (S1), checkbox styling (S2), per-thumb skeleton (CS3), responsive table (CS4), cache GC (IMG2).
6. **`polish`** — the P3 batch (dead phase-chip seam, rep-dot accent, focus outlines, spacing, ARIA consistency) as a final pass.

> Re-run `/impeccable audit` after fixes to see the score climb — A11y is the dimension with the most headroom (1 → 3+ is one focused sweep).
