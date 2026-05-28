# "The Print" round-2 fidelity findings

**Date:** 2026-05-27 · validates `main` @ `690929a` (milestone #3 R0a–R10 + R9 verify all merged)
**Method:** live app + Playwright over all six surfaces (samples / loupe / focus-empty / focus-indexed / folio / scoping-empty / scoping-proposed / builder), with three parallel per-stage audit agents cross-checking shipped (screenshots + source) against `docs/redesign-mockups/*.html` and `DESIGN.md`.
**Inputs:** the round-1 inventory `docs/2026-05-27-the-print-fidelity-findings.md` (used as a disposition checklist), the remediation plan `docs/superpowers/plans/2026-05-27-the-print-remediation-plan.md` (defines what each R-issue was supposed to ship).

---

## One-paragraph headline

R0a–R10 delivered the structural redesign. **"The Print" is the unconditional identity** — no surface still resolves to the old Darkroom defaults, no ice-blue chrome residue, no dark `theme-light` override. Newsreader is loaded and used on every display headline except one (R2-T9 below). The two biggest pieces of work — the keystone R0a token unification and R7's scoping proposal UX — both land. **Every workflow-breaking round-1 finding closed.** What round-2 surfaces are finish-level: two P2 regressions on the destination surface (the series builder), one P2 cluster on the contact-sheet's thumbnail chrome, and ~14 P3 polish items. The remaining work is one tight focused pass — call it R11, no second milestone.

---

## Round-1 disposition — the inventory closed

| Surface | Resolved | Partial (works by inheritance / data not plumbed) | Still-open / Regression | Deferred-as-planned |
|---|---|---|---|---|
| Samples (T-/M-/L-) | 11 | — | 2 (M-4/M-5 — #207) | — |
| Loupe (T-/M-/L-) | 7 | 4 (T-2/T-3/T-6/T-10 — legacy class names) | 1 (L-2) | — |
| Focus (A-/F-/L-) | 13 | 1 (L-11 — verify live) | — | — |
| Folio (F-/N-/X-) | 6 | 1 (F-C — `newMatches` plumbing cold) | — | — |
| Scoping (S-/L-) | 8 | — | 1 (S-G — editable working name) | — |
| Builder (B-/X-) | 8 | 3 (B-B/B-C — work by inheritance) | — | 2 (B-D/B-E — #208) |

Net round-1: ~53 findings, ~46 resolved cleanly, 8 partial-but-functional, **2 still-open of substance** (L-2 = R2-T9 below, S-G = R7-N2 below). Zero regressions.

---

## Round-2 net-new findings — the residuals

**Severities:** P2 degraded-but-functional · P3 polish.

### Builder — the destination surface

- **R8-N1 (P2)** — Duplicated "SERIES / Bundle A" header.
  **Where:** `packages/HimalayaUI/frontend/src/pages/SeriesBuilderPage.tsx:99-124` (outer page header) vs `:284-301` (plate kicker+title).
  **Why it matters:** the same kicker + title renders twice (~80px of vertical stack), diluting the figure-as-plate metaphor that R8 B-J shipped. Mockup `series-builder.html:386-396` has only the plate header.
  **Fix:** delete the outer `<header>`, re-home the `Edit` button into the plate kicker row or the rail head.

- **R8-N2 (P2)** — `COMPARE_PALETTE` not paper-tuned; bySample/distinct traces wash out on paper.
  **Where:** `packages/HimalayaUI/frontend/src/lib/comparison/coloring.ts:47-60` — palette L 0.76–0.80 (docstring still reads "dark-background contrast"), unchanged by R0b.
  **Why it matters:** the live builder defaults to `groupingMode="bySample"` and most series are unindexed, so the bySample trace coloring is the dominant visual on the destination surface. R0b paper-tuned `PHASE_PALETTE` to L 0.50–0.58; the parallel `COMPARE_PALETTE` was missed.
  **Fix:** paper-retune `COMPARE_PALETTE` to L 0.50–0.58 (preserve 12-hue rotation), then re-run `test/coloring.test.ts` no-collision invariant against `PHASE_PALETTE`. One-touch fix; immediate paper feel.

### Samples — the entry surface

- **R2-M11 (P2)** — Contact-sheet thumbnails ship three permanent overlay buttons not in the mockup.
  **Where:** `packages/HimalayaUI/frontend/src/components/ContactSheetRow.tsx:100-114` (selection checkbox), `:127-141` (rep ⊙), `:142-154` (reject ✕).
  **Why it matters:** paper chips bleed across every dark detector window, undermining the "dark window in light paper" metaphor and forcing 1,668 keyboard tab stops (R2-F1) past the contact sheet. Mockup `sample-table.html:203-231` has zero rest-state chrome — selection is a `box-shadow` driven by `.thumb.sel`, rejection is via the CullBar + `X` keystroke (already shipped), rep is a corner state badge on `.is-rep` (already shipped at `:117-124`).
  **Fix:** remove the three permanent buttons; rely on the existing `ring-print-accent` selection ring, the floating CullBar, and the `R`/`X` keystrokes. If keeping per-thumb ✕ for discoverability, `opacity-0 group-hover:opacity-100`.

- **R2-M14 (P2)** — Contact sheet only ever shows the first row of thumbnails.
  **Where:** `ContactSheetRow.tsx:179` (per-row `useExposures(sample.id)` fan-out) + `DetectorImage.tsx:86` (unbounded `fetch`).
  **Why it matters:** 139 samples × per-row `useExposures` + ~556 thumb image fetches blow past Chromium's 6-per-host limit (`ERR_INSUFFICIENT_RESOURCES` in the console). Rows 2-7 visible in `samples-viewport.png` are all stuck on "Loading frames…". Fidelity is moot on a surface that never paints.
  **Fix:** (a) one bulk `useCorpusExposures()` query returning `{sampleId → exposures[]}` so rows don't fan out, (b) `IntersectionObserver`-gate `DetectorImage` rendering so only on-screen thumbs fetch. Either alone unblocks the sheet; both together restore the "scan the whole corpus" promise.

### Loupe — sample focus

- **R2-T9 (P2)** — Loupe heading is `text-2xl` sans (24px Plus Jakarta Sans), not Newsreader serif. **Carries over round-1 L-2; R2 shipped without addressing it.**
  **Where:** `packages/HimalayaUI/frontend/src/pages/LoupePage.tsx:185`.
  **Why it matters:** the only display headline on these two surfaces that misses Newsreader. The contact sheet H1 "The corpus" is serif; the loupe H2 "HEPES Only" is not.
  **Fix:** use the canonical `text-headline` role (Newsreader 19px) or extend the type scale with `text-headline-lg` at 26px in `styles.css`. Do not inline `text-[26px]` (DESIGN.md §3 Fixed-Scale Rule).

- **R2-M12 (P3)** — Loupe rep-box missing active-state border treatment.
  **Where:** `packages/HimalayaUI/frontend/src/components/LoupeSidebar.tsx:156-176`.
  **Why it matters:** the rep-box looks identical whether the active exposure is the representative or not — only the inner dot lights. Mockup `.rep-box.is-rep` adds a terracotta-tinged border (`color-mix(in oklab, var(--accent) 45%, var(--hair))`).
  **Fix:** when `isRepresentative`, swap `border-border` for `border-print-accent/40`.

### Focus — sample workspace

- **R4-N1 (P2)** — Duplicate "Check every phase that is present" copy in the rail.
  **Where:** `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx:273-275` (card-header subtitle) and `:343-346` (rail-note paragraph).
  **Why it matters:** identical framing sentence repeated ~3 lines apart in the same rail. The rail's job is to read as the calm "output"; repetition undermines that.
  **Fix:** drop the header subtitle, keep the longer rail-note (or vice versa — mockup has the header *label* "Index choices" with no inline subtitle, and the explanatory paragraph only below the candidates).

### Polish-tier (P3)

| ID | Surface | Headline | Where |
|---|---|---|---|
| R2-T10 | Loupe | 7 legacy `border-border`/`bg-bg-subtle`/`text-accent` class names survive; work by R0a inheritance | `LoupeSidebar.tsx:121,144-145,154,170-171,189-190`; `LoupePage.tsx:161,166` |
| R2-M13 | Loupe | Sidebar meta rows are Kind+Status (app-internal); mockup ships Integration+Collected (instrument-facing) | `LoupeSidebar.tsx:97-110` |
| R2-F1 | Samples | Keyboard tab order: 4 stops per thumb because of R2-M11's overlay buttons | (tied to R2-M11) |
| R2-M15 | Samples | Cull-bar discoverability legend offscreen on long scrolls — sticky-bottom or fold into CullBar | `SamplesPage.tsx:179-195` |
| R3-N1 | Focus | `card-header` fixed 56px crushes the kicker+serif+subline 3-row stack; mockup `.plate-head` breathes | `PlotCard.tsx:417` + `styles.css:133-141` |
| R3-N2 | Focus | FocusPlotHeader subline `max-w-[60ch]` truncates the exposure label mid-token | `FocusPlotHeader.tsx:54-58` |
| R3-N3 | Focus | PlotCard outer missing the Plate Lift shadow; trace plate reads as a flat rectangle | `PlotCard.tsx:338` (no `.card`) |
| R4-N2 | Focus | "Speculative" section above-the-fold in rail (not in mockup) — verify intentional | `PhasePanel.tsx:349-375` |
| R4-N3 | Focus | q-range numeric inputs retained vs mockup's simpler Auto-fit/+Peak — verify intentional | `PlotCard.tsx:457-480` |
| R5-N1 | Focus | Notes margin missing count badge + per-note rendering (mockup ships threaded notes; shipped is single textarea) — verify intentional | `FocusNotesMargin.tsx:33-49` |
| R6-N1 | Folio | Filter chips "Has transition" / "Cross-experiment" silently zero-result; no source data wired | `SeriesFolioPage.tsx:155-184`; `SeriesFolioCard.tsx:60` |
| R7-N1 | Scoping | "Ordered by" chevron ▾ optically dropped — font-size mismatch with the 15px value text | `ScopingOrderField.tsx:21` |
| R7-N2 | Scoping | No editable serif working name on the plate (round-1 S-G carryover) — verify intentional ("the variable IS the name") | `SeriesScopingPage.tsx:241-243` |

### Production-quality (non-fidelity)

The contact-sheet thumbnail fan-out (**R2-M14**) is the only sustained perf finding. Worth pairing with R2-M11 as a "contact-sheet follow-up" issue rather than its own surface fix.

---

## Verdict per surface

| Surface | Verdict | Why |
|---|---|---|
| Shared chrome (CorpusShell + CorpusTopbar) | **Shipped clean** | R0a's keystone landed; The Print is the unconditional identity. |
| Samples (contact sheet) | **Shipped with residuals** | Structure is right; permanent thumb chrome (R2-M11) and per-row query fan-out (R2-M14) compromise the metaphor and the live-paint promise. |
| Loupe | **Shipped with residuals** | R2 deliverables present; one missed serif headline (R2-T9 = round-1 L-2 unrepaired), one active-state border (R2-M12), legacy token class names (R2-T10). |
| Focus | **Shipped with residuals** | R3/R4/R5 closed the workflow-breakers; residuals are finish-level (duplicate copy, card-header crush, subline truncate, missing plate shadow). |
| Folio (R6) | **Shipped clean** | All F-* + X-1 resolved; only cold-data chip plumbing (R6-N1). |
| Scoping (R7) | **Shipped clean** | The biggest single piece of work delivered fully; only P3 polish (chevron size, working-name question). |
| Builder (R8) | **Shipped with residuals** | Compose-rail controls all land, but the **duplicated header (R8-N1)** and the **dark-tuned `COMPARE_PALETTE` (R8-N2)** are visible regressions on the destination surface. |

---

## Roll-up: the 4 worth filing

1. **R8-N1 + R8-N2** (builder pair, P2) — one-issue twofer: kill the duplicate page header and paper-tune `COMPARE_PALETTE`. Smallest fix with the biggest payoff because the builder is *the* workflow destination.
2. **R2-M11 + R2-M14 + R2-F1** (contact-sheet triple, P2) — strip permanent thumb chrome, fix the per-row query fan-out, and tab-order falls out for free. Restores "this is a contact sheet" + "I can see the whole corpus" on the entry surface.
3. **R2-T9** (P2) — set the loupe heading in Newsreader. One-line fix the round-1 R2 issue missed.
4. **R4-N1** (P2) — drop the duplicate "Check every phase that is present" copy in the focus rail. Trivial; high readability win.

The ~14 P3 items can ride a single polish issue (call it **R11.polish**) or live as a TODO checklist in the canonical findings doc; they aren't workflow-breaking and most are one-liners.
