# "The Print" unification + fidelity remediation — plan

**Date:** 2026-05-27 · validates `main` @ `c5591c1` (redesign fully merged)
**Inputs:** the code-level gap report (`docs/2026-05-27-redesign-gap-report.md`), a live impeccable visual pass (app run against the dev-copy DB, Playwright over all 5 surfaces + focus workspace), and three parallel shipped-vs-mockup fidelity audits (one per workflow stage). Source of truth: `docs/redesign-mockups/*.html` + `docs/redesign-notes.md` (NOT the stale `DESIGN.md`).
**Decision driving this plan:** "The Print" (light warm paper, terracotta accent hue 38, Newsreader serif) is the app's single identity. **The dark "Darkroom" theme is retired.**

---

## 1. The one-paragraph diagnosis

The redesign shipped its *structure* faithfully — every surface's layout, routing, and the load-bearing event-sourced backend are correct, and the gap report verified that. What did **not** ship is the *finish*: the new surfaces wear "The Print" chrome (paper, terracotta kickers, the contact-sheet and trace-plate metaphors) but are built largely from **legacy dark-"Darkroom" components rendered unchanged on a light shell**. Because the app still defaults to `theme:"dark"` (`state.ts:284`), those components resolve their neutral tokens to near-black/near-white and their accent to ice-blue (hue 220) on warm paper. On top of that theme seam sit a set of **missing designed affordances** (the scoping proposal UX, the focus phase-call rail, the folio mini-figures, several contact-sheet and builder controls) that the acceptance audit couldn't see because each surface's *acceptance line* was met.

## 2. The keystone insight (why this is smaller than it looks)

The ~40 findings are dominated by **one root cause with one fix**. Most offending components use the *neutral ramp* tokens (`bg-bg`, `bg-bg-elevated`, `text-fg`, `text-fg-muted`, `border-border`) and the shared `--color-accent`. If we **redefine those token values to Print** at the `:root`/`@theme` level (paper/ink/plate/hair + terracotta accent) and delete the dark defaults + the `theme-light` override + the `theme` toggle, the majority of theme findings resolve **by inheritance**, with no per-component edits. What remains after that sweep is a small set of *semantically* wrong tokens (e.g. a kept-dot using `bg-success` green instead of sage, the phase palette tuned for dark) and the genuinely missing affordances.

So the plan front-loads one keystone workstream (**R0**), after which the per-surface issues are mostly additive feature work that can run in parallel.

## 3. Findings → workstreams (consolidated)

Severity: **P1** workflow-breaking / promised-and-missing · **P2** degraded but functional · **P3** polish.

### R0 — The Print token unification (KEYSTONE, P1, foundational)
Resolves by inheritance: T-1…T-8 (loupe), A-1…A-4 (focus inner), B-A/B-B/B-C (builder gutter + toggles), S-B/S-C (scoping buttons), L-5/L-7 (toggles), and most ice-blue-accent findings.
- **R0a — token remap + serif.** Rewrite `styles.css`: make the Print values (`paper/plate/paper-sunk/ink/ink-soft/ink-faint/hair/hair-strong` + terracotta accent) the *default* `@theme` neutral ramp; delete the dark `@theme` defaults, the `:root.theme-light` block, `color-scheme:dark`, and the fractal-grain overlay. Remove the `theme` field + toggle from `state.ts`/persist migration + any toggle UI. Import **Newsreader** and add a serif display type role.
- **R0b — phase palette retune.** Rework `phases.ts` to the mockup's paper-tuned phase hues (darker L≈0.50–0.58, higher chroma, AA on `--plate`); keep the colour-blind second-channel rule.
- **R0c — residual per-component token cleanup.** After R0a/R0b, sweep the components the audits named for *semantically* wrong tokens that don't auto-resolve: loupe kept-dot (`bg-success`→sage), `ThumbnailGallery` selection rings, `MultiTracePlot`/`MemberMetaGutter` tooltip+selection chrome, `FocusDetectorPanel` inner frame, dark-on-dark caption/badge contrast (T-8).

### R1 — Contact-sheet affordances (P1/P2)
M-1 progress block ("N/M screened" + terracotta bar), M-2 per-sample screened mark + unscreened row tint (coordinate **#162**), M-3 floating cross-row cull bar, M-6 phase-chip status + hollow dot, M-7 frame-number badges, M-10 grease-pencil ✕ SVG (sheet + loupe), L-3 beamtime H1 + descriptive sub, L-4 square thumbnails, L-5..L-11 layout/copy. Tag on-ramp (M-4/F-1) → **#207**.

### R2 — Loupe fidelity (P2/P3)
Mostly absorbed by R0c (tokens). Residual affordances: M-8 signal-strength meta meter, M-9 contact-sheet/loupe segmented view switch, L-9/L-10 subtitle "exposure N of M" + tag value model.

### R3 — Focus header (P1) — the L-6 fix
Give `PlotCard` a focus-variant header (a `headerSlot`/`variant` prop, or a dedicated focus header component) that drops the experiment-picker affordance and renders: terracotta "Integration" kicker → Newsreader serif **sample name** → mono `smp · beamtime · representative exposure` subline, seeded from the route's sample (not the unset global picker). (`PlotCard.tsx:386-416`, `FocusWorkspaceLayout.tsx:60`.)

### R4 — Focus rail as the output (P1/P2)
L-9/L-10/L-11: a distinct **"Phase call"** block above "Candidate indexings" (per-phase serif name, score bar, lattice, series ratio √2:√3:√4, coexistence header), candidates as explicit **multi-select checkboxes**, and a phase-colour candidate hover-preview channel distinct from the accent q-link (with losing-peak dim). Drop the legacy Miller inset from the rail (L-12). Reflections table + 3rd q-link target → **#209**.

### R5 — Focus affordances (P1/P2)
F-11 representative-exposure switcher (mini-detector thumbnails, set-rep), F-12 Notes topbar/drawer fallback below `xl` (currently *unreachable* below the breakpoint), F-13 per-sample stepper in the topbar, F-14 active Index/stage tab on `/sample/:id`.

### R6 — Series folio (P1/P2)
F-A/N-1 live per-series mini-waterfall SVG (reuse `MultiTracePlot`/`MemberTraceLayer` geometry, read-only), F-B per-sample phase strip + transition caption, F-C "+N new match" + draft pills, F-G "+ New series" primary action, F-D filter chips + Variable/Largest sorts, F-H provenance footer, X-1 serif title.

### R7 — Series scoping proposal UX (P1) — the biggest single gap (S-A / L-8)
Surface `proposeOrdering`'s output (the plumbing + build-gate already exist): autogroup summary card, parsed "Ordered by" variable field, per-row **trace sparklines**, amber "check the read" confidence flags with **re-openable ink values** (confirm-not-fill-out), "Himalaya also found" loose-match section, preview phase strip, narrative gate footer, worksheet-plate framing, reorder grips + in-session Undo.

### R8 — Series builder (P2)
B-F **offset slider + log/linear scale toggle** (core compose controls, currently absent), B-G floating offset dock in full-bleed, B-J figure-as-plate container, B-H kicker tag-row + auto figure caption, B-I autogroup read-rail card. Heatmap + peak-tracking → **#208**.

### R9 — Docs: regenerate DESIGN.md (P2) — M-1
Regenerate `DESIGN.md` to "The Print" from the unified token system (after R0) via `impeccable document`. Retire the Darkroom description.

### R10 — Backend: image route robustness (P2) — L-2
`GET /api/exposures/:id/image` should return a graceful 404 (or placeholder) when the source TIFF is missing, not an unhandled `ArgumentError`→500.

## 4. Reconciliation with the already-deferred batch
These are existing issues; fold them in, don't duplicate. **Land-order coupling** (same files/surfaces):
- **#207** corpus tag/rename UI → lands with **R1** (contact-sheet tag on-ramp + loupe tags).
- **#208** heatmap + cross-trace peak-tracking → lands with **R8** (builder).
- **#209** reflections table + 3rd q-link target → lands with **R4/R5** (focus).
- **#162** screened mark → part of **R1**.

## 5. Dependency DAG
```
                 ┌──────────────────────────────────────────────┐
 R0a ─→ R0b ─→ R0c (keystone: token unification + serif + palette)│
                 └──────┬───────────────────────────────────────┘
                        │ (R0 settles the theme; everything else is additive)
        ┌──────┬────────┼────────┬────────┬────────┬────────┐
        ▼      ▼        ▼        ▼        ▼        ▼        ▼
       R1     R3       R4       R5       R6       R7       R8
     (+#162) (focus  (focus   (focus  (folio) (scoping)(builder
     (+#207)  hdr)   rail)    affd.)            BIGGEST) +#208)
                      +#209 ──┘
   R9 (DESIGN.md, after R0)   R10 (backend, independent — can start anytime)
```
- **R0 is the gate.** Doing per-surface token fixes before R0 means fighting the dark defaults; doing them after means most are already resolved.
- **R0a → R0b → R0c** is a short serial chain (palette depends on token base; cleanup depends on both).
- **R1–R8 are parallel** after R0, one owner per surface. R4 and R5 both touch the focus rail/layout — give them a land-order or one owner.
- **R10** (backend 500→404) is fully independent; can land immediately.
- **R7 (scoping)** is the largest single body of new UI work; size it as its own multi-task effort with its own detailed plan.

## 6. Suggested sequencing
1. **R0a + R0b + R10** first (R10 parallel, trivial). R0a/R0b is the single highest-leverage change — re-screenshot every surface immediately after to confirm the inheritance sweep.
2. **R0c** cleanup pass (residual tokens the sweep missed).
3. Then parallel: **R3** (focus header, small + high-visibility), **R6** (folio miniatures), **R8** (builder controls), **R1** (contact-sheet affordances, coordinate #162/#207).
4. **R7** (scoping proposal UX) as a dedicated effort.
5. **R4/R5** (focus rail + affordances, coordinate #209), **R2** residuals.
6. **R9** (regenerate DESIGN.md) last, once the tokens are final.

## 7. Method notes
- Each R-issue should get a detailed per-issue TDD plan (writing-plans) when picked up, exactly as the original redesign did. This document is the *map*, not the per-issue plans.
- Verify each surface live against the screenshots captured in this session (`$CLAUDE_JOB_DIR/*.png`) after R0, and after each surface issue.
- The dev-env path bridge (`/data`→`/Volumes/data` in the copy DB) is needed to render images/traces locally; backup at `$CLAUDE_JOB_DIR/himalaya.db.prebak`.
