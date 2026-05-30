# Plotting mockups — impeccable audit (2026-05-30)

Technical audit of the two high-fidelity mockups (`2026-05-29-focus-plot.html`, `2026-05-29-series-plot.html`). These are prototypes — a11y and responsive are deliberately implementation-phase concerns — so this audit is primarily the **implementation checklist** (§7.4 of the design spec), not a defect list against the mockups.

## Anti-patterns verdict: PASS
Not AI-generated-looking. No AI palette (distinctive warm-paper/ink), no gradient text, no glassmorphism, no hero-metric tiles, no identical card grids, no bounce easing, no generic fonts. One systemic tell: **em dashes in prose** (45 across both files), which the impeccable copy law bans. Export modal is the only "modal-as-first-thought" candidate, and it's defensible (focused preview-then-copy).

## Health score: 13/20 — Acceptable (top edge, near Good)

| # | Dimension | Score | Key finding |
|---|-----------|-------|-------------|
| 1 | Accessibility | 2 | Candidate/member rows are clickable `<div>`s (no button role/tabindex/keyboard); q-link + migration-track are hover-only |
| 2 | Performance | 3 | Hover is class-toggle (no re-render), SVG-only, bounded; minor full re-render on candidate toggle |
| 3 | Theming | 3 | Consistent CSS vars; single light identity by design; deliberate raw hex only in the decoupled export |
| 4 | Responsive | 2 | Zero `@media` breakpoints; fixed-width rails; touch targets < 44px |
| 5 | Anti-patterns | 3 | Distinctive/intentional; em-dash copy is the lone systemic tell |

13 is expected for prototypes: the weak dimensions (a11y, responsive) are built in the implementation phase. The mature dimensions ride the already-hardened Print system.

## Findings → implementation checklist

**P1 (build phase, required)**
- **Clickable divs need button semantics.** `.cand`, `.member`, `.as-foot` fire on click but aren't `<button>`/keyboard-operable. In React, make them `<button>`s. (WCAG 2.1.1, 4.1.2)
- **q-link triple + Series migration track are hover-only.** No keyboard or touch path to the single best feature. Need focus-driven equivalents (focus a peak → same highlight) and a tap model. (WCAG 2.1.1)

**P2 (next pass)**
- **No responsive breakpoints.** Fixed `350px`/`320px` rails + fixed `det-size`; below ~1100px the lower row and rail crowd. Needs a collapse strategy.
- **Touch targets < 44px** (segmented controls ~26px, exposure thumbs 30px, anchor hits 18px). The app's `.icon-button` 44px pattern and `@media (pointer:coarse)` seg rule exist in `styles.css` — reuse them. (WCAG 2.5.5)
- **`text-ink-faint` contrast** (oklch L0.64 on L0.978 paper) likely below 4.5:1 for small mono axis ticks/captions — measure; promote to `ink-soft` if it fails. (WCAG 1.4.3)
- **No landmark structure** (only `<aside>`); topbar/work/stage are plain divs. Add `<header>`/`<main>`/`<nav>`.

**P3 (polish)**
- **Em dashes throughout prose** (45). Swap for commas/colons/parentheses (impeccable copy law).
- **Rail section labels are styled divs, not headings** — `<h2>/<h3>` for screen-reader navigation.
- **Full re-render on candidate toggle** — fine at this scale; in React scope updates to changed marks.

## Positive findings (keep)
- **Second-channel encoding is strong** — colour=phase always paired with silhouette / fill-outline / struck-ghost. Survives colour-blindness (the product's hard rule).
- **Hover is CSS-class-driven**, not re-render — the right performance pattern; carries straight to React.
- **The decoupled export's raw hex is correct-by-design**, not a theming violation.
- **Distinctive, intentional visual system** — passes the slop test cleanly.

## Recommended commands (implementation phase)
1. `impeccable harden` — button semantics, keyboard paths, landmarks, focus management.
2. `impeccable adapt` — responsive breakpoints + 44px touch targets.
3. verify `text-ink-faint` contrast; promote small text if it fails.
4. `impeccable clarify` — strip em dashes.
5. `impeccable polish` — final pass.
