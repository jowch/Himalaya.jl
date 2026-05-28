---
name: HimalayaUI
description: A light, warm-paper "publication plate" thinking tool for indexing SAXS diffraction patterns. Phase colour carries the surface; one terracotta accent is the interaction mark.
colors:
  paper: "oklch(0.978 0.006 85)"
  paper-sunk: "oklch(0.951 0.008 84)"
  plate: "oklch(0.992 0.004 90)"
  ink: "oklch(0.265 0.013 68)"
  ink-soft: "oklch(0.467 0.012 70)"
  ink-faint: "oklch(0.640 0.010 74)"
  hair: "oklch(0.882 0.008 80)"
  hair-strong: "oklch(0.806 0.010 78)"
  accent: "oklch(0.555 0.150 38)"
  frame-edge: "oklch(0.150 0.010 55)"
  success: "oklch(0.520 0.120 162)"
  warning: "oklch(0.620 0.130 70)"
  error: "oklch(0.520 0.170 28)"
  peak-manual: "oklch(0.550 0.200 340)"
  phase-pn3m: "oklch(0.570 0.150 58)"
  phase-im3m: "oklch(0.520 0.120 162)"
  phase-ia3d: "oklch(0.520 0.130 300)"
  phase-fm3m: "oklch(0.550 0.140 18)"
  phase-fd3m: "oklch(0.520 0.120 318)"
  phase-hexagonal: "oklch(0.540 0.130 350)"
  phase-lamellar: "oklch(0.505 0.150 264)"
  phase-square: "oklch(0.550 0.130 132)"
typography:
  display:
    fontFamily: "Newsreader, Georgia, ui-serif, serif"
    fontSize: "27px"
    fontWeight: 500
    lineHeight: 1.15
    letterSpacing: "normal"
  headline:
    fontFamily: "Newsreader, Georgia, ui-serif, serif"
    fontSize: "19px"
    fontWeight: 500
    lineHeight: 1.2
    letterSpacing: "normal"
  title:
    fontFamily: "Plus Jakarta Sans, ui-sans-serif, system-ui, sans-serif"
    fontSize: "15px"
    fontWeight: 600
    lineHeight: 1.3
    letterSpacing: "normal"
  body:
    fontFamily: "Plus Jakarta Sans, ui-sans-serif, system-ui, sans-serif"
    fontSize: "13px"
    fontWeight: 400
    lineHeight: 1.5
    letterSpacing: "normal"
  label:
    fontFamily: "Plus Jakarta Sans, ui-sans-serif, system-ui, sans-serif"
    fontSize: "11.5px"
    fontWeight: 500
    lineHeight: 1.4
    letterSpacing: "0.06em"
  caption:
    fontFamily: "Plus Jakarta Sans, ui-sans-serif, system-ui, sans-serif"
    fontSize: "10.5px"
    fontWeight: 400
    lineHeight: 1.4
    letterSpacing: "normal"
  data:
    fontFamily: "ui-monospace, SFMono-Regular, Menlo, monospace"
    fontSize: "11.5px"
    fontWeight: 500
    lineHeight: 1.4
    letterSpacing: "normal"
rounded:
  sm: "5px"
  md: "7px"
  full: "9999px"
spacing:
  sm: "8px"
  md: "16px"
  stage: "34px"
components:
  button-solid:
    backgroundColor: "{colors.ink}"
    textColor: "{colors.paper}"
    rounded: "{rounded.sm}"
    padding: "5px 12px"
    typography: "{typography.body}"
  button-ghost:
    backgroundColor: "transparent"
    textColor: "{colors.ink-soft}"
    rounded: "{rounded.sm}"
    padding: "5px 12px"
    typography: "{typography.body}"
  button-accent:
    backgroundColor: "{colors.accent}"
    textColor: "{colors.paper}"
    rounded: "{rounded.sm}"
    padding: "5px 12px"
    typography: "{typography.body}"
  input:
    backgroundColor: "{colors.plate}"
    textColor: "{colors.ink}"
    rounded: "{rounded.sm}"
    padding: "5px 9px"
    typography: "{typography.body}"
  plate:
    backgroundColor: "{colors.plate}"
    textColor: "{colors.ink}"
    rounded: "{rounded.md}"
    padding: "16px"
---

# Design System: HimalayaUI — "The Print"

> **Status:** authored 2026-05-27 as the canonical reference for the redesign's fidelity pass (milestone "HimalayaUI — The Print finish"). This describes the **target** Print system the remediation builds to; the code migration to it is tracked by R0a (#221) onward. Until R0a lands, the shipped tokens still default to the retired dark "Darkroom" values, this doc is the spec, not yet the shipped state. R9 (#232) verifies shipped == this doc.

## 1. Overview

**Creative North Star: "The Print"**

HimalayaUI is composed the way a publication plate is composed: on warm paper, in ink, with the figure as the subject. The app's analytical destination is a publication-quality figure (the stacked-offset waterfall), so the whole product is set on the medium that figure will be published on, light, warm paper (`paper`, oklch 0.978, hue 85), not a dark console. Detector images, which are intrinsically dark, sit as **dark windows in the paper** (`frame-edge`, hue 55); everything else is ink on paper.

Where the earlier dark "Darkroom" direction *rationed* a single accent and kept colour scarce, "The Print" lets **phase colour carry the surface**: colour here is the science (the phase identity), used generously and semantically. The chrome stays quiet, ink and hairlines on paper, and a single warm **terracotta** accent (`accent`, oklch 0.555 0.150 38) is the one interaction mark, the grease pencil: the reject ✕, the primary action, the live edit, the q-link highlight.

The personality is unchanged: a **confident expert** who is quiet. The interface makes calls and shows its reasoning; it never raises its voice. Density stays high. The serif display face (**Newsreader**) sets titles, the figure-plate signature; Plus Jakarta Sans does the structural work; monospace is reserved for measured values.

This system rejects its four neighbours: not **legacy scientific software** (no toolbar soup), not a **generic SaaS dashboard** (no hero-metric tiles, no identical card grids, no gradient accents), not a **consumer app** (nothing rounded for charm), and not a **bare academic utility** (visible craft is how the analysis earns trust).

**Key characteristics**
- Warm paper field (hue ~80) so the figure leads; detector images are dark windows set into it.
- Ink + hairline chrome; one terracotta accent as the interaction/grease-pencil mark.
- Phase colour is generous and semantic, paper-tuned so every phase reads at AA on `plate`.
- Newsreader serif for display titles; Plus Jakarta Sans for prose; monospace for measurements.
- The figure is framed as a **plate** (the one elevated object); most surfaces are flat paper + hairlines.

## 2. Colors

A warm paper-and-ink system (every neutral tinted to hue ~70–85) carrying one terracotta interaction accent and a paper-tuned semantic phase palette. Colour is never decorative; a colour is a label.

### Neutrals (paper + ink)
- **Paper** (`paper`, oklch 0.978 0.006 85): the base field. Warm near-white.
- **Paper-sunk** (`paper-sunk`, 0.951): recessed wells, hover fills, the rail margin.
- **Plate** (`plate`, 0.992): the brightest surface, the figure plate and cards float here.
- **Ink** (`ink`, 0.265 0.013 68): primary text and solid-button fills. **Ink-soft** (`ink-soft`, 0.467): secondary text, ghost buttons. **Ink-faint** (`ink-faint`, 0.640): captions, placeholders, disabled.
- **Hair** (`hair`, 0.882) / **Hair-strong** (`hair-strong`, 0.806): hairline separators and card edges; the line does the separating, paired with a soft shadow only on the plate.
- **Frame-edge** (`frame-edge`, 0.150 0.010 55): the warm near-black backing for detector images, the dark window set into the paper.

> **CSS-variable mapping (post-R0a).** The semantic utility tokens are remapped to these values: `--color-bg`→`paper`, `--color-bg-subtle`/`-hover`→`paper-sunk`, `--color-bg-elevated`→`plate`, `--color-fg`→`ink`, `--color-fg-muted`→`ink-soft`, `--color-fg-dim`→`ink-faint`, `--color-border`→`hair-strong`, `--color-border-soft`→`hair`, `--color-accent`→`accent` (terracotta). The dark `@theme` defaults, the `:root.theme-light` override, `color-scheme:dark`, the grain overlay, and the `theme` toggle are removed.

### Accent (the grease pencil)
- **Terracotta** (`accent`, oklch 0.555 0.150 38): the single interaction accent, kept warm so it belongs to the paper. The reject ✕, the primary `button-accent` and ink-solid actions, focus rings, links, the q-link cross-highlight, the live-edit mark. It is the only chrome hue, and it earns that by appearing on a small fraction of any screen.

### Manual-peak magenta
- **Manual-Peak Magenta** (`peak-manual`, oklch 0.550 0.200 340): paper-tuned (darker than the dark-era 0.78) for AA on the plate. Used for one thing only: peak markers the user added by hand. Its hue sits well clear of the terracotta accent (38) and every phase hue.

### Phase palette (paper-tuned, AA on plate)
Eight semantic data colours, one per liquid-crystalline phase, used for index/peak colour, trace overlays, folio phase strips, and Compare/series. All sit at luminance ~0.50–0.58 and chroma ~0.12–0.15, darker and a touch more chromatic than the dark-era values so each reads at WCAG AA on `plate`. Hue spacing clears a ~20° exclusion zone around the terracotta accent (38).
- **Pn3m Amber** (0.570 0.150 58) · **Im3m Sage** (0.520 0.120 162) · **Ia3d Violet** (0.520 0.130 300) · **Fm3m Coral** (0.550 0.140 18) · **Fd3m Rose-Purple** (0.520 0.120 318) · **Hexagonal Rose** (0.540 0.130 350) · **Lamellar Periwinkle-Indigo** (0.505 0.150 264) · **Square Chartreuse** (0.550 0.130 132).

> Im3m / Hexagonal / Lamellar are taken verbatim from the mockups (`:root`, byte-identical). Pn3m was also mockup-verbatim at L 0.585, but AA verification in R0b (#222) measured it at 4.29:1 on `plate` — below the 4.5:1 floor — so its luminance was nudged 0.585 → 0.570 (same chroma 0.150 + hue 58), which clears AA at 4.57:1 with no perceptible hue shift. The other four (Ia3d / Fm3m / Fd3m / Square) were **provisional** values in the same L/chroma band; R0b (#222) finalized them and confirmed all eight clear WCAG AA on `plate` (Square is the tightest of the rest at 4.54:1). See `phases.ts` and `test/phases.test.ts` for the pinned values + per-phase AA assertion.

### Status
- **Success** (oklch 0.520 0.120 162, sage), **Warning** (0.620 0.130 70, amber), **Error** (0.520 0.170 28, red). Paper-tuned for AA. Status-only; never decoration. (Note: the contact-sheet "kept" verdict dot is sage, not a generic green, T-4.)

### Named rules
**The Semantic Colour Rule.** Every colour is a label (interaction, provenance, phase, status). Decorative colour is a bug.
**The Second-Channel Rule.** No meaning rests on hue alone. Every colour-coded element pairs its hue with a second channel, shape, label, position, or pattern, so it survives deuteranopia. Use the Semantic Dot pattern (colour disc + `aria-label`).
**The Phase-Carries-the-Surface Rule.** Unlike the dark era's rationed accent, phase colour is used generously here, it is the science. The terracotta accent remains rationed (a small fraction of any screen); its rarity is what reads as "interactive".

## 3. Typography

**Display font:** **Newsreader** (serif, weight 500), with `Georgia, ui-serif, serif` fallback, the figure-plate signature, used for hero/figure/sample/card titles.
**Body font:** Plus Jakarta Sans (400/500/600), with `ui-sans-serif, system-ui` fallback, all structural prose and chrome.
**Mono font:** the platform monospace stack, reserved exclusively for measured values.

**Character.** Three voices, each a signal: serif = this is a title on a plate; sans = this is chrome/prose; monospace = this is a measurement (q, score, lattice, Miller indices, series ratio).

### Hierarchy
- **Display** (Newsreader 500, ~27px): the figure/sample/series hero title. The plate's name.
- **Headline** (Newsreader 500, 19px): section-level serif titles (folio, scoping).
- **Title** (Plus Jakarta Sans 600, 15px): card titles, primary inputs.
- **Body** (400, 13px): the readable default.
- **Label** (500, 11.5px, letter-spacing 0.06em, UPPERCASE): form/section labels and terracotta kickers.
- **Caption** (400, 10.5px): timestamps, microcopy, status badges.
- **Data** (mono 500, 11.5px): q-values, scores, lattice parameters, Miller indices, phase ratios.

### Named rules
**The Serif-Means-Title Rule.** Newsreader is for titles on a plate (sample name, figure title, card title, the progress numeral). Never set prose or chrome in serif.
**The Monospace-Means-Measurement Rule.** Mono is for values produced by the instrument. Never for prose.
**The Fixed-Scale Rule.** Extend the scale in `styles.css` as a reviewed change; never inline a one-off `text-[Npx]` (e.g. the off-scale `text-2xl` at `LoupePage.tsx:175`, L-2).

## 4. Elevation

Mostly flat: paper and its chrome are flat; tonal steps (`paper` → `paper-sunk` → `plate`) and hairlines handle most separation. Genuine shadow is spent on **the plate**, the figure composed as a print lifts off the paper the way a print lifts off the developing tray.

### Shadow vocabulary
- **Plate Lift** (e.g. `box-shadow: 0 1px 0 rgba(255,255,255,.6) inset, 0 6px 24px -8px rgba(40,30,20,.14)`): soft, warm, shallow. The figure plate, folio cards, the scoping worksheet, the loupe frame surround.
- Detector images are the dark exception: a dark window (`frame-edge`) set into the paper, not a lifted plate.

### Named rule
**The Flat-Except-the-Plate Rule.** Surfaces are flat; depth is tonal (the paper ramp) + hairlines. The plate (and plate-like cards) is the single lifted object. A drop shadow on anything else is wrong.

## 5. Components

### Buttons
- **Shape:** 5px radius (`rounded.sm`), compact padding.
- **Solid (`button-solid`):** ink fill, paper text, the primary/confirm action and active segmented-control segments (`.btn-solid {background:var(--ink); color:var(--paper)}`). This is the mockup's primary-action treatment, not the terracotta fill.
- **Accent (`button-accent`):** terracotta fill, paper text, reserved for the single grease-pencil action where the brand mark is the point (e.g. reject, "+ New series").
- **Ghost (default):** no fill, ink-soft text; hover brings text to `ink` over a `paper-sunk` fill. Most chrome buttons are ghost.
- **Focus:** `focus-visible` outline in `accent`, offset.

### Segmented controls / toggles
- Rest: ink-soft text on transparent. **Active segment: ink fill, paper text** (`.btn.on`/`.seg button.on`). Do **not** use the neutral `bg-bg-subtle` (a dark-era token) for the active state, that was the L-5/B-A defect.

### Cards / plate
- 7px radius (`rounded.md`), `plate` background, a `hair-strong` hairline, and the Plate Lift shadow. The figure plate widens to a max-width and centres; the rail recesses into `paper-sunk`.

### Inputs
- Recessed: `plate`/`paper-sunk` well, `hair` border, 5px radius. Focus shifts the border to `accent` with a 1px accent ring. Scoping values are "confident ink, re-openable" rather than a permanent open field (S-E).

### Detector frame
- A dark window: `frame-edge` background, set into the paper. Rings render in phase-neutral gold on triage surfaces (contact sheet, thumbnails) where the question is image quality, and in phase colour on the focus big detector where the question is identity. Rejected frames carry the hand-skewed terracotta grease-pencil ✕ (M-10).

### Signature components
- **Score bar:** a 2px capsule (`paper-sunk` track) filled in the relevant phase colour to the score percentage.
- **Semantic dot:** an 8px phase/status disc with `role="img"` + `aria-label`, the canonical second channel.
- **Phase chip:** a monospace data badge naming a phase, the instrument voice for a categorical result; on the contact sheet it replaces the static "Not indexed" string (M-6).
- **Mini-waterfall:** a read-only per-series SVG figure on folio cards (F-A), the folio is a wall of distinct figures, not a card grid.
- **Phase strip:** one cell per sample, coloured by its phase call (coexistence = two-phase gradient), with a transition caption.

### Motion
- Fast and nearly invisible: 120ms ease-out on colour/background/border/opacity. No bounce, no elastic, no scroll choreography. (The fractal-noise grain of the dark era is removed.)

## 6. Do's and Don'ts

### Do
- Keep the warm paper field; let the figure (trace plate, waterfall, detector window) be the subject.
- Treat every colour as a label; let phase colour be generous (it is the science) while the terracotta accent stays rationed.
- Set titles in Newsreader, prose in Plus Jakarta Sans, measurements in monospace.
- Make the active segmented-control segment ink-on-paper; make the primary action ink-solid, the grease-pencil action terracotta.
- Pair every colour encoding with a second channel; verify phase colours at AA on `plate`.

### Don't
- Don't use the dark-era neutral tokens (`bg-bg`, `bg-bg-elevated`, `text-fg`, `border-border`) or the ice-blue accent (hue 220) on a Print surface, they are retired (the whole dark↔Print seam, §0).
- Don't look like legacy scientific software, a SaaS dashboard, a consumer app, or a bare notebook.
- Don't carry meaning by hue alone.
- Don't give a drop shadow to anything but the plate (and plate-like cards).
- Don't inline one-off type sizes; extend the scale.
- Don't add bounce, elastic, scroll-driven motion, or the dark-era grain overlay.
