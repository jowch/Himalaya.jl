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
  frame-tag: "oklch(0.82 0.01 80)"
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
    fontSize: "31px"
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
  md: "5px"
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

> **Status: this is the living visual spec for "The Print", not a finished remediation log.** It is the canonical, extendable reference: every new screen, primitive, and token is expected to be authored against the rules below and to extend the named scales rather than scatter one-offs. The token primitives in the frontmatter are normative; the prose contextualizes how to apply them.
>
> History, for accuracy: the system was first authored 2026-05-27 for the redesign's fidelity pass; the token migration landed in R0a (#221) + R0b (#222) + R0c (#223), the shipped tokens (`styles.css` `@theme`, `phases.ts`) were verified to match this doc at R9 (#232), and the retired dark "Darkroom" defaults are gone. The system is a single light "Print" identity; there is no dark theme. The primitives were later relocated to `src/print/ui/` (the old `src/components/ui/` directory no longer exists) as the whole frontend moved under `src/print/`. When you add a screen, component, state, or token, add it here too, keeping the doc and the shipped system in lockstep.
>
> **Update (component-library extraction).** The Print's recurring patterns were extracted into closed-look primitives under `src/print/ui/` (Button, Card, SegmentedControl, PhaseChip, PhaseStrip, ModalShell, Kicker, IconButton, ScoreBar, Dot, Toast, HintText, and the base-layer Input, Field, Menu, Chip, Tooltip, and more) and the system is now **mechanically enforced**: `scripts/check-design.mjs` runs as a pure-absolute `lint:design` build step that fails the build on any inline appearance utility (`text-[…]`, `rounded-[…]`, raw colour literals, side-stripes) outside the appearance-authoring dirs (`src/print/ui/**`, plus the renderer dirs `src/print/plot/**`, `src/print/detector/**`, `src/print/comb/**`, `src/print/export/**`). Consumers pass placement-only `className`; appearance lives in the primitives. Two token truths changed at extraction: `rounded.md` 7px → 5px (radius collapsed to one 5px step), and `--color-print-accent` now sources from `--color-accent`. See `docs/superpowers/plans/2026-05-29-component-library-extraction.md`.

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

### Spacing & density

This is a dense instrument, not a marketing page; whitespace is a tool for grouping, not for air. The named spacing scale is intentionally thin:

- **sm** (8px): the gap inside a control row, between an icon and its label, between a chip and its neighbour.
- **md** (16px): the default card body padding and the gap between sibling controls. The `Card` `padding="md"` step is `p-4` (16px).
- **stage** (34px): the breathing gap that sets a figure plate or a page section apart from its surroundings, the one place the layout is allowed to open up.

Around those three named steps, surfaces use Tailwind's default spacing scale for the in-between values (card padding ranges 13–28px across the `Card` `padding` prop's `sm`/`md`/`lg` steps; page frames pad with `px-10 py-8`). The shared `.card-header` strip is a fixed 56px so card top edges line up across a workspace.

**The Thin-Scale Rule.** The named scale (`sm`/`md`/`stage`) is deliberately small and is an open extension point. When a recurring spacing value earns a name, add it to the scale as a reviewed change rather than scattering one-off paddings; until then, prefer the nearest Tailwind step over inventing a magic number.

### Copy and UX writing

The voice matches the brand: an **assured, precise, unshowy** colleague. Copy is direct and specific, never tentative, never cute, never a wizard's hand-holding. It states what is true and what the next action is.

- **The No-Em-Dash Rule.** Never use an em dash (or a "--" stand-in) in interface copy or in this system's docs. Break the thought into two sentences, or use a comma or colon.
- **Buttons** name the action as a verb in the user's terms: "Reanalyze", "New series", "Save figure", "Reject". Not "OK" / "Submit" / "Go".
- **Empty states** tell the user why the surface is empty and what to do next, in the `EmptyState` title + body pair (e.g. title "No series match", body "Clear the search or filter to see the whole folio.").
- **Error copy** is honest and specific, never a masquerade: a failed fetch says it failed and offers a recovery ("Couldn't load the folio" / "The series list failed to load. Try reloading the page."), and is never allowed to read as an empty "no results" state.
- **Loading copy** is quiet and momentary ("Loading series…") behind the skeleton, not a spinner narration.

**The Domain-Vocabulary-Is-Correct Rule.** The instrument's own terms are *precise*, not jargon to be softened: Pn3m, Im3m, Ia3d, Fm3m, Fd3m, Hexagonal, Lamellar, Square, Miller, q-value, hkl, reanalyze. These are the scientist's working vocabulary; "simplifying" them is a regression, not an accessibility win. Write to the expert, not around them.

## 2. Colors

A warm paper-and-ink system (every neutral tinted to hue ~70–85) carrying one terracotta interaction accent and a paper-tuned semantic phase palette. Colour is never decorative; a colour is a label.

### Neutrals (paper + ink)
- **Paper** (`paper`, oklch 0.978 0.006 85): the base field. Warm near-white.
- **Paper-sunk** (`paper-sunk`, 0.951): recessed wells, hover fills, the rail margin.
- **Plate** (`plate`, 0.992): the brightest surface, the figure plate and cards float here.
- **Ink** (`ink`, 0.265 0.013 68): primary text and solid-button fills. **Ink-soft** (`ink-soft`, 0.467): secondary text, ghost buttons. **Ink-faint** (`ink-faint`, 0.640): captions, placeholders, disabled.
- **Hair** (`hair`, 0.882) / **Hair-strong** (`hair-strong`, 0.806): hairline separators and card edges; the line does the separating, paired with a soft shadow only on the plate.
- **Frame-edge** (`frame-edge`, 0.150 0.010 55): the warm near-black backing for detector images, the dark window set into the paper.
- **Frame-tag** (`frame-tag`, 0.82 0.01 80): the light caption tint for the mono labels set *over* the dark `frame-edge` (the loupe frame-tag, R0c #223 / T-8); a dark-on-dark caption would vanish on the detector frame.

> **No legacy shim (post-R3-F, #259).** Earlier the dark-era neutral-ramp names (`--color-bg*`/`--color-fg*`/`--color-border*`) survived as a value-remapped shim; that shim is now excised. Use the canonical Print tokens directly — `paper`/`paper-sunk`/`plate`, `ink`/`ink-soft`/`ink-faint`, `hair`/`hair-strong`, `accent` (terracotta). Reintroducing an old name is self-revealing: Tailwind won't generate the utility and `var(--color-fg)` won't resolve. The dark `@theme` defaults, the `:root.theme-light` override, `color-scheme:dark`, the grain overlay, and the `theme` toggle were removed at R0a.

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
- **Display** (Newsreader 500, ~31px): the figure/sample/series hero title. The plate's name.
- **Headline** (Newsreader 500, 19px): section-level serif titles (folio, scoping).
- **Title** (Plus Jakarta Sans 600, 15px): card titles, primary inputs.
- **Body** (400, 13px): the readable default.
- **Label** (500, 11.5px, letter-spacing 0.06em, UPPERCASE): form/section labels and terracotta kickers.
- **Caption** (400, 10.5px): timestamps, microcopy, status badges.
- **Data** (mono 500, 11.5px): q-values, scores, lattice parameters, Miller indices, phase ratios.

### Named rules
**The Serif-Means-Title Rule.** Newsreader is for titles on a plate (sample name, figure title, card title, the progress numeral). Never set prose or chrome in serif.
**The Monospace-Means-Measurement Rule.** Mono is for values produced by the instrument. Never for prose.
**The Fixed-Scale Rule.** Extend the scale in `styles.css` as a reviewed change; never inline a one-off `text-[Npx]`. (Enforced: `scripts/check-design.mjs` fails the build on any inline `text-[…]` outside `src/print/ui/**`.) The display step is one such extension point: it sits at 31px (restored to the mockup hero size from the greenfield shell's conservative 27px — TY-DISPLAY), a deliberate, render-verified decision made by editing the `--text-display` token, never by inlining a size.

## 4. Elevation

Mostly flat: paper and its chrome are flat; tonal steps (`paper` → `paper-sunk` → `plate`) and hairlines handle most separation. Genuine shadow is spent on **the plate**, the figure composed as a print lifts off the paper the way a print lifts off the developing tray.

### Shadow vocabulary
- **Plate Lift** (e.g. `box-shadow: 0 1px 0 rgba(255,255,255,.6) inset, 0 6px 24px -8px rgba(40,30,20,.14)`): soft, warm, shallow. The figure plate, folio cards, the scoping worksheet, the loupe frame surround.
- Detector images are the dark exception: a dark window (`frame-edge`) set into the paper, not a lifted plate.

### Named rule
**The Flat-Except-the-Plate Rule.** Surfaces are flat; depth is tonal (the paper ramp) + hairlines. The plate (and plate-like cards) is the single lifted object. A drop shadow on anything else is wrong.

## 5. Components

For each component, lead with a short character line, then specify shape, colour assignment, states, and any distinctive behaviour.

### State taxonomy

Every interactive primitive draws from one shared vocabulary of states. The treatments below are what the shipped primitives in `src/print/ui/` actually implement; where a state is not yet systematized it is flagged as an extension point.

- **Rest.** The quiet default. Ghost buttons read ink-soft on transparent; inputs are a `plate` well with a `hair-strong` hairline; cards are flat `plate` + hairline.
- **Hover.** A small tonal lift, never a colour shift that implies meaning. Ghost-button hover brings text to `ink` over a `paper-sunk` fill; solid and accent buttons brighten (`hover:brightness-110`); outline buttons and Field rows fill to `paper-sunk`. All of it rides the global 120ms ease-out colour transition.
- **Active / pressed.** For toggles, the "on" state is the canonical **ink fill, paper text** (segmented-control active segment); an armed tool button switches to the terracotta accent fill with `aria-pressed`. Never a recessed `paper-sunk` fill for "on" (the L-5/B-A defect).
- **Focus-visible.** A terracotta `accent` ring, keyboard-only (`focus-visible`, offset). Buttons outline 2px accent; ghost and inverse variants use a 1px accent outline; inputs land the ring on the whole field via `focus-within:border-accent`; SVG hit targets (custom-index snap points, slider thumbs) stroke accent on `:focus-visible` only. Mouse clicks never show a ring. Bespoke `<button>`/`<a>` elements outside the primitives inherit the 2px accent ring by default from a base-layer rule in `styles.css`, so the state is systematized everywhere; primitive refinements and deliberate suppressions override it from higher layers (per-property — the 1px refiners redeclare `outline-offset-0` explicitly).
- **Disabled.** `disabled:opacity-45 disabled:cursor-not-allowed` on buttons; individual segmented-control segments can be disabled (e.g. Loupe with no sample selected). Disabled is a dimming, not a colour change.
- **Selected.** `Card selected` adds accent selection chrome: an accent-tinted border plus a subtle inset accent ring (flat, no shadow), and exposes `data-selected="true"` for tests and targeted CSS. This is the candidate-pick / row-selection signal.
- **Busy / loading.** Two distinct treatments, chosen by surface: a **skeleton** (boneyard `Skeleton`, the slow `.bone` pulse, fed a committed `*.bones.json` capture or an inline fixture) for the first paint of a data-driven region; a **progress fill** (`ProgressBar`, accent capsule, `role="progressbar"`) for a known-duration job. A bare spinner is reserved for indeterminate inline waits and is not yet a systematized primitive (extension point).
- **Error.** Inline and honest. Inputs take `invalid` → `border-error` plus `aria-invalid`, second-channeled by consumer error text. Page-level fetch failures render an `EmptyState` with failure copy (never a masquerading "no results"). `Toast` (`role="status"`) carries transient feedback with a leading status glyph + word (severity is the word + `aria-label`, not an edge hue).
- **Read-only.** A genuinely inert value renders as a static, non-interactive row, no chevron, no pointer cursor, not a `<button>` (the `Field` static branch). Controls must not lie: an affordance that does nothing is removed, not styled live.

**The Honest-Surface-State Rule.** Skeleton for the first load of a data region; spinner or progress for an in-flight action on already-loaded content; `EmptyState` (with a next-step CTA in the body) for a legitimately empty or failed surface; inline `border-error` + error text for a single invalid field. Never let a failure read as an empty state, and never show a control that cannot act.

### Buttons
- **Shape:** 5px radius (`rounded.sm`), compact padding.
- **Solid (`button-solid`):** ink fill, paper text, the primary/confirm action and active segmented-control segments (`.btn-solid {background:var(--ink); color:var(--paper)}`). This is the mockup's primary-action treatment, not the terracotta fill.
- **Accent (`button-accent`):** terracotta fill, paper text, reserved for the single grease-pencil action where the brand mark is the point (e.g. reject, "+ New series").
- **Ghost (default):** no fill, ink-soft text; hover brings text to `ink` over a `paper-sunk` fill. Most chrome buttons are ghost.
- **Focus:** `focus-visible` outline in `accent`, offset.

### Segmented controls / toggles
- Rest: ink-soft text on transparent. **Active segment: ink fill, paper text** (`.btn.on`/`.seg button.on`). Use `bg-paper-sunk` for any sunk surface, never a dark-era neutral, the active state must not read as a recessed fill (the L-5/B-A defect).

### Cards / plate
- 5px radius (`rounded.md`), `plate` background, a `hair-strong` hairline, and the Plate Lift shadow. The figure plate widens to a max-width and centres; the rail recesses into `paper-sunk`. (Radius collapsed to a single 5px step in the 2026-05-29 component-library extraction: `rounded.sm` and `rounded.md` are both 5px, so cards and controls share one corner; the two names survive only to allow a future re-differentiation by editing one token. The `Card` primitive is the single source of the plate look.)

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
- **Reflections card:** the lower-row table beside the detector window on the focus workspace. Column headers are `phase · order · q · d`. The mockup (`docs/redesign-mockups/focus-workspace.html`) labels the second column `hkl`; the shipped header is **`order`** by deliberate rename — SAXS peaks in our data model carry a 1-based ratio-series position, not an (h,k,l) Miller triple (`lib/seriesRatio.ts`, `src/phase.jl` `phaseratios`). "order" reads as "1st reflection, 2nd reflection…", which is what the value is. Treat the live `order` header as correct; the mockup `hkl` is the intentional divergence, not a regression.

### Motion

Motion is polish measured in milliseconds, never spectacle. It exists to soften state changes so the eye keeps its place, not to draw attention to itself.

- **Default transition.** A global rule transitions `color`, `background-color`, `border-color`, and `opacity` over **120ms ease-out** (`styles.css` `@layer base`). That is the whole motion budget for hovers, focus, selection, and tonal shifts.
- **Easing.** Ease-out only. No bounce, no elastic, no spring. Movement decelerates into rest and stops.
- **Functional animations, kept short.** A few one-shot keyframes carry genuine feedback: overlay fade-ins (`.anim-pal-in` 120ms, `.anim-pal-scale` 140ms popover entry, `overlay-fade-in` 140ms for newly-built trace marks) and the slow `.bone` skeleton pulse (1.8s). These convey state (something appeared, something is loading), so they earn their motion.
- **Motion conveys state only.** Nothing animates for delight. There is no page-load choreography, no scroll-driven sequence, no staggered entrance beyond the skeleton's brief reveal stagger. If a motion does not report a state change, it does not belong.
- **The Reduced-Motion Rule.** Under `prefers-reduced-motion: reduce`, every transition and animation is near-zeroed to 0.01ms (not removed, so `transitionend` / `animationend` listeners still fire) and `scroll-behavior` snaps to auto. The result has no perceptible spatial movement while functional feedback that depends on completing (spinners, progress, focus) is preserved.

> The fractal-noise grain, palette fades, and trace-overlay choreography of the retired dark "Darkroom" era are gone. Do not reintroduce them.

## 6. Do's and Don'ts

### Do
- Keep the warm paper field; let the figure (trace plate, waterfall, detector window) be the subject.
- Treat every colour as a label; let phase colour be generous (it is the science) while the terracotta accent stays rationed.
- Set titles in Newsreader, prose in Plus Jakarta Sans, measurements in monospace.
- Make the active segmented-control segment ink-on-paper; make the primary action ink-solid, the grease-pencil action terracotta.
- Pair every colour encoding with a second channel; verify phase colours at AA on `plate`.

### Don't
- Don't reach for the dark-era neutral names (`bg-bg`, `bg-bg-elevated`, `text-fg`, `border-border`) or the ice-blue accent (hue 220), the neutral shim was excised (R3-F, #259) so those utilities no longer exist; use `bg-paper`/`bg-plate`/`text-ink`/`border-hair-strong` (the whole dark↔Print seam, §0).
- Don't look like legacy scientific software, a SaaS dashboard, a consumer app, or a bare notebook.
- Don't carry meaning by hue alone.
- Don't give a drop shadow to anything but the plate (and plate-like cards).
- Don't inline one-off type sizes; extend the scale.
- Don't add bounce, elastic, scroll-driven motion, or the dark-era grain overlay.
