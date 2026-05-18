---
name: HimalayaUI
description: A dark, instrument-like thinking tool for indexing SAXS diffraction patterns.
colors:
  bg: "oklch(0.095 0.004 250)"
  bg-elevated: "oklch(0.145 0.006 250)"
  bg-subtle: "oklch(0.168 0.007 250)"
  bg-hover: "oklch(0.185 0.007 250)"
  fg: "oklch(0.95 0.005 250)"
  fg-muted: "oklch(0.64 0.009 250)"
  fg-dim: "oklch(0.46 0.009 250)"
  border: "oklch(0.27 0.008 250)"
  border-soft: "oklch(0.21 0.007 250)"
  accent: "oklch(0.82 0.08 220)"
  success: "oklch(0.72 0.08 155)"
  warning: "oklch(0.76 0.11 75)"
  error: "oklch(0.70 0.14 28)"
  peak-manual: "oklch(0.78 0.22 340)"
  phase-pn3m: "oklch(0.80 0.13 62)"
  phase-im3m: "oklch(0.78 0.12 160)"
  phase-ia3d: "oklch(0.76 0.12 300)"
  phase-fm3m: "oklch(0.78 0.13 18)"
  phase-fd3m: "oklch(0.76 0.12 318)"
  phase-hexagonal: "oklch(0.79 0.12 185)"
  phase-lamellar: "oklch(0.80 0.10 248)"
  phase-square: "oklch(0.79 0.12 132)"
typography:
  headline:
    fontFamily: "Plus Jakarta Sans, ui-sans-serif, system-ui, sans-serif"
    fontSize: "19px"
    fontWeight: 600
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
    letterSpacing: "0.04em"
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
  sm: "6px"
  md: "12px"
  full: "9999px"
spacing:
  sm: "8px"
  md: "16px"
  card-header: "56px"
components:
  button-primary:
    backgroundColor: "{colors.accent}"
    textColor: "#ffffff"
    rounded: "{rounded.sm}"
    padding: "4px 10px"
    typography: "{typography.body}"
  button-ghost:
    backgroundColor: "transparent"
    textColor: "{colors.fg-muted}"
    rounded: "{rounded.sm}"
    padding: "4px 10px"
    typography: "{typography.body}"
  button-danger:
    backgroundColor: "transparent"
    textColor: "{colors.fg-muted}"
    rounded: "{rounded.sm}"
    padding: "4px 10px"
    typography: "{typography.body}"
  input:
    backgroundColor: "{colors.bg}"
    textColor: "{colors.fg}"
    rounded: "{rounded.sm}"
    padding: "4px 8px"
    typography: "{typography.body}"
  card:
    backgroundColor: "{colors.bg-elevated}"
    textColor: "{colors.fg}"
    rounded: "{rounded.md}"
    padding: "16px"
---

# Design System: HimalayaUI

## 1. Overview

**Creative North Star: "The Darkroom"**

HimalayaUI works the way a darkroom works. The room is kept dark on purpose, not for mood, but so the thing that matters can be seen: the diffraction signal. The detector image, the 1D integration trace, and the Miller plot are luminous objects floating on a near-black field (`bg`, oklch 0.095). Everything else, the chrome, the cards, the labels, recedes into shadow so the eye goes straight to the data. Even the surface itself carries a faint fractal-noise grain at low opacity, which reads less like a texture and more like film.

The personality is a **confident expert** who happens to be quiet. The interface makes calls (the best phase, the best index, ranked, with a clear primary choice) and shows the reasoning behind them, but it never raises its voice to do so. Density is high and deliberate: a five-step type scale that bottoms out at 10.5px, a 250-hue cool-gray neutral ramp, 12px-radius cards that float on a deep ambient shadow. Color is rationed and semantic: a single ice-blue accent (hue 220) for interaction, a magenta (hue 340) reserved for manual peak edits, and an eight-phase data palette tuned so no two phases, and no phase and no accent, can be confused.

This system explicitly rejects its four neighbors. It is not **legacy scientific software** (no toolbar soup, no grey system chrome). It is not a **generic SaaS dashboard** (no hero-metric tiles, no identical card grids, no decorative gradients). It is not a **consumer app** (nothing is rounded for charm, nothing is bright for delight). And it is not a **bare academic utility** (the craft is visible, because visible craft is how a scientist decides the analysis can be trusted).

**Key Characteristics:**
- Near-black cool-gray field (hue 250) so luminous data leads; chrome recedes.
- One ice-blue accent, rationed; magenta reserved for manual edits; an 8-hue semantic phase palette.
- Dense, small type: a five-step scale from 10.5px to 19px, never extended ad hoc.
- Cards float on deep ambient shadow; everything else is flat.
- Fast, near-invisible motion: 120ms ease-out, no choreography, no bounce.

## 2. Colors

A near-monochrome cool-gray system (every neutral tinted to hue 250) carrying one interaction accent and a tightly-engineered semantic data palette. Color is never decorative; a color is a label.

### Primary
- **Ice Blue** (`accent`, oklch 0.82 0.08 220): the single interaction accent. Buttons, focus rings, links, the auto-detected peak markers in the trace overlay. Low chroma (0.08) keeps it calm against the dark field. It is the only chrome color with a hue, and it earns that by appearing on a small fraction of any screen.

### Secondary
- **Manual-Peak Magenta** (`peak-manual`, oklch 0.78 0.22 340): high-chroma, near-neon, used for one thing only: peak markers the user added by hand. Its hue sits a deliberate 120 degrees off the ice-blue accent and 22 degrees off the nearest phase color, so a manual peak can never be visually confused with an auto peak or an index.

### Tertiary
The **Phase Palette**: eight semantic data colors, one per known liquid-crystalline phase, used to color indices, trace overlays, and Compare-page series. All sit at luminance 0.76 to 0.80 and chroma 0.10 to 0.13 so they read clearly on the dark field without any one shouting. The hue spacing is engineered: each phase hue avoids a 20-degree exclusion zone around the accent (220) and the manual-peak magenta (340), and clears the high-chroma warning band near hue 75.
- **Pn3m Amber** (oklch 0.80 0.13 62) · **Im3m Sage** (oklch 0.78 0.12 160) · **Ia3d Violet** (oklch 0.76 0.12 300) · **Fm3m Coral** (oklch 0.78 0.13 18) · **Fd3m Rose-Purple** (oklch 0.76 0.12 318) · **Hexagonal Seafoam** (oklch 0.79 0.12 185) · **Lamellar Periwinkle** (oklch 0.80 0.10 248) · **Square Chartreuse** (oklch 0.79 0.12 132).

### Neutral
- **Darkroom Black** (`bg`, oklch 0.095 0.004 250): the base field. Near-black, faintly cool. Light theme inverts to oklch 0.985.
- **Lifted Surface** (`bg-elevated`, oklch 0.145 0.006 250): card backgrounds, one step up from the field.
- **Subtle / Hover Fills** (`bg-subtle` 0.168, `bg-hover` 0.185): inset wells and hover states.
- **Foreground** (`fg` oklch 0.95): primary text. **Muted** (`fg-muted` 0.64): labels, secondary text. **Dim** (`fg-dim` 0.46): the faintest tier, placeholders and disabled.
- **Borders** (`border` oklch 0.27, `border-soft` 0.21): hairline separators; soft is the near-invisible default for card edges.

### Status
- **Success Green** (oklch 0.72 0.08 155), **Warning Amber** (oklch 0.76 0.11 75), **Error Red** (oklch 0.70 0.14 28). Status-only; never used as decoration or accent.

### Named Rules
**The Semantic Color Rule.** Every color is a label. If a color does not encode meaning (interaction, provenance, phase, status), it does not belong on the screen. Decorative color is a bug.

**The Second-Channel Rule.** No meaning is ever carried by hue alone. Every color-coded element pairs its hue with a second channel: a shape, a text label, a position, or a pattern. An index distinguished only by color is broken, because it disappears under deuteranopia.

**The Rationed Accent Rule.** The ice-blue accent appears on a small fraction of any screen. Its calm and its rarity are what make it read as "interactive" rather than "decorated."

## 3. Typography

**Display / Body Font:** Plus Jakarta Sans (weights 400, 500, 600), with `ui-sans-serif, system-ui, sans-serif` fallback.
**Mono Font:** the platform monospace stack (`ui-monospace, SFMono-Regular, Menlo, monospace`).

**Character:** One humanist sans does all the structural work; a monospace face is reserved exclusively for numbers and scientific values. There is no display face and no serif. The pairing is sans-for-prose, mono-for-data, and the contrast between them is itself a signal: if it is monospaced, it is a measurement.

### Hierarchy
A deliberately small, dense five-step scale. The scale is fixed; sizes are never inlined ad hoc.
- **Headline** (600, 19px, line-height 1.2): page-level and rare. The largest type in the app.
- **Title** (600, 15px, line-height 1.3): card titles, primary input default state.
- **Body** (400, 13px, line-height 1.5): chat messages, hints, paragraph text. The readable default.
- **Label** (500, 11.5px, line-height 1.4, letter-spacing 0.04em, UPPERCASE): form and section labels.
- **Caption** (400, 10.5px, line-height 1.4): timestamps, microcopy, status badges, empty-state placeholders.
- **Data** (mono, 500, 11.5px; strong variant 600): q-values, scores, lattice parameters, phase chips. The numeric voice.

The seven semantic roles in CSS (`text-caption`, `text-label`, `text-meta`, `text-body`, `text-title`, `text-data`, `text-data-strong`) are the consumable surface; components express the role, not a size/weight/color triple.

### Named Rules
**The Fixed-Scale Rule.** Five steps, full stop. If a size seems missing, extend the scale in `styles.css` as a deliberate, reviewed change. Never scatter a one-off `text-[14px]`.

**The Monospace-Means-Measurement Rule.** The mono face is for values produced by the instrument: q, score, lattice constants, Miller indices. Never use it for prose, and never set prose data in the sans face.

## 4. Elevation

A hybrid system, but mostly flat. The dark field and its chrome are flat and depth-free; tonal layering (`bg` to `bg-elevated` to `bg-subtle` to `bg-hover`) handles most separation. Genuine shadow is spent on exactly one thing: the floating cards. A card carries a deep ambient drop shadow plus a 1px inset highlight rim, so it reads as physically lifted off the darkroom field, the way a print lifts off the developing tray. The light theme keeps the same shape at lower intensity.

### Shadow Vocabulary
- **Card Lift** (`box-shadow: 0 1px 0 rgba(255,255,255,.03) inset, 0 12px 36px -10px rgba(0,0,0,.7)`): the only structural shadow. Dark theme. The inset top highlight catches a faint edge of light; the large soft drop separates the card from the field.
- **Card Lift, Light** (`box-shadow: 0 1px 0 rgba(255,255,255,.6) inset, 0 6px 24px -8px rgba(20,24,40,.15)`): the light-theme equivalent, softer and shallower.

### Named Rules
**The Flat-Except-Cards Rule.** Surfaces are flat. Depth is conveyed by tonal steps in the neutral ramp, not by shadow. The single exception is the card, which floats. If something other than a card has a drop shadow, it is wrong.

## 5. Components

### Buttons
- **Shape:** gently rounded (6px radius, `rounded.sm`). Compact padding (4px vertical, 10px horizontal); buttons are small and dense.
- **Primary:** ice-blue fill (`accent`) with a matching border and white text; hover lifts brightness 10%. The assertive, reserved-for-the-main-action variant.
- **Ghost (default):** no fill, no border, muted text (`fg-muted`); hover brings text to full `fg` over a `bg-hover` fill. Most buttons in the app are ghost; the chrome stays quiet.
- **Danger:** like ghost, but hover shifts text to `error` red. No red fill at rest; destructive intent is shown only on approach.
- **Focus:** `focus-visible` outline in `accent`, 2px for primary, 1px for ghost/danger, always offset.

### Cards / Containers
- **Corner Style:** 12px radius (`rounded.md`), the softest corner in the system.
- **Background:** `bg-elevated`, one tonal step above the field.
- **Border:** a 0.5px hairline in `border-soft`, nearly invisible; the shadow does the separating, not the line.
- **Shadow Strategy:** the Card Lift shadow (see Elevation). Cards float; nothing else does.
- **Card Header:** an optional fixed 56px header strip with a 1px `border` bottom rule and 1rem horizontal padding, so the top edges of every headed card align across the workspace.
- **Internal Padding:** 16px (`spacing.md`) as the comfortable default.

### Inputs / Fields
- **Style:** recessed. `bg` background (darker than the card it sits in, reading as an inset well), 1px `border`, 6px radius, compact 4px/8px padding.
- **Focus:** border shifts to `accent` and a 1px `accent` outline appears. A quiet, precise focus, no glow.

### Navigation
- **App Header:** a thin 44px row, contents pushed to the right (a utility cluster); the left is intentionally empty so the page title can claim the plot card's top strip.
- **Page Nav:** a separate "TabRocker" control on its own row below the header. Keeping page-nav out of the global header is deliberate; it stops the control from crowding the page title.

### Signature Components
- **Score Bar:** a 2px-tall capsule track (`bg-hover`) with a fill colored by the relevant phase color, width set to the score percentage. A minimal, horizontal evidence bar; it shows confidence without a number shouting.
- **Semantic Dot:** an 8px colored disc carrying an `aria-label` and `role="img"`. The canonical second channel: it pairs a phase or status color with an accessible text label so the meaning survives without color.
- **Trace Overlay:** SVG peak and tick marks drawn over the diffraction trace. Auto peaks render in ice-blue, manual peaks in magenta; new nodes fade in over 140ms. Shape distinguishes mark kinds, so provenance never rests on color alone.
- **Phase Chip:** a monospace strong-data badge naming a phase; the numeric, instrument voice for a categorical result.

### Motion
Motion is fast and nearly invisible: a global 120ms ease-out transition on color, background, border, and opacity. Two named entry animations exist (`pal-in`, a 120ms fade for overlays and backdrops; `pal-scale`, a 140ms translate-and-scale for popovers) plus a 140ms overlay fade-in for freshly-drawn trace marks. No bounce, no elastic, no scroll choreography. Motion is polish, never spectacle.

## 6. Do's and Don'ts

### Do:
- **Do** keep the dark field as the default and let the diffraction data (trace, detector image, Miller plot) be the brightest thing on screen.
- **Do** treat every color as a label. Use the `accent` for interaction, `peak-manual` magenta only for hand-added peaks, the phase palette only for phase encoding, the status colors only for status.
- **Do** pair every color encoding with a second channel: a shape, a label, a position. Use the Semantic Dot pattern (color disc + `aria-label`) whenever a color carries meaning.
- **Do** stay inside the five-step type scale. If a size is missing, extend the scale in `styles.css` as a reviewed change.
- **Do** set scientific values (q, score, lattice, Miller indices) in the monospace face; set prose in Plus Jakarta Sans.
- **Do** make buttons ghost by default; reserve the ice-blue primary fill for the single main action on a surface.
- **Do** verify contrast at WCAG AA in both the dark and light themes.

### Don't:
- **Don't** look like **legacy scientific software**: no toolbar soup, no grey system chrome, no widget pile.
- **Don't** look like a **generic SaaS dashboard**: no hero-metric tiles, no identical card grids, no decorative gradient accents.
- **Don't** look like a **consumer or playful app**: nothing rounded for charm, no bright primary colors, no mascots or gamified flourishes.
- **Don't** look like a **bare academic utility**: an unstyled, function-only Jupyter-notebook-as-a-product. Visible craft is how the analysis earns trust.
- **Don't** carry meaning by hue alone; an index or peak distinguished only by color is broken under color blindness.
- **Don't** use `#000` or `#fff` as a surface or text color; every neutral is tinted to hue 250.
- **Don't** give a drop shadow to anything that is not a card. Depth elsewhere comes from the tonal neutral ramp.
- **Don't** add bounce, elastic, or scroll-driven motion. Transitions are 120ms ease-out and nearly invisible.
- **Don't** inline one-off type sizes (`text-[14px]`); extend the scale instead.
