# Impeccable primitive checklist (greenfield src/print/ui)

Every `src/print/ui/` primitive must satisfy this before commit. Distilled from the
`impeccable` skill (register: **product**) + `DESIGN.md` named rules + `PRODUCT.md`.
Source of truth for the per-task "Impeccable pass" step.

## A. Second-Channel rule (HARD — color-blind requirement)
- No meaning rests on hue alone. Any color-coded element pairs hue with a second
  channel: shape, label, position, or pattern. Phase/status discs use the Semantic
  Dot pattern (color + `role="img"` + `aria-label`). Verify the encoding still reads
  in grayscale.

## B. Semantic color + Phase-carries-surface
- Every color is a label (interaction / provenance / phase / status). No decorative color.
- Phase color via `phaseColor()` only (never a literal). Terracotta `accent` is the one
  rationed interaction mark — used for the single primary/grease-pencil action, focus
  rings, the live-edit/reject mark. Not for decoration or inactive states.

## C. Interaction states (product register — ship them all)
- Every interactive primitive defines: default, hover, focus-visible, active, and
  disabled where applicable. Focus is a visible `focus-visible` outline in `accent`,
  never removed without replacement. No heavy/full-saturation accent on inactive states.

## D. Touch target
- Interactive targets are ≥44×44px effective, OR the primitive documents (in a code
  comment) that it is a known-dense control whose hit area is enlarged by its container.
  (The compact toggles were a flagged audit defect — do not reproduce silently.)

## E. Voice discipline (DESIGN.md type rules)
- Serif (Newsreader, via `.text-display`/`.text-headline`/`.text-headline-lg`) = titles
  on a plate only. Sans = chrome/prose. Monospace = measured values only (q, score,
  counts, ids, lattice). No serif in labels/buttons/data; no off-scale `text-[Npx]`.

## F. Tokenized appearance (Fixed-Scale + radius)
- Sizes/colors/radii come from named scale roles / `@theme` tokens. Radius is
  `rounded-sm` (5px) / `rounded-full`. No arbitrary `rounded-[Npx]`, no raw color literal
  in a consumer string. (`lint:design` enforces outside `print/ui/`; we hold `print/ui/`
  to the same bar by choice.)

## G. Absolute bans (impeccable shared + product)
- No side-stripe accent border (`border-l/r` >1px as a colored accent) — full border +
  leading icon/word instead. No gradient text. No decorative glassmorphism/`backdrop-blur`.
  No display font in UI labels. No bounce/elastic motion; transitions are 120ms ease-out
  on color/background/border/opacity only — never on layout properties. Respect
  `prefers-reduced-motion`. No em dashes in any rendered copy.

## H. Flat-except-the-plate
- Only `Card` (the plate) and plate-like surfaces carry the Plate-Lift shadow. Other
  primitives are flat: tonal steps (paper → paper-sunk → plate) + hairlines do separation.

## I. Closed-look API hygiene
- Appearance internal; consumer `className` is placement-only and appended last. Variants
  are a `Record<Variant,string>` map + local `cx()`. Public union types exported from the
  barrel. `data-*` attributes set for every variant/state so tests avoid class strings.

## J. Deterministic scan (impeccable detect)
- `npx impeccable detect src/print/ui/<Name>.tsx` returns no anti-pattern hits. This is
  the deterministic floor (side-stripes, gradient text, glass, off-scale values, provider
  tells); a clean scan is necessary but NOT sufficient — A–I still apply by inspection.

## Per-primitive sign-off
For each primitive, confirm: ☐ A ☐ B ☐ C ☐ D ☐ E ☐ F ☐ G ☐ H ☐ I ☐ J, and that the
Storybook story renders every state listed in C. A primitive that legitimately has no
interactive state (e.g. a static badge) marks C/D N/A with a one-line reason.
