# Component-library extraction — design spec (2026-05-29)

> **Status:** design, ready for implementation-planning. Second project in the program
> **audit → component library → trace-plotting redesign**.
> **Inputs:** `docs/2026-05-29-frontend-audit-impeccable.md` (Theme A), the approved brainstorming
> design, and the exhaustive call-site inventory in
> `docs/2026-05-29-library-extraction-inventory.json` (8 concerns, ~300 call-site entries — overlapping
> where a site touches multiple concerns — used as the migration checklist rather than reproduced here).

## 1. Goal & principle

The frontend's design system is well-authored but **unenforced**: color stays consistent because it
has a single source (`phaseColor()`) and a regression test, while type, radius, and component look drift
because nothing stops a consumer from slapping appearance utilities on top. The audit proved the failure
mode empirically (a `ScoreBar` primitive shipped, nobody adopted it, three divergent inline bars exist
instead).

**Principle: the library owns *appearance*; consumers own *placement*. A guard makes it stick.**

This project does three things, in this order:
1. **Found** the system: tokenize the radius scale (currently dead), add the scrim token, de-dupe the
   accent, adopt-or-delete the existing primitives, and stand up the enforcement guard + a static catalog.
2. **Promote** the high-duplication inline patterns to closed-appearance primitives.
3. **Sweep** the codebase onto them and **ratchet the guard to fully enforcing.**

Out of scope (chosen during brainstorming): a publishable/decoupled package, Storybook, refactoring
feature-component state/query coupling, the visual peak/strip glyph encodings (→ plotting redesign), and
any plotting primitives (→ plotting redesign).

## 2. The styling contract — "closed look, open placement"

Every primitive follows the existing `Button.tsx` pattern (a `Record<Variant, string>` class map, base
classes first, `...props` spread, typed props extending the DOM element). Concretely:

- **Appearance is expressed only through semantic props** (`variant`, `size`, `tone`, and domain props
  like `phase`). The primitive computes its own classes/colors internally; `phaseColor()` and
  `color-mix` recipes live inside the primitive, never at a call site.
- **`className` is accepted but placement-only**: margin, width, `flex-1`/`shrink-0`/`self-*`, grid/flex
  position, transforms. It is appended last. Appearance utilities in a consumer `className`
  (`bg-*`, `text-*` color/size, `border-*` color/width, `rounded-*`, `shadow-*`, `backdrop-*`, raw color
  literals) are **banned by the guard** outside `src/components/ui/**`.
- **Primitives that currently leak appearance via props are fixed during extraction.** Canonical example:
  `ScoreBar(score, color)` → `ScoreBar({ value, phase, size })`, deriving color from `phaseColor(phase)`
  internally and replacing the hardcoded `h-0.5 w-full` with a `size` enum.
- **Tests assert rendered semantics, never class strings** — this is the *target* convention and the
  rule for all new primitive tests (the primitive exposes `data-variant`/`data-tone`/`data-size`/
  `data-state` and the test asserts *that*). It is **not yet universally true**: a handful of existing
  tests assert class strings or a leaked prop and **must be rewritten in the same PR that changes the
  underlying primitive** — `test/ui/ScoreBar.test.tsx` (asserts `style.background` + `score`/`color`
  props), `test/scoping.test.tsx:117-119` and `test/ScopingFoot.test.tsx:31-33` (assert
  `bg-ink`/`text-paper`/`.not.toContain("bg-accent")`). See §7. Existing `data-testid`s ride through
  `...props` and must be preserved exactly (they back E2E specs).
- **No new dependency.** A 3-line local `cx(...parts)` join helper is allowed; no cva, clsx, or
  tailwind-merge.

`src/components/ui/` is the library home; everything is exported from `ui/index.ts` (the barrel). The
guard's exclusion of `ui/**` is precisely what lets primitives be the one place appearance is authored.

## 3. Token foundation (Phase 0, `src/styles.css` `@theme`)

```css
/* --- Radius scale: make the radius tokens REAL --- */
/* Tailwind v4 generates rounded-* from the --radius-* NAMESPACE, not singular --radius.
   `--radius: 6px` is DEAD (no var(--radius) consumer) and rounded-sm/md currently resolve to
   Tailwind STOCK defaults (sm 4px, md 6px), NOT any design value. The 5 mockups scatter 4/5/6/7px
   across surfaces; we normalize the primary surfaces to 5px (the dominant mockup value, and the
   approved decision — see §8). DESIGN.md prose still says rounded.md=7px; it is now wrong and should
   be corrected to 5px alongside this change. */
--radius-sm:   5px;     /* buttons, inputs, chips, segmented controls */
--radius-md:   5px;     /* cards, plates, sheets, modals (matches mockups; was DESIGN.md 7px) */
--radius-full: 9999px;  /* pills, dots */
/* DELETE --radius: 6px */

--color-scrim: oklch(0.05 0 0 / 0.65);   /* the literal 3 modal backdrops hand-inline → bg-scrim */
--color-print-accent: var(--color-accent); /* was a duplicated terracotta literal → one source */
```

`--radius-sm` and `--radius-md` deliberately both resolve to 5px: this keeps the utility-name semantics
(controls use `rounded-sm`, cards use `rounded-md`) so a future contributor could re-differentiate by
editing one token, while every current surface lands on the single 5px corner the mockups favor.

Consequences (all deliberate, reviewed visual changes — see §8 sign-off):
- `rounded-md`/`rounded-sm` sites snap from stock Tailwind (6px/4px) to 5px.
- `.card` (`styles.css:175`) hardcodes `border-radius: 12px` and modal frames ship `rounded-xl` (12px);
  both fold onto `rounded-md` (5px). **Sign-off item** (touches the figure plate, the most prominent surface).
- `--unindexed` is **not** added: per the approved decision, unindexed cells/peaks/rings fold onto the
  existing `ink-faint` (0.640) token (the mockups' dedicated `--unindexed` 0.660 is imperceptibly different
  and `ink-faint` is DESIGN.md's "placeholder/disabled" role).
- **Bare `rounded` (185 sites) is out of scope**: v4's bare `rounded` is its own 4px default independent
  of `--radius-*`; sweeping it is a separate low-value mechanical pass. Documented, not done here.
- The `@theme` test must assert the **resolved utility output** (a `rounded-md` element computes 5px), not
  just the custom-property presence — a stray `@layer`/specificity issue could shadow the utility while the
  property exists.

## 4. The enforcement guard — `scripts/check-design.mjs`

A zero-dependency Node script that globs `packages/HimalayaUI/frontend/src/**/*.{ts,tsx}`, **excludes
`src/components/ui/**`**, and flags banned patterns. It is a deliberately *approximate* line/token
detector: sound for the patterns below against today's surface, with a documented escalation path to an
AST pass (Babel is already a transitive dep via `@vitejs/plugin-react`) if a future false-negative appears.

**Ban set:**
1. `\btext-\[` — inline arbitrary type size or color (`text-[10px]`, `text-[var(--color-…)]`). Fixed-Scale.
2. `\brounded-\[` — inline arbitrary radius (9 sites today: `rounded-[1px]`/`[1.5px]`/`[3px]`/`[7px]`/`[10px]`). Mirrors #1 so the new radius scale can't be bypassed the way `text-[` was. (Bare `rounded` is *not* banned — see §3.)
3. Raw appearance color in a color utility: `\b(bg|text|border|fill|stroke|ring|outline|from|via|to|decoration)-\[…(oklch|oklab|rgba?|#[0-9a-fA-F]{3,8}|hsla?)`. **`shadow-[…]` is excluded** — there are 4 real `shadow-[…rgba…]` Plate-Lift sites in tsx (the worksheet inline shadows) and they are legitimate per DESIGN.md §Elevation. Decide which rule owns `bg-[oklch(…)]` (this one) so it isn't double-counted into rule #5's baseline.
4. Side-stripe: `\bborder-(l|r)-(?:4|[2-9]|\[)` — `InfrastructureBanner.tsx:55` (Toast at `Toast.tsx:85` is inside `ui/`, excluded; both are converted when the side-stripe fix lands — now pulled forward, see §6 Phase 0).
5. Raw color literal **anywhere** in a line (covers `style={{…}}`, which routinely spans multiple lines so it can't be anchored on the `style={{` token): first **strip all `var(--color-*)` tokens from the line**, then flag a remaining raw `oklch(`/`oklab(`/`rgb[a](`/`#hex` literal. This passes every legit computed-color site (their non-`var` color is an interpolated `${…}` expression with no raw literal: `SampleStatusChip.tsx:25`, `MemberMetaRow.tsx:341`, and the bare-string-with-`var()` cases `FocusReflectionsTable.tsx:156`, `SeriesFolioCard.tsx:144`) and flags the raw `oklch(0.05…)` scrim literals. **No file-level allowlist is needed for this rule** — the token-stripping does the discrimination. Residual gap (documented, accepted): a color laundered through a `const x = "oklch(…)"` then `style={{background:x}}` is not caught; the catalog + review cover it.

**Allowlist (rule #3 only — arbitrary-value color *utilities* in code that legitimately needs them):**
`src/phases.ts`, `src/lib/comparison/coloring.ts`, `src/lib/figure-export/**`, `MemberHeatmapLayer.tsx`,
`DetectorImage.tsx`, `FocusDetectorPanel.tsx`, `src/main.tsx`. Rules #1/#2/#4/#5 take **no** file allowlist.

**Baseline schema:** a committed `scripts/design-baseline.json` keyed on **(rule, normalized-violation-text)
content hash**, *not* (file, line) and *not* (file, rule)-count. Rationale: content-keying makes moving a
violation between files (routine during a sweep, and during the repo's "one responsibility per file"
splits) a no-op for the baseline, so only *net-new* appearance is ever blocked; file-count keying would
either false-block a legitimate relocation into a new file or let a move-and-add escape via a flat global
sum. The seed baseline is generated from the inventory JSON's violation lists (the ~96 `text-[`, the 9
`rounded-[`, the raw-color/scrim sites, the one enforced side-stripe).

**Wiring:**
```jsonc
"lint:design": "node scripts/check-design.mjs",
"build": "npm run lint:design && tsc --noEmit -p tsconfig.build.json && vite build"
```
The `build`/CI `lint:design` step is the **blocking** enforcement (`exit 2` on any hash not in baseline).
A `settings.json` PostToolUse-style hook on `frontend/src/**/*.tsx` edits runs the guard **warn-only**
(prints new violations to stderr, `exit 0`) — it surfaces appearance-in-className at edit time without
false-blocking the author mid-refactor while they are lowering the baseline. (A hard PreToolUse block was
considered and rejected: at PreToolUse the on-disk file is pre-edit, so the hook would need to parse
proposed `new_string`/`content` from `tool_input` — more than the existing hooks do — and would fight the
ratchet PR. CI is the hard gate; the hook is the fast feedback.)

**Ratchet protocol:**
1. **Phase 0:** land tokens + guard + the seed baseline. CI fails on any violation hash *not* in baseline;
   nothing existing blocks, nothing new lands.
2. **Each adoption PR (Phases 1–3):** deletes its call sites and **must remove exactly those hashes from
   the baseline.** CI asserts the baseline set only ever shrinks (`baseline_new ⊆ baseline_old`); a PR that
   removes violations without shrinking the baseline is rejected.
3. **Per-rule hardening:** as soon as a *rule's* baseline entries reach zero, remove that rule's
   baseline-allowance branch so it becomes absolute, even while other rules still carry backlog. This
   prevents one stubborn straggler (the 96 `text-[` are the largest) from stranding the whole enforcement
   flip. When all rules are absolute, delete `design-baseline.json`.

Guard unit tests (fixtures): `text-[10px]` flagged; `text-base` clean; `rounded-[3px]` flagged;
`rounded-md` clean; `style={{ color }}` clean; `style={{ background: \`...${x}...\` }}` clean;
`"color-mix(in oklab, var(--color-accent) 9%, transparent)"` clean (var-only); `background: "oklch(…)"`
flagged; `border-l-4` flagged; `shadow-[…rgba…]` clean; a `ui/` path excluded; an allowlisted path's
`bg-[oklch]` clean; baseline subset-diff exit codes.

## 5. The catalog (static HTML)

A single hand-authored `docs/design-system.html` (sibling to `docs/redesign-mockups/*.html`, matching that
existing convention), built on the Print tokens, rendering every primitive's variants × sizes × states on
warm paper. It is a visual reference and adoption surface; the **guard + DESIGN.md remain the enforced
source of truth** (the HTML can drift and that's acceptable — it isn't load-bearing). Each primitive
section below lists the states the catalog must show. A static HTML page (over Storybook or a generated
gallery) is a deliberate choice: the enforcement comes from the guard, so the catalog can stay lightweight
documentation; the drift risk is accepted knowingly.

## 6. Primitives & migration waves

Each wave is its own TDD'd PR (failing test → minimal impl → verify → commit), keeps `npm run build` +
full Vitest green, and lowers the guard baseline by its deleted call sites. The exhaustive per-site
checklist for every primitive lives in `docs/2026-05-29-library-extraction-inventory.json`.

### Phase 0 — Foundation
Tokens (§3) + guard (§4) + catalog scaffold (§5) + the adopt-or-delete pass over the existing `ui/`
primitives (table below) + **the side-stripe fix, pulled forward** (it's an impeccable absolute ban *and*
the worst color-only-severity a11y failure, with no dependency on the primitive waves, so it ships first,
not last): convert `InfrastructureBanner.tsx:55` and `Toast.tsx:85` from `border-l-4` + hue-only severity
to a **leading status icon + word + full hairline border** (kind conveyed by icon, not edge hue). This
closes the audit's #1 AI-tell and the Toast a11y finding immediately (sign-off §8). The adopt-or-delete pass:

| Primitive | Verdict | Action |
|---|---|---|
| **Button** | ADOPT + fix | Rename variant axis `primary\|ghost\|danger` → `solid\|accent\|ghost\|danger`; `solid` = ink-fill + **`text-paper`** (fixes the lone `text-white` leak, `Button.tsx:11`); `accent` = terracotta. Migrate the **5 `variant="primary"` usages across 3 files** (SpeculativeBuilder, StaleIndicesBanner, OnboardingFlow) → `solid` (**terracotta→ink, sign-off §8**). Same PR rewrites the class-string assertions in `test/scoping.test.tsx:117-119` + `test/ScopingFoot.test.tsx:31-33` (they pin `bg-ink`/`text-paper`/`.not bg-accent` on the inline scoping confirm buttons) to `data-variant` assertions, and folds those two inline confirm buttons (`ScopingConfirmModal.tsx:53`, `ScopingFoot.tsx:43`) into `<Button variant="solid">` (preserving their testids). |
| **ScoreBar** | ADOPT by rewrite | `{ value, phase, size: "bar"\|"compact" }`; color from `phaseColor(phase)` internally; track `bg-hair`. Adopt the 2 inline phase bars (`PhasePanel:58` → `bar`, `PhasePanel:153` → `compact`). Exclude `SamplesPage:133` (it's a progress meter, not a score bar — future `Meter`). **Breaking API (`score`→`value`, `color`→`phase`)**: `test/ui/ScoreBar.test.tsx` is a real importer (asserts `style.background`/`score`/`color`) — it is a `tsc` error the moment the props change, which fails the per-PR `build` gate, so it **must be rewritten in the same commit** (assert `style.width` for clamping + a `data-phase`, not a raw color). |
| **Dot** | ADOPT by extend | Add `tone: accent\|success\|muted\|neutral` + `size: xs\|sm`; keep `role="img"`+`label`, allow `aria-hidden` decorative dots to skip the label. Adopt the 7 inline semantic dots. |
| **Toast** | ADOPT, fix barrel + stripe | Add `export { ToastContainer }` to `index.ts`; and (as part of the pulled-forward side-stripe work above) convert its `border-l-4` hue-only severity to a leading status icon + word + full hairline border, closing the color-only-severity a11y finding. |
| **HintText** | ADOPT as-is | Correct, 6 importers. No change. |
| **SectionLabel** | DELETE | 0 importers, no inline `text-label <h3>` sites to serve; subsumed by `Kicker`. |
| **Input** | DELETE | 0 importers; the real text-field cases (adornment, borderless-inline, steppers) need a designed `TextField` — a deferred later wave, not a trivial wrapper. |

Add `data-variant` to Button, `data-tone` to Dot, keep `data-score-bar` on ScoreBar — so tests assert
semantics, not classes. A barrel-surface test locks the public exports (catches the Toast omission / any
Input/SectionLabel re-add).

### Phase 1 — Highest-duplication primitives

**`SegmentedControl<T>`** — unifies 5 single-select segmented controls.
- API: `{ options: SegmentOption<T>[], value, onChange, role?: "group"|"radiogroup", variant?: "bordered"|"plain", size?, "aria-label" (required), testId? }`. `role` drives child semantics (`button[aria-pressed]` vs `radio[aria-checked]`).
- Canonical active state = **ink-on-paper** (`bg-ink text-paper`). Converges the two DESIGN.md "L-5/B-A defect" recessed fills (`AnnotationToggles`, `PlotCard`'s inline `XScaleToggle` — both `bg-paper-sunk text-ink`) onto it (**sign-off §8**).
- Carve-outs: `AnnotationToggles` is **multi-select** (two independent booleans) — does NOT fit `SegmentedControl<T>`; fix its banned `bg-paper-sunk` active fill → ink-on-paper in place this wave, defer a `ToggleButton` primitive. Document the in-place ink fill as the *interim canonical* multi-select treatment (and baseline-note it) so the guard/future readers don't mistake the corrected toggle for un-migrated drift with nowhere to migrate to. `SeriesFolioPage` filter chips are `rounded-full` pills → defer to a future `Chip`. `CorpusTopbar` view-seg mixes `<Link>`/`<span>`/disabled-button and is route-driven → leave as a bespoke nav switch borrowing the shared class map (not `onChange`-driven).
- Preserve exactly the **segment-level** contract the tests pin, not just the container: per-segment `data-active`/`data-value`, `role` (`radio` for radiogroup), and `aria-pressed`/`aria-checked` (`test/GroupingModeToggle.test.tsx:51-64`, `test/AnnotationToggles.test.tsx:45-72`), plus the container `data-mode` (`GroupingMode` + `SeriesBuilderPage.test.tsx:271-275`), `aria-current` (view-seg), and every `data-testid`. Fold in the 44px touch target gated behind `@media (pointer:coarse)` so dense desktop toolbars don't balloon (**sign-off §8**).
- Add WAI-ARIA roving arrow-key nav for the `radiogroup` variant (a11y upgrade; new behavior to test).

**`PhaseChip`** — the monospace phase-tinted data badge that always renders the phase name.
- API: `{ phase, variant?: "tint"|"solid", size?: "sm"|"md" }`; owns `phaseColor` + `color-mix`. Canonical tint = 13% fill (the M-6 value, mockup-derived), `rounded-sm`, mono bold; `solid` = phase fill + `text-paper`. The mockup HTML chips are **borderless** — the 35% hairline border is a *new, deliberate* addition for edge legibility on `plate` (sign-off §8), not a mockup value; pick it consciously rather than treating it as canonical.
- Adopt the 2 true chip sites: `SampleStatusChip` (stays as a domain wrapper owning the `phase|null` "Not indexed" empty state, delegates appearance to `PhaseChip`) and `MemberMetaRow` (`size="md"`). Adds a hairline border to `SampleStatusChip` (**minor sign-off §8**).
- **Explicitly NOT chips** (do not migrate — over-reach trap): `PhasePanel:141` candidate label and `:42` serif title (colored text), `FocusReflectionsTable:178` dot+label, `SpeculativeBuilder:243` "anchor" tag, `SeriesFolioCard` swatches. Reasons in the inventory.

**`PhaseStrip`** — consolidates `ScopingPhaseStrip` + `SeriesPhaseStrip` (near-duplicate files, both deleted).
- API: `{ segments: { phase: string|null; coexistWith?: string|null }[], size?, emptyLabel? }`; both callers map their domain shape into `segments`; the primitive derives its own caption.
- Canonical: unindexed cell = `ink-faint` (0.640) — the approved decision, folding onto DESIGN.md's placeholder/disabled token rather than adding the mockups' near-identical dedicated `--unindexed` (0.660). Scoping's pale `hair` was the shipped outlier (**sign-off §8, visibly darker**). Caption "throughout" uses the truthful distinct-count rule, not first===last (**behavior sign-off §8**); decorative arrow `aria-hidden`.
- **Second channel (partial now, rest deferred — be honest about which):** each segment gets `aria-label`/`title` naming its phase ("Pn3m", "Pn3m + Lamellar (coexistence)", "Unindexed"). `label` is itself an enumerated DESIGN.md second channel, so this **closes the AT/pointer channel**. But the audit's actual P1 — the *grayscale-glance* readability of adjacent segments for a sighted deuteranope — needs a **visual** glyph/pattern channel, which is **deferred to the plotting redesign** because it must co-design with the peak-shape encoding (don't ship two pattern vocabularies). State plainly: this wave partially discharges the Second-Channel requirement for strips; the visible half is consciously deferred. `buildPhaseStrip` stays as the Series data builder; the Scoping `heading` kicker moves to a sibling `<Kicker>` at the call site (drop `heading` from the primitive).

### Phase 2 — Chrome primitives

**`ModalShell`** — generalizes the proven `ConflictModalShell` over 7 hand-rolled overlays.
- API: `{ open, onClose, size?: "sm"|"md"|"lg", align?: "center"|"top", closeOnEsc?, closeOnOutsideClick?, variant?: "dialog"|"drawer", aria-* , testId? }`. Owns scrim (`bg-scrim`) + frame (`plate`/`hair-strong`/Plate-Lift) + Esc + focus-trap + outside-click; returns null when `!open`.
- **`useFocusTrap` fix** (load-bearing): add `textarea:not([disabled])` and `a[href]` to `FOCUSABLE` — unblocks the Notes drawer whose only focusable is a textarea (the empty-`focusable` early-return is *why* it currently leaks). This is a **global change to all 6 trap consumers** (NavModal, OnboardingFlow ×2, ConflictModalShell, SpeculativeBuilder, ScopingConfirmModal, Notes drawer), not a local Notes fix; none assert focusable-set size so risk is low, but ship the textarea-only trap test to lock it, and land the fix *with or before* any ModalShell drawer migration.
- **ModalShell must NOT impose an initial-focus target** — leave focus management to children. `NavModal` synchronously focuses its input on open; if the shell stole focus to the frame, that regresses (and the jsdom-synchronous focus a test relies on breaks).
- Per-site notes: `NavModal` passes `closeOnEsc={false}` (Esc lives in its input handler with `preventDefault` + Backspace-pops-chips); `OnboardingFlow` passes `closeOnOutsideClick={false}` (non-dismissable) and collapses its double-dialog into one frame — verify NameStep's Enter-submit still reaches its handler via bubbling once the per-step `role=dialog` is gone (TutorialStep already has `tabIndex=-1`); `SpeculativeBuilder` migration **fixes a real a11y bug** (`role=dialog` currently on the scrim); the Notes drawer uses `variant="drawer"` (lower-z scrim). Keep `ConflictModalShell` as a thin composition over `ModalShell` (preserves its test surface). **Modal frame radius = 5px** (follows the `.card` decision, resolving the inventory's open question; shipped `rounded-xl` 12px → `rounded-md` 5px).

**`Kicker`** — the uppercase-tracked eyebrow/label; subsumes `SectionLabel`. Per the approved decision it
gets its **own type role matching the mockups**, *not* `.text-label`: the mockup eyebrows render at
**weight 700**, `accent` tone at `tracking 0.13em`, `faint` tone at `tracking 0.09em` (≈11px). Collapsing
onto `.text-label` (500/0.06em) would visibly lighten every eyebrow, so Kicker defines its own role and
`tone: "accent"|"faint"` switches color. ~41 inline kickers collapse to it (normalizing the per-site
9/10/10.5/11px size + tracking drift onto the one role, preserving the 700 weight). `as?: "div"|"span"|"h2"|"h3"`
for the heading cases. Borderline sites flagged in the inventory:
table column-header strips (`SamplesPage:144`, `FocusReflectionsTable` headers — keep table semantics),
`SeriesRecipeEditor:131` (`text-ink` "Recipe" — needs a third `ink` tone or stays inline), and the
rotated heatmap axis title (placement-heavy but legal). Preserve `data-testid`s (`focus-plot-kicker`,
`series-builder-editing-badge`, `heatmap-axis-title`).

**`IconButton`** — standardizes 11+ icon-only/dismiss buttons.
- API: `{ label (required → aria-label), dismiss?, tone?: "ghost"|"danger", disabled?, ... }`. Standard
  44px hit area (`min-h/min-w`, glyph stays compact), the canonical accent `focus-visible` ring (only
  `PhasePanel:163` has it today), canonical dismiss glyph `×` (U+00D7). The grease-pencil reject **✕** on
  detector frames is a domain mark, **not** an IconButton. Fixes two sites with no `aria-label` (only
  `title`). Resolve the dismiss-hover convention (accent vs error) — lean accent for "remove chip",
  `danger`/error for destructive delete.

**`Card`/`Plate`** — `{ elevated? }` (default flat). `elevated` reuses the `.card` Plate-Lift recipe
(single source of the shadow); flat = `bg-plate` + hairline + **no shadow** (Flat-Except-the-Plate). Only
the four genuine lifted surfaces (figure plate, folio card, scoping worksheet, builder worksheet) use
`elevated`. The bg-plate grep over-captures: inputs, buttons, chips, banners, modals, the segmented
container, and floating objects (dock, popover, CullBar) route to their own primitives, not `Card` (the
inventory enumerates which). Adds `data-elevated`.

### Phase 3 — Sweep & ratchet to enforcing
Migrate the remaining inline kickers/plates/score bars/`text-[`/`rounded-[` onto primitives + the scale;
drain `design-baseline.json` to empty per rule, hardening each rule to absolute as its backlog clears (§4).
(The side-stripe ban was already fixed in Phase 0.)

## 7. Testing conventions

- Per-primitive Vitest asserting rendered semantics: roles, `aria-label`/`aria-checked`/`aria-pressed`,
  text content, and `data-*` attributes. **Never** assert Tailwind class strings (for new tests).
- **Existing class-string / leaked-prop tests to rewrite in the same PR as their primitive** (the spec's
  "stay green" guarantee does NOT hold for these): `test/ui/ScoreBar.test.tsx` (props change → tsc error,
  ships with the ScoreBar rewrite); `test/scoping.test.tsx:117-119` + `test/ScopingFoot.test.tsx:31-33`
  (`bg-ink`/`text-paper`/`bg-accent` assertions, ship with the Button-solid migration). Verify (likely no
  rewrite) `test/GroupingModeToggle.test.tsx` + `test/AnnotationToggles.test.tsx` still pass against the
  preserved segment-level `data-active`/`role`/`aria-*`.
- `useFocusTrap.test.tsx` gains a textarea-only trap case (locks the fix).
- All other migrated components' tests stay green against the same `data-testid`/role nodes.
- A `@theme` assertion test (mirroring `test/phases.test.ts`) pins `--color-scrim` present,
  `--color-print-accent === --color-accent` (resolved-value equality), `--radius` (singular) gone, AND the
  **resolved utility output** (a `rounded-md` element computes 5px) — not just custom-property presence.
- Guard fixture tests per §4.
- `npm run build` (now `lint:design` + tsc + vite) and full Vitest green each PR.

## 8. Decisions taken (all resolved)

**Mockup-vs-DESIGN forks, resolved with the user (the mockup-fidelity review surfaced these):**
- **Radius → 5px** for cards/plates/modals (match the mockups; DESIGN.md's 7px prose is now wrong and to be corrected). `.card` 12px and modal `rounded-xl` 12px both → 5px.
- **Unindexed → `ink-faint` (0.640)** — fold onto DESIGN.md's placeholder token, *not* the mockups' near-identical `--unindexed` (0.660).
- **Kicker → its own 700-weight role** (accent 0.13em / faint 0.09em), matching the mockups, *not* the lighter `.text-label` (500/0.06em).

**Other defaults (no sign-off needed):** bare `rounded` left alone; baseline **content-hash-keyed**;
guard CI-blocking + **PreToolUse hook warn-only**; **per-rule** enforcement flip; `ConflictModalShell`
kept as a composition; `Dot` two sizes; `PhaseChip` ships a `solid` variant + a deliberate 35% border;
touch-target via `@media (pointer:coarse)`; `SegmentedControl` single-select only
(AnnotationToggles/filter-chips/view-seg carved out); static-HTML catalog (drift accepted).

**Visible changes shipping (all approved, all DESIGN.md/mockup-correct):**
1. **Button confirm actions go terracotta → ink** (the `solid` variant). 5 usages / 3 files. (The mockups confirm primary buttons are ink, not terracotta; "+ New series" is ink-solid too, terracotta is reserved for the reject/grease-pencil mark.)
2. **`.card` + modal radius 12px → 5px** — touches the figure plate, the most prominent surface.
3. **PhaseStrip unindexed cells get darker** (`hair` → `ink-faint`) and the caption stops saying "X throughout" for non-monotone series (truthful distinct-count).
4. **`SampleStatusChip` gains a hairline border** (minor; mockup chips are borderless, this is a deliberate add).
5. **Touch targets grow on coarse pointers** (desktop unchanged).
6. **Toast + InfrastructureBanner lose the edge stripe**, gaining a leading status icon + word (Phase 0).

## 9. Build sequence

Phase 0 (tokens + guard + catalog + adopt/delete + side-stripe fix) → Phase 1 (`SegmentedControl`,
`PhaseChip`, `PhaseStrip`) → Phase 2 (`ModalShell`+focus-trap, `Kicker`, `IconButton`, `Card`/`Plate`) →
Phase 3 (sweep + per-rule ratchet to enforcing). Phase 0 is the only hard prerequisite for the rest;
Phases 1 and 2 primitives are independent of each other and can interleave; Phase 3 depends on all.
