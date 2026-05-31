# Greenfield UI Base Layer + Review-Fix Implementation Plan (Phase 1.5)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Act on the user's Storybook review of the 24 Phase-1 `src/print/ui/` primitives: extract a base-primitive layer (Input, Chip, Menu, Tooltip), refactor the specialized primitives onto it, add the missing "filled-out-properly" variants (Button `danger`, key-value tags), and clear the review's bug list — so Phase-2 page composition draws from a complete, consistent, DESIGN.md-faithful primitive set.

**Architecture:** Same closed-look idiom as the Phase-1 plan (`docs/superpowers/plans/2026-05-31-greenfield-ui-primitive-sweep.md` — read its "Conventions every task follows" section; it applies verbatim: local `cx()`, appearance internal, consumer `className` placement-only + last, `data-*` for every state, tests assert `data-*`/roles never class strings, per-primitive impeccable pass with `npx impeccable detect` + the checklist `docs/greenfield-rebuild/impeccable-primitive-checklist.md` + `lint:design` + `tsc -p tsconfig.build.json`). New base primitives live in `src/print/ui/`; specialized primitives refactor to compose them while PRESERVING their public `data-testid`/role/prop contract so existing `test/print-ui/` specs stay green.

**Standing principle (user, 2026-05-31):** The Print mockups (`docs/redesign-mockups/`) are a **reference/sketch, not absolute** — fill out every basic variant **properly per `DESIGN.md`**, not only what a mockup happened to show.

**Tech Stack:** React 18, TS strict + `exactOptionalPropertyTypes`, Tailwind v4 (`@theme` in `src/styles.css`), Vitest + RTL, Storybook v10. Guard: `scripts/check-design.mjs`.

**Working dir for all commands:** `packages/HimalayaUI/frontend` (docs at repo root `docs/`).

**Key DESIGN.md anchors:** §2 Status (Error `oklch 0.520 0.170 28`, status-only, sanctioned for destructive); §5 Inputs (recessed plate/paper-sunk well, hair border, 5px, focus→accent border + 1px ring); §5 Buttons (solid/accent/ghost + the danger we add); Semantic-Colour + Second-Channel + Flat-Except-the-Plate rules. Menu is plate-like → may carry Plate-Lift; Tooltip is the `frame-edge` dark exception.

---

## Group A — Review bug-fixes (independent; no base-layer dependency; do first)

### Task A1: TopBar wordmark + 700 weight
**Why:** Review note 5. The mockup wordmark is SANS, uppercase, weight 700, letter-spacing 0.16em ("HIMALAYA · SAXS") — my story used serif, and only weights 400/500/600 are `@import`ed so 700 falls back.
**Files:** Modify `src/styles.css` (add `@import "@fontsource/plus-jakarta-sans/700.css";` after the existing 600 import, line ~4); Create `src/print/ui/Wordmark.tsx` + `Wordmark.stories.tsx` + `test/print-ui/Wordmark.test.tsx`; Modify `src/print/ui/index.ts`; Modify `src/print/ui/TopBar.stories.tsx`.
- [ ] Test (`Wordmark.test.tsx`): renders the brand text; renders an optional `tail` (e.g. "SAXS") in a `[data-role="tail"]`; `data-testid="wordmark"`.
- [ ] Run `npm test -- Wordmark` → FAIL.
- [ ] Implement `Wordmark`: `interface WordmarkProps { children: ReactNode; tail?: ReactNode; className?: string }`. Render `<span data-testid="wordmark" className={cx("font-sans font-bold uppercase tracking-[0.16em] text-ink text-xs", className)}>{children}{tail && <span data-role="tail" className="text-ink-faint font-semibold"> · {tail}</span>}</span>`. (`tracking-[0.16em]` is an allowed in-`print/ui` arbitrary tracking for the brand mark — comment it; `text-xs` ≈ the mockup 12px wordmark.) Add `font-bold` requires the 700 import (this task adds it). Barrel: `export { Wordmark } from "./Wordmark";`.
- [ ] Run `npm test -- Wordmark` → PASS.
- [ ] Update `TopBar.stories.tsx`: replace the `font-serif` wordmark with `<Wordmark tail="SAXS">Himalaya</Wordmark>`.
- [ ] Impeccable pass (detect + checklist E: the wordmark is a SANS brand mark, not a serif title; F: tracking/size on-scale-or-commented) + guards.
- [ ] Commit: `feat(print-ui): Wordmark (sans/upper/700) + import jakarta 700 + fix TopBar story`.

### Task A2: Kicker tracking parity + drop duplicate story
**Why:** Review note 13. VERIFIED root cause: `.text-kicker` (styles.css:269) owns the size (`--text-sm`); there is NO font-size difference between tones. The accent reads "larger" because `.text-kicker-accent` (line 271) uses `letter-spacing: 0.13em` while `.text-kicker-faint` (line 272) uses `0.09em` — the wider tracking spreads it. Also default tone IS `faint`, so Default and Faint stories are identical.
**Files:** Modify `src/styles.css` (lines ~271–272); Modify `src/print/ui/Kicker.stories.tsx` (create if absent).
- [ ] Edit `src/styles.css`: unify the letter-spacing so the two tones differ ONLY by color. Set both `.text-kicker-accent` and `.text-kicker-faint` to the SAME `letter-spacing` (use `0.10em` as the shared value, between the current 0.13/0.09). Result: `.text-kicker-accent { color: var(--color-print-accent); letter-spacing: 0.10em; }` and `.text-kicker-faint { color: var(--color-ink-faint); letter-spacing: 0.10em; }`. Do NOT touch `.text-kicker`'s geometry.
- [ ] Add/adjust `Kicker.stories.tsx`: export `Accent` (tone accent) and `Faint` (tone faint) only — NO redundant `Default` (default===faint). Each renders sample eyebrow copy (e.g. "BEAMTIME 2024-03").
- [ ] Verify in the file that accent/faint now share both size and tracking (grep the rules; quote them in the report).
- [ ] Impeccable pass (detect on any touched `.tsx`; lint:design; tsc) + guards.
- [ ] Commit: `fix(print-ui): Kicker accent/faint size parity + dedupe story`.

### Task A3: Dot "Labeled" story fix
**Why:** Review note 11. `Dot`'s `label` is aria-only (a status disc has no visible text); the "Labeled" story misleads.
**Files:** Modify `src/print/ui/Dot.stories.tsx`.
- [ ] Replace the misleading `Labeled` story with `WithCaption`: a `render` placing the Dot beside VISIBLE text, e.g. `<span className="inline-flex items-center gap-1.5"><Dot tone="success" label="Kept" /> <span className="text-sm text-ink-soft">Kept</span></span>`, with a comment that `label` is the screen-reader name (the visible word is the sibling text). Keep the per-tone stories.
- [ ] Impeccable pass (story-only; just run lint:design + tsc + build check sanity) + commit: `fix(print-ui): Dot story shows label as aria + visible caption`.

### Task A4: Card `padding` prop
**Why:** Review note 3. Mockup padding lives in `.card-body` and varies (13–16px folio, 26–30px loupe).
**Files:** Modify `src/print/ui/Card.tsx`; Modify `test/print-ui/Card.test.tsx`; Modify `src/print/ui/Card.stories.tsx`.
- [ ] Test: `<Card padding="md">` sets `data-padding="md"` and applies a padding utility; default (no prop) sets NO `data-padding` and no padding (preserves current flexible behavior).
- [ ] Run → FAIL.
- [ ] Implement: add `padding?: "sm" | "md" | "lg"` to `CardOwnProps`. Module map `const paddingClass: Record<"sm"|"md"|"lg",string> = { sm: "p-3", md: "p-4", lg: "p-7" }` (≈13/16/28px). When set, append `paddingClass[padding]` to `appearance` and set `data-padding`. Default unset = no padding. Preserve all existing behavior (draft/elevated/border/as/type/rest).
- [ ] Run → PASS.
- [ ] Story: add `Padded` (`elevated padding="md"` with title+body content).
- [ ] Impeccable pass + commit: `feat(print-ui): Card optional padding prop`.

### Task A5: Hot PeakGlyph → darker outline stroke (no ring)
**Why:** Review note 14. The surrounding accent `<circle>` ring is disliked; use a darker/thicker outline on the glyph instead.
**Files:** Modify `src/print/ui/PeakGlyph.tsx`; check `src/print/ui/peakMark.ts` for the `ring` descriptor field; Modify `src/print/ui/PeakGlyph.stories.tsx`.
- [ ] Read `peakMark.ts` `peakGlyph()` + the `PeakGlyphDescriptor` (`ring`, `stroke`, `strokeWidth`). Decide: for the "hot" descriptor, instead of `ring`, paint the glyph polygon with a DARKER, THICKER outline (e.g. stroke `var(--color-accent)` at `strokeWidth + 1.5` over the normal halo). Keep the descriptor's `ring` field but stop rendering the surrounding `<circle>`; instead, when `ring` is truthy, render the polygon `stroke={var(--color-accent)}` with an increased `strokeWidth` (the "hot" emphasis as a heavier outline, not a halo ring).
- [ ] Edit `PeakGlyph.tsx`: remove (or gate off) the `<circle data-role="hot-ring">` block; in the polygon branch, when `ring` is set, use `stroke="var(--color-accent)"` and `strokeWidth={strokeWidth + 1.5}` (keep `data-role="peak-glyph"`; keep a `data-hot="true"` on the polygon so tests/legend can detect it). Preserve caret/diamond/triangle geometry.
- [ ] Update `PeakGlyph.stories.tsx` `Hot` story comment to describe the darker-stroke emphasis. Verify no other call site depends on the removed `hot-ring` (grep `hot-ring`).
- [ ] Run `npm test -- PeakGlyph` (and any peakMark test) → PASS; if a test asserted the ring `<circle>`, update it to assert the heavier stroke / `data-hot`.
- [ ] Impeccable pass + commit: `feat(print-ui): Hot peak glyph uses heavier outline stroke, not a ring`.

### Task A6: IconButton hover previews + correct the close-× semantics note
**Why:** Review note 12 (tones look identical at rest — true & mockup-faithful, differ only on hover) and note 9 (close/remove × is NEUTRAL per mockup, never accent — the seeded IconButton's "accent = chip-remove" comment is wrong).
**Files:** Modify `src/print/ui/IconButton.tsx` (comment only — do NOT change tone behavior); Modify `src/print/ui/IconButton.stories.tsx`.
- [ ] Fix the doc comment on `IconButtonTone`/the `accent` tone: accent is for the rationed terracotta INTERACTION mark (e.g. a primary inline action), NOT chip-remove; chip-remove/close uses `ghost` (neutral ink-faint→ink), per DESIGN.md §2 + mockup `.cs-x`/`.tag` remove (ink-faint→ink). Do not alter the toneClass map.
- [ ] Story: add a note (story description/parameters) that the three tones are intentionally identical at REST and differ only on hover; if `storybook-addon-pseudo-states` is NOT installed (check package.json — likely not), instead add a `HoverIntent` story rendering the three tones in a row with a caption explaining the hover colors, so the reviewer understands the distinction without a hover addon.
- [ ] Impeccable pass + commit: `docs(print-ui): correct IconButton close-x semantics + hover-intent story`.

---

## Group B — Base primitives (new)

### Task B1: Input — base recessed field
**Why:** DESIGN.md §5 Inputs. The shared base for SearchInput, the scoping order-field, note entry, the tag editor.
**Files:** Create `src/print/ui/Input.tsx` + `.stories.tsx` + `test/print-ui/Input.test.tsx`; barrel.
- **Interface:** `interface InputProps extends Omit<InputHTMLAttributes<HTMLInputElement>, "size"> { value: string; onValueChange: (v: string) => void; inputSize?: "sm" | "md"; leading?: ReactNode; trailing?: ReactNode; invalid?: boolean }`. (`onValueChange` gives the string directly; raw `onChange` still passes through via `...rest` if needed. `leading`/`trailing` are slot adornments — the search magnifier becomes a `leading`.)
- **Render:** a wrapper `<div data-testid="input" data-invalid={invalid?"true":undefined} className={cx("inline-flex items-center gap-2 bg-plate border rounded-sm transition-colors focus-within:border-accent", invalid ? "border-error" : "border-hair-strong", sizeClass[inputSize], className)}>` with `{leading}` then `<input value={value} onChange={(e)=>onValueChange(e.target.value)} className="flex-1 bg-transparent border-none outline-none text-base text-ink placeholder:text-ink-faint" {...rest} />` then `{trailing}`. `sizeClass = { sm: "px-2.5 py-1", md: "px-3 py-1.5" }`.
- **Tests:** reflects value; typing fires `onValueChange` with the new string; `data-invalid` only when invalid; renders leading/trailing slots; focus-within is the accent ring (assert structurally — the wrapper has `focus-within:border-accent`; test the behavior by asserting the input is focusable and the wrapper testid present — do NOT assert the class string; instead assert `data-invalid` toggling for the colored-border path).
- **Story:** `Default`, `WithLeadingIcon`, `Invalid`, `Small`.
- **Impeccable:** C (focus-within accent ring; invalid = error border + (consumer) message — note error border is second-channeled by the consumer's error text), F tokens, recessed plate well per §5.

### Task B2: Chip — base pill with variants
**Why:** Collapse FacetChip/FilterChip/TagPill geometry into one base; review notes 1/6/9.
**Files:** Create `src/print/ui/Chip.tsx` + `.stories.tsx` + `test/print-ui/Chip.test.tsx`; barrel (+ export `ChipVariant`).
- **Interface:** `type ChipVariant = "static" | "removable" | "add" | "toggle" | "trigger";` `interface ChipProps { variant?: ChipVariant; children?: ReactNode; active?: boolean; onClick?: () => void; onRemove?: () => void; className?: string }`. Default variant `"static"`.
- **Render** (one element per variant, shared pill base `pillBase = "inline-flex items-center gap-1 rounded-full text-xs font-semibold whitespace-nowrap transition-colors"`):
  - `static`: `<span data-testid="chip" data-variant="static" className={cx(pillBase, "px-2 py-px text-ink-soft bg-plate border border-hair", className)}>{children}</span>`.
  - `removable`: static + a trailing remove `<button type="button" aria-label="Remove" onClick={onRemove} className="text-ink-faint hover:text-ink focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-sm">×</button>`. **NEUTRAL × (hover→ink), per review note 9 / mockup.**
  - `add`: `<button type="button" data-variant="add" aria-label="Add" onClick={onClick} className={cx(pillBase,"px-2 py-px text-ink-faint border border-dashed border-hair-strong hover:text-accent hover:border-accent focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent", className)}>{children ?? "+ add"}</button>`. (Accent-on-hover IS correct for the add invite — mockup `.tag-add`.)
  - `toggle`: `<button type="button" data-variant="toggle" aria-pressed={active} data-active={active?"true":"false"} onClick={onClick} className={cx(pillBase,"px-3 py-1 border focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent", active?"bg-ink text-paper border-ink":"bg-plate text-ink-soft border-hair-strong hover:border-ink-faint", className)}>{children}</button>`.
  - `trigger`: `<button type="button" data-variant="trigger" aria-haspopup="menu" onClick={onClick} className={cx(pillBase,"px-3 py-1 text-ink font-semibold bg-plate border border-hair-strong hover:bg-paper-sunk focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent", className)}>{children}<span aria-hidden className="text-ink-faint">▾</span></button>`.
- **Tests:** each variant sets `data-variant`; removable's × has `aria-label="Remove"` + fires onRemove; toggle reflects `aria-pressed`/`data-active` + fires onClick; trigger has `aria-haspopup="menu"`; add fires onClick.
- **Story:** one per variant (Static, Removable, Add, ToggleOff, ToggleOn, Trigger).
- **Impeccable:** A (toggle = ink/paper inversion not hue); B (accent only on add-hover + focus rings; remove × neutral); C; F.

### Task B3: Menu — plate dropdown popover (net-new)
**Why:** Review notes 1/9; not in mockups — design to the system.
**Files:** Create `src/print/ui/Menu.tsx` + `.stories.tsx` + `test/print-ui/Menu.test.tsx`; barrel (+ `MenuOption`).
- **Interface:** `interface MenuOption<T extends string> { value: T; label: ReactNode; disabled?: boolean }` `interface MenuProps<T extends string> { open: boolean; options: ReadonlyArray<MenuOption<T>>; onSelect: (v: T) => void; onClose: () => void; "aria-label": string; activeValue?: T; className?: string }`.
- **Render:** when `open`, a `<div role="menu" aria-label={ariaLabel} data-testid="menu" className={cx("card absolute z-20 mt-1 min-w-[180px] py-1", className)}>` (`.card` = the Plate-Lift surface — Menu is plate-like, allowed). Each option `<button role="menuitem" type="button" disabled={o.disabled} data-value={o.value} data-active={o.value===activeValue?"true":undefined} onClick={()=>{onSelect(o.value); onClose();}} className={cx("flex w-full items-center px-3 py-1.5 text-sm text-left transition-colors", o.value===activeValue?"text-ink bg-paper-sunk":"text-ink-soft hover:text-ink hover:bg-paper-sunk", "disabled:opacity-40 disabled:cursor-not-allowed focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent")}>{o.label}</button>`. Keyboard: Escape calls `onClose`; ArrowUp/Down move focus among `menuitem`s (roving focus via refs, mirroring SegmentedControl's radiogroup pattern). Render nothing when `!open`.
- **Tests:** renders N `menuitem`s when open; nothing when closed; clicking an option fires `onSelect(value)` then `onClose`; Escape fires `onClose`; ArrowDown moves focus to the next item; disabled option doesn't fire.
- **Story:** `Open` (a few options, one active), `WithDisabled`. (Render inside a `relative` decorator so the absolute popover sits correctly.)
- **Impeccable:** H (Menu is plate-like → Plate-Lift is allowed, note it); C (menuitem hover/focus/active); keyboard a11y (role=menu/menuitem + Escape + arrows).

### Task B4: Tooltip — caption popover (net-new)
**Why:** Review "other" note; mockups only use browser `title`.
**Files:** Create `src/print/ui/Tooltip.tsx` + `.stories.tsx` + `test/print-ui/Tooltip.test.tsx`; barrel.
- **Interface:** `interface TooltipProps { label: string; children: ReactElement; side?: "top" | "bottom"; className?: string }` (wraps a single focusable/hoverable trigger child; shows on hover + focus).
- **Render:** a `<span className="relative inline-flex">` wrapping `children` (with `aria-describedby` wired to the tip id) + a `<span role="tooltip" id={id} data-testid="tooltip" data-side={side} className={cx("pointer-events-none absolute z-30 whitespace-nowrap rounded-sm px-2 py-1 text-xs", side==="top"?"bottom-full mb-1":"top-full mt-1", className)} style={{ background: "var(--color-frame-edge)", color: "var(--color-frame-tag)" }}>{label}</span>`. Visibility: track `open` via `onMouseEnter/Leave` + `onFocus/Blur` on the wrapper; render the tooltip span only when open (so the test can assert show/hide). Honor `prefers-reduced-motion` is automatic (global rule); no entrance animation needed beyond opacity.
- **Tests:** tooltip not in DOM at rest; appears on `fireEvent.mouseEnter`/`focus` of the trigger; hides on leave/blur; has `role="tooltip"` + the trigger gets `aria-describedby` matching the tip id; renders the label.
- **Story:** `OnButton` (wraps a `<Button>`), `Bottom`.
- **Impeccable:** Tooltip is the `frame-edge` dark exception (DESIGN.md §2 Frame-edge / §4 detector dark window) — dark caption on the warm field, `frame-tag` text; note it. A11y: role=tooltip + aria-describedby. G: no bounce.

---

## Group C — Button danger + refactor specialized onto base

### Task C1: Button `danger` variant (proper, per DESIGN.md §2 Status)
**Files:** Modify `src/print/ui/Button.tsx`; Modify `test/print-ui/Button.test.tsx`; Modify `Button.stories.tsx`.
- **Behavior:** `danger` = error-red text at REST (quiet but identifiable), error-red FILL + paper text on hover/active (committed destructive moment), accent focus ring. Replace the current `danger` variantClass:
  `danger: "text-error border border-transparent hover:bg-error hover:text-paper hover:border-error focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"`.
- [ ] Test: a `variant="danger"` button has `data-variant="danger"`; renders its label; fires onClick. (Behavior/contract only — not color strings.) Keep existing armed/variant tests.
- [ ] Run → PASS (the data-variant already flows; the test mainly locks the variant name + onClick).
- [ ] Story: update `Danger` story to show the destructive read (label "Delete sample"); keep it.
- [ ] Impeccable: B/§2 — red is status/destructive, rationed, second-channeled by the label; confirm `--color-error`/`bg-error`/`text-error`/`border-error` utilities resolve (grep `--color-error` in styles.css; if the Tailwind utility `bg-error` isn't generated, use the `@theme` token name that is — verify against an existing `text-error` usage, which the seed Button already used).
- [ ] Commit: `feat(print-ui): Button danger variant reads destructive (error-red, per DESIGN.md §2)`.

### Task C2: Refactor SearchInput onto Input (+ optional suggestions via Menu)
**Files:** Modify `src/print/ui/SearchInput.tsx`; keep `test/print-ui/SearchInput.test.tsx` GREEN (same `data-testid="search-input"` + placeholder + value + onChange contract); Modify story.
- [ ] Reimplement `SearchInput` to render `<Input data-testid="search-input" leading={<MagnifierSvg/>} value={value} onValueChange={onChange} placeholder={placeholder} ... />`. Preserve the existing prop signature (`value/onChange/placeholder/className`) and the `data-testid="search-input"` (pass it through to Input's wrapper, OR wrap). The magnifier svg moves to the `leading` slot. Keep the icon `aria-hidden`.
- [ ] OPTIONAL suggestions (keep minimal): add `suggestions?: string[]` + `onPick?: (s: string) => void`; when non-empty render a `Menu` below. If this complicates the green-test contract, defer suggestions to Phase 2 and just do the Input refactor — note the decision in the report.
- [ ] Run `npm test -- SearchInput` → PASS (unchanged contract). Update the story to match.
- [ ] Impeccable pass + commit: `refactor(print-ui): SearchInput composes Input (+ optional Menu suggestions)`.

### Task C3: Refactor FilterChip onto Chip
**Files:** Modify `src/print/ui/FilterChip.tsx`; keep `test/print-ui/FilterChip.test.tsx` green (`data-testid="filter-chip"`, `aria-pressed`, `data-active`, onClick, label).
- [ ] Reimplement `FilterChip` as a thin wrapper: `<Chip variant="toggle" active={active} onClick={onClick} className={className} data-testid="filter-chip">{label}</Chip>` — BUT the Chip's toggle sets `data-testid="chip"`; to preserve the FilterChip contract, either (a) have Chip forward an optional `testId` prop, or (b) keep FilterChip rendering its own button but using Chip's class recipe via a shared exported helper. Prefer (a): add an optional `testId?: string` to `ChipProps` that overrides `data-testid` (default per-variant `"chip"`). Then FilterChip passes `testId="filter-chip"`. Update Chip + its test accordingly.
- [ ] Run `npm test -- FilterChip Chip` → PASS.
- [ ] Impeccable pass + commit: `refactor(print-ui): FilterChip composes Chip toggle variant`.

### Task C4: Refactor FacetChip onto Chip trigger + wire Menu
**Files:** Modify `src/print/ui/FacetChip.tsx`; keep `test/print-ui/FacetChip.test.tsx` green (`data-testid="facet-chip"`, label, onClick, aria-hidden chevron). 
- [ ] Reimplement `FacetChip` to use `<Chip variant="trigger" testId="facet-chip" onClick={onClick}>{label}</Chip>` (the trigger variant already provides the chevron + `aria-haspopup`). OPTIONAL: add `options?` + `onSelect?` so FacetChip manages an internal `open` state and renders a `Menu`; if it complicates the green contract, keep FacetChip trigger-only and leave Menu wiring to the consumer (note the decision). Preserve the existing prop signature.
- [ ] Run `npm test -- FacetChip` → PASS.
- [ ] Impeccable pass + commit: `refactor(print-ui): FacetChip composes Chip trigger variant`.

### Task C5: ScreenedMark → CheckCircle (generalize)
**Why:** Review "other" note — repurpose the 13px disc+check as a multi-select selection check (series creation) AND screened-status. Mockup `.screened-mark` confirms the visual.
**Files:** Rename `src/print/ui/ScreenedMark.tsx` → `src/print/ui/CheckCircle.tsx`; rename test + story; update barrel; grep for any `ScreenedMark` usage.
- [ ] `git mv` the three files (ScreenedMark.tsx/.stories.tsx, test/print-ui/ScreenedMark.test.tsx) to CheckCircle.
- [ ] Reimplement as `interface CheckCircleProps { checked: boolean; label?: string; className?: string }` — same 13px disc + paper check when `checked`, hairline ring when not. `aria-label` defaults to `checked ? "Selected" : "Not selected"` but a `label` prop overrides (so a screened-status use can pass `label={done?"Screened":"Not screened"}`). `data-checked` when checked. Keep `role="img"`.
- [ ] Update tests to the `checked`/`label` API (was `done`); update the barrel `export { CheckCircle } from "./CheckCircle";` (remove the ScreenedMark export). Grep `ScreenedMark` across `src/print` + `test/print-ui` → no stragglers.
- [ ] Run `npm test -- CheckCircle` → PASS.
- [ ] Impeccable pass (A: selection carries fill+check+aria-label, not hue) + commit: `refactor(print-ui): ScreenedMark -> CheckCircle (multi-select + screened)`.

---

## Group D — Key-value tags

### Task D1: Tag model + TagPill (key faint sans + value mono ink), composing Chip
**Files:** Create `src/print/ui/tag.ts` (the `Tag` type); Modify `src/print/ui/TagPill.tsx`; Modify `test/print-ui/TagPill.test.tsx`; Modify story; barrel (`export type { Tag }`).
- **Type:** `export interface Tag { key: string; value?: string }`.
- **TagPill interface:** `interface TagPillProps { tag: Tag; onRemove?: () => void; className?: string }` (was `children`). Render via `Chip`: key-only → `<Chip variant={onRemove?"removable":"static"} onRemove={onRemove} testId="tag-pill"><span data-role="tag-key">{tag.key}</span></Chip>`; key+value → same but content is `<span data-role="tag-key" className="text-ink-soft">{tag.key}</span><span data-role="tag-value" className="font-mono text-ink ml-1">{tag.value}</span>`. (Key faint/soft sans label + value mono ink — the chosen rendering.) Set `data-has-value={tag.value?"true":undefined}`.
- **Tests:** key-only renders just the key, `data-has-value` absent; key+value renders key in `[data-role="tag-key"]` and value in `[data-role="tag-value"]` (mono), `data-has-value="true"`; remove × present only with onRemove, neutral, fires onRemove; `data-testid="tag-pill"`.
- **Story:** `KeyOnly` (`tag:{key:"LL37"}`), `KeyValue` (`tag:{key:"temperature", value:"37C"}`), `Removable`.
- **Impeccable:** E (value = mono measurement, key = sans label — the chosen voice split); B (neutral remove ×); composes Chip (DRY).

### Task D2: TagEditor — key + optional value entry
**Files:** Create `src/print/ui/TagEditor.tsx` + `.stories.tsx` + `test/print-ui/TagEditor.test.tsx`; barrel.
- **Interface:** `interface TagEditorProps { onCommit: (tag: Tag) => void; onCancel?: () => void; knownKeys?: string[]; className?: string }`.
- **Render:** an inline row using `Input`: a key `Input` (with an optional `Menu` of `knownKeys` as a `trailing` trigger or a datalist-style suggestion — keep it simple: a key Input + an optional value Input), a value `Input` (optional), and a commit affordance (Enter or a small `Button`). On commit, call `onCommit({ key, value: value || undefined })` (omit empty value → key-only tag). Escape → `onCancel`. Manage local `key`/`value` state.
- **Tests:** typing a key and pressing Enter (or clicking commit) fires `onCommit({key})` (no value → key omitted); typing key+value commits `{key, value}`; empty key does NOT commit; Escape fires `onCancel`.
- **Story:** `Default`, `WithKnownKeys`.
- **Impeccable:** composes Input; C (focus rings via Input); value optional → key-only path works.

### Task D3: TagList → Tag[] + Chip add invite + TagEditor
**Files:** Modify `src/print/ui/TagList.tsx`; Modify `test/print-ui/TagList.test.tsx`; Modify story.
- **Interface:** `interface TagListProps { tags: Tag[]; onAdd?: (tag: Tag) => void; onRemove?: (tag: Tag) => void; editable?: boolean; className?: string }` (tags now `Tag[]`).
- **Render:** map `tags` → `<TagPill tag={t} onRemove={editable&&onRemove?()=>onRemove(t):undefined} />` (conditional spread for the optional, per exactOptionalPropertyTypes). The add affordance: a `Chip variant="add"` that toggles an inline `TagEditor` (local `adding` state); on `TagEditor.onCommit`, call `onAdd(tag)` and close. Keep the reveal-on-hover behavior for the add chip when tags exist (group-hover opacity), always-visible when empty.
- **Tests:** renders N TagPills (keyed by `key` — assume unique keys; if dup keys possible, key by index); add chip present with `onAdd`, opens the editor, committing fires `onAdd(tag)`; editable+onRemove makes pills removable; empty list shows the add invite.
- **Story:** `WithTags` (mix of key-only + key+value), `Empty`, `Editable`.
- **Impeccable:** composes TagPill + Chip + TagEditor (DRY); G (opacity reveal only).

---

## Group E — Storybook review refresh

### Task E1: Update index, whole-set gate, rebuild, relaunch
**Files:** Modify `docs/greenfield-rebuild/primitive-storybook-index.md`.
- [ ] Confirm every primitive (incl. new Input, Chip, Menu, Tooltip, Wordmark, TagEditor, CheckCircle) has a story (the `comm` check from the Phase-1 plan Task 26).
- [ ] `npx impeccable detect src/print/ui/` → exit 0 across the set.
- [ ] Full `test/print-ui/` run green; `npm run build-storybook` PASS; `lint:design` + `tsc -p tsconfig.build.json` exit 0.
- [ ] Update the index doc: add the base primitives + their variants, the danger button, key-value tags, CheckCircle; note the refactors (FacetChip/FilterChip/SearchInput/TagPill now compose base primitives). Add a "Phase-1.5 changes" section.
- [ ] Relaunch Storybook (`npm run storybook -- --ci`) and report the URL.
- [ ] Commit: `docs(rebuild): Storybook index — Phase-1.5 base layer + refactors`.

---

## Self-Review

**1. Coverage of the user's 14+ notes:** 1 FacetChip dropdown → B3 Menu + C4. 2 danger → C1. 3 Card padding → A4. 4 >2 phases → documented (binary, domain max 2; no task needed). 5 TopBar fonts → A1. 6 TagPill dashed add → B2 `add` variant + D3. 7 tag key/value → D1–D3. 8 GripHandle → kept (no task). 9 search suggestions → C2 (optional Menu); × color → B2 neutral remove + A6 comment. 10 ToggleSwitch no-label → ToggleSwitch already requires `label` as `aria-label`; ADD a task if a visually-label-less variant is needed — **note:** add `hideLabel?: boolean` to ToggleSwitch (visually hide the text, keep aria) as a small step folded into A-group if time; otherwise track. 11 Dot label → A3. 12 IconButton identical → A6 (faithful; hover story). 13 Kicker → A2. 14 Hot glyph → A5. ScreenedMark repurpose → C5. Base primitives + refactor → B + C2/C3/C4. ✓
**2. Placeholders:** Each task gives interface + render + tests + story + impeccable notes. The two OPTIONAL items (SearchInput suggestions, FacetChip menu wiring) explicitly allow deferral with a reported decision — not vague TODOs.
**3. Type consistency:** `Tag` (D1) consumed by TagPill/TagEditor/TagList. `ChipVariant`/`testId` (B2) consumed by FilterChip/FacetChip (C3/C4). `MenuOption` (B3) by FacetChip/SearchInput. `Input`'s `onValueChange` used by SearchInput/TagEditor. `CheckCircle` API (`checked`/`label`) replaces ScreenedMark `done`. Refactors PRESERVE existing `data-testid`s so Phase-1 tests stay green (Chip `testId` override is the mechanism).
**4. Ordering:** A (independent fixes) first. B base primitives before C/D refactors that compose them (Input→SearchInput/TagEditor; Chip→FilterChip/FacetChip/TagPill; Menu→FacetChip/SearchInput). D after B (TagPill needs Chip; TagEditor needs Input). E last.
**5. ToggleSwitch no-label (note 10):** fold a `hideLabel?: boolean` step into the ToggleSwitch as an A-group addendum (A7) if building; keep aria-label intact.

### Task A7: ToggleSwitch optional visually-hidden label (note 10)
**Files:** Modify `src/print/ui/ToggleSwitch.tsx` + test + story.
- [ ] Add `hideLabel?: boolean`; when true, the text `<span>` is visually hidden (still the aria-label + screen-reader text) so the switch can stand alone in a dense toolbar. Default false. NOTE: `sr-only` is NOT defined in styles.css — it's a Tailwind v4 built-in. First confirm `sr-only` actually hides (quick check / build); if it does not generate, hide via an inline style object `{ position:"absolute", width:1, height:1, padding:0, margin:-1, overflow:"hidden", clip:"rect(0 0 0 0)", whiteSpace:"nowrap", border:0 }` (the canonical visually-hidden recipe) applied when hideLabel.
- [ ] Test: with `hideLabel`, the switch still has `role="switch"` + `aria-label={label}` and the visible text span is `sr-only` (assert the span has the sr-only class via a `data-role` or that the accessible name still resolves). Keep existing tests green.
- [ ] Story: `NoVisibleLabel`.
- [ ] Impeccable pass + commit: `feat(print-ui): ToggleSwitch hideLabel (icon-dense standalone) `.
