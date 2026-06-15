# Greenfield Batch 3 — finish the tier-1 composites

> **For agentic workers:** executed via superpowers:subagent-driven-development (fresh implementer + review per task). Steps use checkbox (`- [ ]`).
> Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Gates per task: the new `test/print-*` test green · `npm run lint:design` clean · `npx tsc --noEmit -p tsconfig.build.json` clean. Commit only the named files (never `git add -A`).

**Goal:** Build the 5 remaining Layer-2 tier-1 composites (`FolioHeader`, `AutoGroup`, `PhaseBlock`, `ReadingRow`, `MemberRow`), plus the 2 refactor-on-contact primitives they require (`Swatch`, `PhaseLabel`). This closes the primitive-only composite layer and unlocks the downstream Layer-3 panels (AssignmentCart, MemberList, ScopePlate, ReadingPanel, Gallery).

**Architecture:** Bottom-up. The two new primitives land first (appearance authored in the design-guard-exempt `src/print/ui/`), then each composite composes primitives + placement/named-token utilities only. `src/print/components/**` is NOT guard-exempt → composites carry zero inline appearance literals.

**Tech stack:** React + TypeScript (strict), Tailwind v4 (named-token + type-role utilities), Vitest + RTL, Storybook (CSF, auto-globbed from `src/print/**/*.stories.tsx`).

---

## The load-bearing discovery (why the 2 primitives)

Three composites need to render a **phase color** (`phaseColor(phase)` → an OKLCH string). The design guard (`scripts/check-design.mjs`) is a *source-text* scanner: rule #5 flags a literal `oklch(`/`rgba(`/hex on a line (after stripping `var(--color-*)`), so an inline `style={{ color: phaseColor(x) }}` would technically *pass* (it's a runtime call, no literal). **But that leaks appearance authoring into the placement-only layer** and breaks the closed-look contract. The established pattern proves the right move: `ScoreBar`/`PhaseChip` (both in exempt `ui/`) call `phaseColor()` internally and paint via inline `style`; the composite `CandidateRow` only passes `phase={...}`. So phase-colored dots/swatches/text get a primitive.

- **`Swatch`** — phase-colored color sample (square for MemberRow, circle for ReadingRow). New primitive.
- **`PhaseLabel`** — phase-colored text wrapper (PhaseBlock name + score, ReadingRow phase name). New primitive.
- `Dot` is intentionally left untouched (its `tone` enum stays status-pure; a circular color sample is a `Swatch`, not a status Dot — and this avoids regression on Dot's existing tests).

## Verified facts (do not re-derive)

- `phaseColor(phase: string): string` — `src/phases.ts:45`. Returns an OKLCH string.
- Accent CSS var: `--color-print-accent` (aliases `--color-accent`). Inline `style`/SVG `fill="var(--color-print-accent)"` is guard-clean (rule #5 strips `var(--color-*)`; SVG `fill="…"` is not a `fill-[…]` utility so rule #3 does not match).
- Type roles (`src/styles.css`): `--text-xs 10.5` · `--text-sm 11.5` · `--text-base 13` · `--text-lg 15` · `--text-xl 19` · `--text-display 27` (serif) · `--text-headline-lg 26` (serif). Role classes: `text-display` (serif 27), `text-headline` (serif 19), `text-headline-lg` (serif 26), `text-body` (13), `text-caption` (xs), `text-kicker` + `text-kicker-accent`/`-faint`, `text-meta` (sm), `text-label` (sm), `text-data` (mono sm), `text-data-strong` (mono sm badge).
- Token utilities: `bg-paper` `bg-paper-sunk` `bg-plate` · `text-ink` `text-ink-soft` `text-ink-faint` · `border-hair` `border-hair-strong` · `text-print-accent`/`bg-print-accent`. Radius collapsed to one step: use `rounded`/`rounded-sm` (both 5px); circles use `rounded-full`. Layout/size arbitraries (`h-[9px]`, `w-[…]`, `gap-*`, inline `style` for layout) are allowed in composites; appearance arbitraries (`text-[…]`, `rounded-[…]`, raw `oklch()`/hex, `border-l/r-2..9`) are NOT.

### Verified primitive signatures

- `Kicker({ tone?: "accent"|"faint" = "faint", as?: "div"|"span"|"h2"|"h3", className?, children })`
- `IconButton({ label: string /*required, aria-label*/, tone?: "ghost"|… = "ghost", dismiss?: boolean, children?, …button })`
- `GripHandle({ className? })` — renders `⋮⋮`, color shifts on `group-hover` (parent row must be a `group`).
- `PhaseChip({ phase: string, coexistWith?: string[], variant?: "tint"|"solid" = "tint", size?: "sm"|"md", className? })`
- `ScoreBar({ value: number, phase: string, size?: "bar"|"compact" = "bar", className? })`
- `Dot({ tone?, size?: "xs"|"sm"|"md", bordered?, label?, …span })`

### House conventions

- Tests assert via `data-*` / roles / text — NEVER class strings (`test/AGENTS.md`).
- Composite file pattern: `data-testid` on root; `className` is appended-last placement-only; document it `/** PLACEMENT-ONLY. */`.
- Colocate stories: `src/print/<…>/<Name>.stories.tsx`. Unit tests: `test/print-ui/<Name>.test.tsx` (primitives) / `test/print-components/<Name>.test.tsx` (composites).
- House-style references to read: `src/print/components/CandidateRow.tsx`, `PlateHeader.tsx`; primitive: `src/print/ui/ScoreBar.tsx` (the phaseColor-via-style pattern), `Dot.tsx`.

---

## Task 1 — `Swatch` primitive (phase-colored color sample)

**Files:** Create `src/print/ui/Swatch.tsx`, `src/print/ui/Swatch.stories.tsx`, `test/print-ui/Swatch.test.tsx`. Modify `src/print/ui/index.ts` (export).

Mirror `ScoreBar`'s phaseColor-via-`style` pattern exactly. A small color sample that fills with the phase's color.

```tsx
import { phaseColor } from "../../phases";

export type SwatchShape = "square" | "circle";

interface SwatchProps {
  /** Phase whose color fills the swatch (via phaseColor). */
  phase: string;
  /** "square" (rounded color chip, default) | "circle" (reading-row dot). */
  shape?: SwatchShape;
  /** PLACEMENT-ONLY. */
  className?: string;
}

const shapeClass: Record<SwatchShape, string> = {
  square: "rounded-sm",
  circle: "rounded-full",
};

export function Swatch({ phase, shape = "square", className = "" }: SwatchProps): JSX.Element {
  return (
    <span
      data-swatch
      data-phase={phase}
      aria-hidden="true"
      className={`inline-block shrink-0 h-[9px] w-[9px] ${shapeClass[shape]} ${className}`}
      style={{ background: phaseColor(phase) }}
    />
  );
}
```

- [ ] **Step 1 — failing test** (`test/print-ui/Swatch.test.tsx`): renders `[data-swatch]`; `data-phase` equals the given phase; default shape carries `rounded-sm`, `shape="circle"` carries `rounded-full` (assert via `className.includes` is OK for *shape geometry only* — but prefer asserting `data-phase` + that the inline `style.background` equals `phaseColor(phase)` exactly, imported from `../../src/phases`); `aria-hidden="true"` present. Do NOT assert color class strings.
- [ ] **Step 2** run → fail.
- [ ] **Step 3** implement `Swatch.tsx` as above; add `export { Swatch } from "./Swatch";` (+ type) to `index.ts`.
- [ ] **Step 4** stories: `AllShapes` (square + circle across 3–4 phases from fixtures), plus a default story.
- [ ] **Step 5** gates: `npm test -- print-ui/Swatch` · `npm run lint:design` · `npx tsc --noEmit -p tsconfig.build.json`.
- [ ] **Step 6** commit `src/print/ui/Swatch.tsx Swatch.stories.tsx src/print/ui/index.ts test/print-ui/Swatch.test.tsx`.

> Scope note: gradient-for-coexistence + form-factor dashed-outline modes are **deferred** to the later series-plot `.member` work — do NOT build them now. Leave a one-line `// later: gradient/coexistence mode (series-plot member)` comment.

---

## Task 2 — `PhaseLabel` primitive (phase-colored text)

**Files:** Create `src/print/ui/PhaseLabel.tsx`, `.stories.tsx`, `test/print-ui/PhaseLabel.test.tsx`. Modify `index.ts`.

Colors its children with the phase color; the consumer sizes/weights/fonts it via placement `className`. (PhaseBlock uses it twice — serif 22px name, mono 12px score; ReadingRow uses it for the 12px bold name.)

```tsx
import type { ReactNode } from "react";
import { phaseColor } from "../../phases";

interface PhaseLabelProps {
  /** Phase whose color tints the text (via phaseColor). */
  phase: string;
  as?: "span" | "div";
  /** PLACEMENT-ONLY (font/size/weight live here, e.g. text-display, text-data font-bold). */
  className?: string;
  children: ReactNode;
}

export function PhaseLabel({ phase, as: Tag = "span", className = "", children }: PhaseLabelProps): JSX.Element {
  return (
    <Tag data-phase-label data-phase={phase} className={className} style={{ color: phaseColor(phase) }}>
      {children}
    </Tag>
  );
}
```

- [ ] **Step 1 — failing test**: renders children text; `[data-phase-label]` present with `data-phase` == phase; inline `style.color` === `phaseColor(phase)`; `as="div"` renders a `div`. No class-string color assertions.
- [ ] **Step 2** fail · **Step 3** implement + export from `index.ts` · **Step 4** stories (`AllPhases`, and one showing it sized as a serif title vs a mono value) · **Step 5** gates · **Step 6** commit the 4 paths.

---

## Task 3 — `FolioHeader` composite

**Files:** Create `src/print/components/FolioHeader.tsx`, `.stories.tsx`, `test/print-components/FolioHeader.test.tsx`.
**Mockup:** `docs/redesign-mockups/series-folio.html` (the `.head` block) — kicker + serif title + subtitle on the left, a right-aligned count block.

Structure: root `flex items-end justify-between gap-8`; left column = `<Kicker tone="accent">` + serif `<h1 className="text-display">` + `<p className="text-body text-ink-soft mt-2 max-w-[60ch]">`; right `shrink-0 text-right` = count number `<div className="text-headline-lg text-ink">` + label `<div className="text-caption text-ink-faint mt-0.5">` (uppercase via existing utility / `uppercase` layout class is fine).

Props: `{ kicker: string; title: string; subtitle?: ReactNode; count: number; countLabel: string; className?: PLACEMENT-ONLY }`.

> Fidelity snap (documented): the mockup title is 31px and the count 26px. The named scale has no 31px step, and `text-[31px]` is guard-banned. Use `text-display` (27px) for the title (nearest named role) and `text-headline-lg` (26px — exact) for the count. Note this in a code comment; if fidelity review objects, the fix is a new `--text-*` scale step in `styles.css`, not an arbitrary.

- [ ] **Step 1 — failing test**: root `data-testid="folio-header"`; kicker text, title text (in an `h1`), subtitle text, count number, and count label all render and are findable by role/text; subtitle omitted → not rendered. No class assertions beyond presence/role.
- [ ] **Step 2** fail · **Step 3** implement (compose `Kicker` from `../ui`) · **Step 4** stories (`Default` with subtitle + count; `NoSubtitle`) against fixtures · **Step 5** gates (`npm test -- print-components/FolioHeader`) · **Step 6** commit 3 paths.

---

## Task 4 — `AutoGroup` composite

**Files:** Create `src/print/components/AutoGroup.tsx`, `.stories.tsx`, `test/print-components/AutoGroup.test.tsx`.
**Mockups:** `docs/redesign-mockups/series-scoping.html` (`.autogroup` summary variant — recessed `bg-paper-sunk`, no actions) AND `series-builder.html` (`.autogroup` compose variant — `bg-plate`, title + link-button actions).

Confident-expert summary callout: a star icon (★) + body text (with bolded emphasis spans) + optional title + optional action link-buttons.

Structure: root box `rounded border border-hair p-3` with bg by variant (`bg-paper-sunk` for `"summary"`, `bg-plate` for `"compose"`). Inline star SVG `className="w-[15px] h-[15px] shrink-0"` `fill="var(--color-print-accent)"` (guard-clean). Optional bold title (`text-meta text-ink font-bold` or `font-semibold`). Body `text-body text-ink-soft` with emphasis rendered as `<strong className="text-ink font-semibold">` (caller passes `ReactNode` body so it can embed the bold spans). Actions row `flex gap-3 mt-2`: link-buttons are plain `<button>` with named-token classes only — accent: `className="text-data-strong text-print-accent hover:underline"` (or `text-sm font-semibold text-print-accent`), muted: `text-sm font-semibold text-ink-faint hover:underline`.

Props: `{ variant?: "summary"|"compose" = "summary"; title?: string; children: ReactNode /*body*/; actions?: Array<{ label: string; onClick?: () => void; muted?: boolean }>; className?: PLACEMENT-ONLY }`.

- [ ] **Step 1 — failing test**: root `data-testid="auto-group"` with `data-variant`; star renders (`[data-role="autogroup-star"]` or an `svg` — give the svg a `data-role`); title shows only when passed; body text renders; each action renders a `<button>` with its label and fires `onClick` on click; `variant="summary"` vs `"compose"` reflected in `data-variant`. No class-string assertions (the bg-by-variant is verified via `data-variant`, not the class).
- [ ] **Step 2** fail · **Step 3** implement (give the star `<svg data-role="autogroup-star" …>`; use the 5-point star path from the mockup) · **Step 4** stories: `Summary` (scoping copy, no actions), `Compose` (builder copy: title "Auto-grouped" + Confirm/Adjust actions) · **Step 5** gates · **Step 6** commit 3 paths.

> The link-buttons are intentionally NOT the `Button` primitive (those carry border/padding). Plain `<button>` + named-token text classes is contract-clean (no appearance literal). If review prefers a primitive, that's a follow-up `LinkButton` in `ui/`, not this task.

---

## Task 5 — `PhaseBlock` composite

**Files:** Create `src/print/components/PhaseBlock.tsx`, `.stories.tsx`, `test/print-components/PhaseBlock.test.tsx`. Depends on Task 2 (`PhaseLabel`).
**Mockup:** `docs/redesign-mockups/2026-05-29-focus-plot.html` — the assignment-cart assigned-phase row (`.as-block`/`.asb-*`).

One assigned phase in the cart. Three stacked lanes:
1. Top row `flex items-baseline gap-2`: serif name `<PhaseLabel phase={phase} className="text-display">` (snap: mockup 22px; nearest serif role is `text-display`/27 — document, or use `text-headline`/19 if 27 is too big; **pick `text-headline` (19) is too small, 22 sits between — use `text-display` and note the snap**), grow with `flex-1`; score value `<PhaseLabel phase={phase} className="text-data-strong font-bold">{score.toFixed(2)}</PhaseLabel>`; remove `<IconButton label={`Remove ${phase}`} dismiss tone="ghost" className="hover:text-print-accent" onClick={onRemove} />`.
2. Meta line `text-data text-ink-soft mt-1`: `{lattice} · {reflections} reflections` (caller supplies a `meta` string/node).
3. `<ScoreBar value={score} phase={phase} size="bar" className="mt-2" />`.

Props: `{ phase: string; score: number; meta: ReactNode; onRemove?: () => void; className?: PLACEMENT-ONLY }`.

- [ ] **Step 1 — failing test**: root `data-testid="phase-block"`; phase name text present (in a `[data-phase-label]`); score text present (`"0.92"`); meta text present; a remove button with `aria-label` `Remove <phase>` present and fires `onRemove` on click; `[data-score-bar]` present with `data-phase` == phase. No class assertions.
- [ ] **Step 2** fail · **Step 3** implement (import `PhaseLabel`, `ScoreBar`, `IconButton` from `../ui`) · **Step 4** stories (`Default`; `Coexistence`/high-score + low-score variants across phases) · **Step 5** gates · **Step 6** commit 3 paths.

> Decision: the remove × uses `IconButton dismiss tone="ghost"` + a `hover:text-print-accent` placement class to hit the mockup's accent-on-hover (named token, guard-clean). If review prefers a real `tone="accent-ghost"`, that's a `ui/` follow-up, not this task.

---

## Task 6 — `ReadingRow` composite

**Files:** Create `src/print/components/ReadingRow.tsx`, `.stories.tsx`, `test/print-components/ReadingRow.test.tsx`. Depends on Task 1 (`Swatch`) + Task 2 (`PhaseLabel`).
**Mockup:** `docs/redesign-mockups/2026-05-29-series-plot.html` — the reading panel `.rd-row`.

A phases-present summary line. Two stacked parts:
1. Top `flex items-center gap-[7px]`: `<Swatch phase={phase} shape="circle" />` + phase name `<PhaseLabel phase={phase} className="text-meta font-bold">` + a right-pushed span `<span className="ml-auto text-data text-ink-soft">{span}</span>` (the variable range, e.g. `1:0 → 1:0.75`).
2. Lattice line `<div className="text-data text-ink-faint ml-4">{lattice}</div>` (e.g. `a 205 → 195 Å`).
Root `flex flex-col gap-0.5`.

Props: `{ phase: string; span: ReactNode; lattice: ReactNode; className?: PLACEMENT-ONLY }`.

- [ ] **Step 1 — failing test**: root `data-testid="reading-row"`; `[data-swatch]` with `data-phase` == phase; phase name in `[data-phase-label]`; span text + lattice text render. No class assertions.
- [ ] **Step 2** fail · **Step 3** implement (import `Swatch`, `PhaseLabel`) · **Step 4** stories (`AllPhases` — a small stacked list across phases, with realistic spans/lattice from fixtures) · **Step 5** gates · **Step 6** commit 3 paths.

---

## Task 7 — `MemberRow` composite (series-builder `.trow` variant)

**Files:** Create `src/print/components/MemberRow.tsx`, `.stories.tsx`, `test/print-components/MemberRow.test.tsx`. Depends on Task 1 (`Swatch`).
**Mockup:** `docs/redesign-mockups/series-builder.html` — the trace-list row `.trow` (grip + square swatch + name/dose + phase chip).

> Scope: build ONLY the series-builder reorderable `.trow` variant. The series-plot `.member` variant (gradient coexistence swatch + inline colored phase names + lattice data) is a DIFFERENT row that belongs with the later MemberList/SeriesPlate work — do NOT build it here. Record this split in the ledger.

Row: root is a `group flex items-center gap-2 px-2 py-1.5 rounded` with `hover:bg-plate hover:border hover:border-hair` (resting transparent). Children in order:
1. `<GripHandle />` (root must be `group` so its color reveals on hover).
2. `<Swatch phase={dominantPhase} shape="square" />` (the 9px rounded color chip — dominant phase color).
3. Info column `flex-1 min-w-0`: name `<div className="text-meta font-semibold text-ink truncate">{name}</div>` + sub `<div className="text-data text-ink-faint truncate">{sub}</div>` (sample id / dose / variable value).
4. `<PhaseChip phase={dominantPhase} coexistWith={coexistWith} variant="tint" size="sm" className="shrink-0" />`.

Props: `{ name: string; sub?: ReactNode; phase: string; coexistWith?: string[]; className?: PLACEMENT-ONLY }`. (Drag wiring is the consumer's job; the grip is the affordance only.)

- [ ] **Step 1 — failing test**: root `data-testid="member-row"` and carries the `group` class only as needed (assert behavior, not class — instead assert the grip renders: `[aria-hidden]` grip or its `⋮⋮` text); `[data-swatch]` `data-phase` == phase; name + sub text render; `[data-testid="phase-chip"]` present with the phase; sub omitted → not rendered. No appearance-class assertions.
- [ ] **Step 2** fail · **Step 3** implement (import `GripHandle`, `Swatch`, `PhaseChip`) · **Step 4** stories (`Default`; `Coexistence` with `coexistWith`; a stacked `List` of 4 rows showing the hover/grip reveal) · **Step 5** gates · **Step 6** commit 3 paths.

---

## After the batch

- [ ] Full `npm test` (Vitest one-shot) green; `npm run lint:design` clean; `npx tsc --noEmit -p tsconfig.build.json` clean; `npm run build-storybook` exit 0.
- [ ] Flip ledger (`docs/greenfield-component-ledger.md`): Layer-1 `Swatch` stays ✅ (already), add the 2 new primitives to Layer-0 count (41 total); flip `PhaseBlock`, `FolioHeader`, `AutoGroup`, `MemberRow`, `ReadingRow` → ✅ in Layer-2; record the MemberRow builder-vs-plot split + the FolioHeader/PhaseBlock type-snap decisions in the registry. Update coverage summary.
- [ ] Manual fidelity pass: serve `storybook-static` and compare the 5 composites + 2 primitives against the cited mockup regions.
- [ ] Update memory (`project_greenfield_composite_layer.md` + `MEMORY.md`) with Batch 3.
- [ ] Do NOT merge/push — branch stays unmerged per standing constraint.
