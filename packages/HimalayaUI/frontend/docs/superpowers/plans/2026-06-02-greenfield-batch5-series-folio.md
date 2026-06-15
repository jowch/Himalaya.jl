# Greenfield Batch 5 — Series-folio vertical slice (Layer 3)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Fresh implementer per task + spec review + quality review. Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Typecheck with `npx tsc --noEmit -p tsconfig.build.json` (the default root tsc has unrelated errors). Commit ONLY the named files (never `git add -A`). Branch `worktree-greenfield-ui-rebuild` stays UNMERGED — do not push, do not offer PR/finish.

**Goal:** Build the first Layer-3 panel cluster — the Series folio (`SeriesCard` → `Gallery`) — so the greenfield rebuild has a full, demoable page-fragment in Storybook.

**Architecture:** Bottom-up. One refactor-on-contact primitive (`NoticePill`) unblocks the card; `SeriesCard` composes the already-built `CardFigure` (mini-waterfall) + `PhaseStrip` (phase strip + caption) + `NoticePill` + placement-only text/tokens; `Gallery` is a CSS-columns masonry of cards with an empty state.

**Tech Stack:** React + TypeScript, Tailwind (named-role + token classes only outside `ui/`), Vitest + RTL, Storybook 8 (`@storybook/react-vite`).

**Source of truth:** `docs/redesign-mockups/series-folio.html` (repo root). Read `.card`, `.card-fig`, `.card-body`, `.card-kick`, `.fig-n`, `.pill*`, `.card-title`, `.card-meta`, `.ps*` (handled by `PhaseStrip`), `.card-foot`, `.cf-prov`, `.cf-edit`, `.gallery`, `.empty`, and the `cardHTML`/`phaseStrip` JS.

---

## Verified API facts (checked against live source 2026-06-02)

- **`CardFigure`** (`src/print/waterfall/CardFigure.tsx`): `({ rows: WaterfallRow[]; width?: number /*340*/; className?: string })`. Root `data-testid="card-figure"`, `data-row-count`. Renders one `<div data-row-key data-phase>` per row, phase color via inline `style.color`. Reuses `TracePlot` per row (frozen: `interaction={false}`, `axes={false}`).
- **`PhaseStrip`** (`src/print/ui/PhaseStrip.tsx`): `({ segments: PhaseSegment[]; size?: "sm"|"md" /*md*/; emptyLabel?: string; orientation?: "horizontal"|"vertical" /*horizontal*/; className?: string })`. `PhaseSegment = { phase: string|null; coexistWith?: string[]|null; state?: "form_factor"|"null" }`. Horizontal renders the `.ps` cell row PLUS the `.ps-cap` caption (`data-testid="ps-cap"`) automatically — the caption derives "A → B" vs "throughout" vs `emptyLabel` from the COUNT OF DISTINCT indexed phases. Each cell `data-testid="ps-seg"`. **SeriesCard passes `segments` straight through — it does NOT re-implement the strip or caption.**
- **`WaterfallRow`** (`src/print/waterfall/waterfallModel.ts`): `{ key; label; trace; phase: string|null; state; anchors; bandHeight; yOffset }`. NOTE: a row carries only its DOMINANT `phase`, not coexisting phases — so the card's `PhaseStrip` segments are passed SEPARATELY (the page derives them), not reduced from `rows`.
- **Fixtures** (`src/print/waterfall/waterfall.fixtures.ts`): `FULL` (5 real rows: Ia3d → Im3m+Ia3d coexist → Im3m + dense cubic + sparse lamellar), `TRANSITION` (3 rows: Ia3d → coexist → Im3m), `MIXED_STATES`. Import these for `rows`.
- **`phaseColor(phase: string): string`** (`src/phases`) returns an `oklch(...)` string — a RUNTIME call, invisible to the design guard. Primitives apply it via inline `style`. JSDOM canonicalizes `oklch(0.570 …)`→`oklch(0.57 …)` on style write, so phase-color tests MUST round-trip through a reference element's `style` and compare read-backs (see `test/print-waterfall/CardFigure.test.tsx` `readBackColor` helper), NEVER `toBe(phaseColor(...))`.
- **Existing primitives for reference (do NOT use for the pill):** `Chip` (5 variants, full-pill border+bg — wrong look), `Badge` (flat mono count, no fill — wrong), `TagPill`. The folio pills are a distinct 10px tinted/dashed look → new primitive.
- **Design guard** (`scripts/check-design.mjs`): `ui/` is exempt (`isExcluded` → `components/ui/`, `print/ui/`, `print/plot/`, `print/detector/`, `print/comb/`, `print/export/`). `src/print/components/**` is NOT exempt — no inline appearance literals there (`text-[…]`, `rounded-[…]`, raw `oklch(`/`rgba(`/quoted-hex, side-stripes). `npm run lint:design` must stay green. Named roles (`text-base`, `text-caption`, `text-display`) and token classes (`bg-plate`, `border-hair`, `text-ink-faint`, `text-accent`) ARE allowed in composites.
- **Story idiom:** `import type { Meta, StoryObj } from "@storybook/react-vite";` + loose `const meta: Meta<typeof X> = {…}` (NOT `satisfies`). `@storybook/test` is NOT installed — use a plain `const noop = () => {};`, never `fn()`. Stories auto-glob `src/print/**/*.stories.tsx`.
- **Test idiom** (`test/AGENTS.md`): assert via `data-*` / roles / text / attributes — NEVER Tailwind/SVG class strings. Component test dirs mirror source: `test/print-ui/`, `test/print-components/`.

---

## Task A: `NoticePill` primitive (refactor-on-contact)

The card-kick pill in two tones: `new` (accent-tinted "+N new match", the recipe-found-samples signal) and `draft` (dashed faint "Draft"). Lives in `ui/` because it authors appearance (accent tint + dashed border) the placement-only card cannot.

**Files:**
- Create: `src/print/ui/NoticePill.tsx`
- Create: `src/print/ui/NoticePill.stories.tsx`
- Create: `test/print-ui/NoticePill.test.tsx`
- Modify: `src/print/ui/index.ts` (add the export — verify the barrel exists; if primitives are exported elsewhere, follow that pattern)

**Spec (derive exact values from `.pill`, `.pill-new`, `.pill-draft` in the mockup):**
- Props: `{ tone: "new" | "draft"; children: ReactNode; className?: string }`. (`new` = accent text on `color-mix(in oklab, var(--color-accent) 14%, transparent)` bg; `draft` = `text-ink-faint`, `bg-paper-sunk`, `1px dashed border-hair-strong`.)
- Geometry (mockup `.pill`): `font-size:10px; font-weight:700; letter-spacing:0.03em; padding:2px 8px; border-radius:999px`. Use the nearest named role / `rounded-full`; if a literal is unavoidable it is OK HERE (ui/ is guard-exempt), but prefer tokens/roles where they exist.
- Root `data-testid="notice-pill"`, `data-tone={tone}`. Accent color comes from `var(--color-accent)` / `--color-print-accent` tokens (NOT a raw literal) where possible.
- `className` is PLACEMENT-ONLY, appended last.

**TDD steps:**
- [ ] **A1** Write `test/print-ui/NoticePill.test.tsx`: (1) renders children text; (2) `data-tone` reflects the `tone` prop for both `"new"` and `"draft"`; (3) the `draft` variant is visually distinct from `new` (assert via `data-tone`, not class strings). Run → fail.
- [ ] **A2** Build `src/print/ui/NoticePill.tsx` from the mockup. Run test → pass.
- [ ] **A3** Add the `ui/index.ts` export (match the existing barrel pattern). Write `NoticePill.stories.tsx` (a `New` story `+2 new match` and a `Draft` story).
- [ ] **A4** `npx tsc --noEmit -p tsconfig.build.json` + `npm test -- print-ui/NoticePill` + `npm run lint:design` → all green.
- [ ] **A5** Commit ONLY the 4 named files. Message: `feat(print): NoticePill primitive (folio card-kick pill — new/draft tones)`.

---

## Task B: `SeriesCard` composite

One folio card: frozen mini-waterfall + body (kick row with fig-label + optional pill, serif title, meta line, phase strip, hairline footer with provenance + edit attribution). Draft variant gets a dashed border + dimmed figure.

**Files:**
- Create: `src/print/components/SeriesCard.tsx`
- Create: `src/print/components/SeriesCard.stories.tsx`
- Create: `test/print-components/SeriesCard.test.tsx`

**Spec (derive from `.card`, `.card-fig`, `.card-body`, `.card-kick`, `.fig-n`, `.card-title`, `.card-meta`, `.card-foot` + `cardHTML` JS):**
```ts
export interface SeriesCardProps {
  /** Mini-waterfall rows (low→high). Drives <CardFigure>. */
  rows: WaterfallRow[];
  /** Phase-strip segments (low→high). Passed straight to <PhaseStrip>;
   *  carries coexistence the rows don't, so it's separate from `rows`. */
  segments: PhaseSegment[];
  /** Kick label: "Fig. 1" for a saved series, "Recipe" for a draft. */
  figLabel: string;
  title: string;
  /** Meta line: count + ordering variable → "5 samples · by LL37 : lipid ratio". */
  sampleCount: number;
  variable: string;
  /** Footer left: beamtime string, or a cross-experiment node. */
  provenance: ReactNode;
  /** Footer right: "edited 2 days ago · JC" — pass the resolved edited label + author. */
  editedLabel: string;
  author: string;
  /** Optional kick pill. */
  notice?: { tone: "new"; count: number } | { tone: "draft" };
  /** Draft styling: dashed border + dimmed figure (mockup `.card.is-draft`). */
  draft?: boolean;
  onClick?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}
```
- Root `<article data-testid="series-card" data-draft={draft ? "true" : "false"}>`, clickable (`onClick`, `role`/keyboard per the existing composite idiom — check `CandidateRow`/`MemberRow` for the project's clickable-row convention; if cards are plain `onClick` divs in the mockup, keep parity but expose `data-testid`).
- Figure region wraps `<CardFigure rows={rows} />` with the mockup's `.card-fig` top-padding + sunk-line bottom border (placement/token only). Draft dims it.
- Kick row: `.fig-n` accent kicker text (`figLabel`) on the left; `notice` → `<NoticePill tone="new">+{count} new match</NoticePill>` or `<NoticePill tone="draft">Draft</NoticePill>` on the right (`justify-between`).
- Title: serif. Mockup is 19px — `text-headline` (19) is exact; use it (no snap needed). `data-testid` not required on the title; assert by text.
- Meta line: `<b>{sampleCount} samples</b> · by {variable}` — bold count uses `text-ink-soft`, middot + "by …" in `text-ink-faint`.
- `<PhaseStrip segments={segments} />` (renders cells + caption).
- Footer: hairline-top row, `justify-between`, `provenance` left + `edited {editedLabel} · {author}` right, all `text-caption`/`text-ink-faint`.

**TDD steps:**
- [ ] **B1** Write `test/print-components/SeriesCard.test.tsx`: (1) renders title text + `figLabel` + meta (`{n} samples`, `variable`); (2) embeds a `card-figure` with the right `data-row-count`; (3) embeds `ps-seg` cells (count === segments.length) and a `ps-cap`; (4) `notice={{tone:"new",count:2}}` → a `notice-pill` with `data-tone="new"` and "+2 new match"; (5) `draft` → root `data-draft="true"`; (6) `onClick` fires on click; (7) footer shows `editedLabel` + `author` + provenance. Assert via `data-*`/text/roles only. Run → fail.
- [ ] **B2** Build `src/print/components/SeriesCard.tsx`. Run test → pass.
- [ ] **B3** Write `SeriesCard.stories.tsx`: `Transition` (rows=`TRANSITION`, segments=`[{phase:"Ia3d"},{phase:"Im3m",coexistWith:["Ia3d"]},{phase:"Im3m"}]`, figLabel "Fig. 1", notice `{tone:"new",count:1}`), `Full` (rows=`FULL`, 5 segments matching, figLabel "Fig. 3"), `Draft` (rows=`TRANSITION` subset, `draft`, figLabel "Recipe", notice `{tone:"draft"}`), `CrossExperiment` (provenance = a node with the ⇄ + "April + July · q normalized"). Use `noop` for `onClick`.
- [ ] **B4** `npx tsc --noEmit -p tsconfig.build.json` + `npm test -- print-components/SeriesCard` + `npm run lint:design` → all green.
- [ ] **B5** Commit ONLY the 3 named files. Message: `feat(print): SeriesCard folio card (mini-waterfall + phase strip + provenance)`.

---

## Task C: `Gallery` composite

The masonry wall: a CSS-columns layout of `SeriesCard`×N, with an empty state when the (already-filtered) list is empty. The folio's search/sort/filter chrome is page-level (Layer 4) and is OUT OF SCOPE here — `Gallery` renders whatever cards it's given.

**Files:**
- Create: `src/print/components/Gallery.tsx`
- Create: `src/print/components/Gallery.stories.tsx`
- Create: `test/print-components/Gallery.test.tsx`

**Spec (derive from `.gallery` and `.empty` in the mockup):**
```ts
export interface GalleryProps {
  /** Already-filtered/sorted; Gallery only lays them out. */
  children: ReactNode;        // SeriesCard elements
  /** Shown (instead of children) when the list is empty. */
  empty?: ReactNode;
  className?: string;         // PLACEMENT-ONLY
}
```
- `.gallery`: `column-count:3; column-gap:22px` with responsive `@media` (2 cols ≤1180px, 1 col ≤720px). Tailwind columns utilities (`columns-1`/`columns-2`/`columns-3` + responsive prefixes) are NAMED utilities, allowed in a composite. Each card must `break-inside: avoid` + `margin-bottom` — the mockup puts this on `.card`; since `Gallery` wraps externally-passed `SeriesCard` children, apply the break-inside via a per-child wrapper `<div className="break-inside-avoid mb-...">` (mockup `.card { break-inside: avoid; margin-bottom: 22px }`).
- Root `data-testid="gallery"`. When `children` is empty (an empty array / no nodes) and `empty` is provided, render `<div data-testid="gallery-empty">{empty}</div>` instead.
- Decide "empty" by `React.Children.count(children) === 0`.

**TDD steps:**
- [ ] **C1** Write `test/print-components/Gallery.test.tsx`: (1) renders all passed children (e.g. 3 `series-card`s) inside `gallery`; (2) with zero children + an `empty` node → renders `gallery-empty` and NOT `gallery` cards; (3) with children present, does NOT render `gallery-empty`. Run → fail.
- [ ] **C2** Build `src/print/components/Gallery.tsx`. Run test → pass.
- [ ] **C3** Write `Gallery.stories.tsx`: `Wall` (3–5 `SeriesCard`s of varying row counts → shows the masonry's distinct heights), `Empty` (no children + an `EmptyState`/text empty node). Reuse the SeriesCard story fixtures.
- [ ] **C4** `npx tsc --noEmit -p tsconfig.build.json` + `npm test -- print-components/Gallery` + `npm run lint:design` → all green.
- [ ] **C5** Commit ONLY the 3 named files. Message: `feat(print): Gallery masonry wall of series cards + empty state`.

---

## Task D: Ledger + memory update

- [ ] **D1** In `docs/greenfield-component-ledger.md`: flip `SeriesCard` and `Gallery` (Layer 3 table, lines ~118–119) ⬜→✅ with file paths; add `NoticePill` to the Layer 0 primitive list (now 42) and to the watch-list resolution note; bump the coverage summary (Layer 3: 3 ✅ · ~14 ⬜; Layer 0: 42); rewrite the closing frontier paragraph to point at the next unbuilt panel cluster.
- [ ] **D2** Commit ONLY `docs/greenfield-component-ledger.md`. Message: `docs(print): ledger — Batch 5 Series-folio slice (NoticePill/SeriesCard/Gallery)`.
- [ ] **D3** (controller, after the batch) update the `project_greenfield_composite_layer` memory + `MEMORY.md` index line.

## Verification (after the batch)

`npm test` (full) green · `npx tsc --noEmit -p tsconfig.build.json` exit 0 · `npm run lint:design` clean · `npm run build-storybook` exit 0 · visual fidelity check of `components/SeriesCard` + `components/Gallery` stories against `series-folio.html`.
