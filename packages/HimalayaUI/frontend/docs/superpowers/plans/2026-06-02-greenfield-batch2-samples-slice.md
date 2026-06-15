# Greenfield Batch 2 — Samples-page vertical slice + Sparkline renderer

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development (fresh implementer per task + two-stage review). Steps use `- [ ]`.
> Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Commit ONLY the files each task names — never `git add -A` (unrelated pre-existing modified files exist).
> Typecheck with `npx tsc --noEmit -p tsconfig.build.json` (the default root tsc has unrelated errors).
> Branch `worktree-greenfield-ui-rebuild` stays UNMERGED / UNPUSHED. Do not offer merge/finish/PR.

**Goal:** Build the components that deliver two of the three user-named gaps — the **ThumbnailGallery** (used on 3 surfaces) and the **Samples-page table rows** — bottom-up, plus the high-reuse **`Sparkline`** Layer-1 micro-renderer.

**Architecture:** Bottom-up within the slice. Layer-1 `Sparkline` (independent) → Layer-2 `Thumbnail` → Layer-2 `ThumbnailGallery` (deps Thumbnail) → Layer-2 sample-table cell composites (independent) → Layer-3 `SampleTableRow` (composes the gallery + cells + `CheckCircle` + `TagList`). Every component is derived from its mockup region, NOT ported from legacy `src/components/`.

**Tech stack:** React + TS, Vitest + RTL, Storybook (auto-globs `src/print/**/*.stories.tsx`), Tailwind v4 with `@theme` tokens, d3-flavored projection helpers in `src/print/plot/`.

---

## Conventions (the per-component recipe — same as Batch 1)

- **Location:** `src/print/<dir>/<Name>.tsx` + colocated `<Name>.stories.tsx`. Unit tests in `test/print-<dir>/<Name>.test.tsx`.
  - Renderers (`Sparkline`) live in **`src/print/plot/`** (design-guard-EXEMPT) → may own pixels/SVG literals. Test dir `test/print-plot/`.
  - Composites (`Thumbnail`, `ThumbnailGallery`, cells, `SampleTableRow`) live in **`src/print/components/`** (NOT exempt) → **placement + named-token classes only**. Test dir `test/print-components/`.
- **Closed-look, placement-only** for composites: compose primitives + layout utilities (flex/grid/gap/margin/min-width) + named type roles (`text-data`, `text-meta`, `text-xs`, `text-caption`…) + token utilities (`bg-paper-sunk`, `border-hair`, `text-ink-faint`…). **No** inline appearance literals (`text-[…]`, `rounded-[…]`, raw `oklch(...)`/hex, side-stripes). `npm run lint:design` enforces this and MUST stay green.
- **Refactor-on-contact:** if a composite needs a size/treatment no role or primitive exposes, extend the role (`styles.css`) or primitive (`src/print/ui/`, exempt) as a reviewed sub-step — never scatter a one-off literal in the composite.
- **Token mapping (mockup → design system):** `--paper`→`bg-paper`, `--paper-sunk`→`bg-paper-sunk`, `--plate`→`bg-plate`, `--ink`→`text-ink`, `--ink-soft`→`text-ink-soft`, `--ink-faint`→`text-ink-faint`, `--hair`→`border-hair`, `--hair-strong`→`border-hair-strong`, `--accent`→`bg-accent`/`text-accent`, `--frame-edge`→`bg-frame-edge`. Phase colors → `phaseColor()` from `src/phases.ts` (renderer) or the `PhaseChip` primitive (composite).
- **Tests assert structure/behaviour via `data-*` / roles / text — never class strings.**
- **Per component (TDD):** write failing test → run, see it fail → build from primitives/mockup → write story against `src/print/fixtures/` → `npm test -- <name>` + `npm run lint:design` + `npx tsc --noEmit -p tsconfig.build.json` all green → commit only the named files.
- **Ledger:** flip the component ⬜→✅ in `docs/greenfield-component-ledger.md` (with file path) in the SAME commit that lands it.

Fixtures available (`src/print/fixtures/`): `realTraces.ts` (`realTraces[exposureId]` → `Trace` for 37/65/66/67/93), `realMembers` (`realSeriesMembers.ts`, real phases/lattices), `thumbs/<id>.png` (5 detector thumbnails), `traces/<id>.json`. Extend with small synthetic data as needed for table-row states.

---

## Task 1: `Sparkline` (Layer-1 micro-renderer)

The tiny 76×28 trace thumbnail in `series-scoping.html` sample rows and candidate rows. A dedicated micro-renderer (NOT a `TracePlot` reuse — no axes, no peaks, no interaction; just a filled area + stroke + baseline).

**Files:**
- Create: `src/print/plot/Sparkline.tsx`
- Create: `src/print/plot/Sparkline.stories.tsx`
- Test: `test/print-plot/Sparkline.test.tsx`

**Spec (from `series-scoping.html` `.spark` + the `sparkline(comps)` JS):**
- Container `<span data-testid="sparkline">`: `76px × 28px`, `bg-paper-sunk`, `border border-hair`, `rounded-sm`, `overflow-hidden`, `flex-shrink-0`, `inline-block`. (This dir is exempt, so an exact-px box via inline style or `w-[76px] h-[28px]` is allowed — prefer inline `style={{ width: 76, height: 28 }}` for the fixed figure dimension.)
- Inner `<svg viewBox="0 0 76 28">`, `width=100% height=100%`, `aria-hidden`, `display:block`.
- Geometry: left/right pad `padX = 4`, bottom pad `padB = 4` (baseline at `y = 24`), top headroom so the peak amplitude spans ~17px (`y` range `[24, 5]`).
- x is **log-q**: use `makeAxis([qMin, qMax], [4, 72], "log")` from `./projection` (import `makeAxis`). Map each sample's q→x via `axis.to(q)`.
- Marks, in z-order: (1) baseline `<line x1=4 y1=24 x2=72 y2=24 stroke="var(--color-hair)" stroke-width=1 data-role="spark-baseline">`; (2) filled area `<path data-role="spark-area" fill={color} fill-opacity={0.1}>` closing down to the baseline; (3) trace stroke `<path data-role="spark-line" fill="none" stroke={color} stroke-width={1.4} stroke-linejoin="round">`.
- **Color:** `color = phase ? phaseColor(phase) : "var(--color-ink-faint)"` (the mockup's `--unindexed` ochre-gray maps to `text-ink-faint`'s value; use the token). Import `phaseColor` from `../../phases`.

**Props interface:**
```ts
export interface SparklineProps {
  /** Measured trace (q ascending, I parallel). */
  trace: { q: number[]; I: number[] };
  /** Dominant phase for hue; null/undefined → unindexed (ink-faint). */
  phase?: string | null;
  /** PLACEMENT ONLY. */
  className?: string;
}
```
Build the path by sampling the provided trace directly (map each `(q,I)` to `(axis.to(q), yScale(I))`); a linear y-scale over `[0, max(I)]` → pixel `[24, 5]`. Guard empty/degenerate traces (≤1 point or flat) by rendering just the baseline.

**TDD steps:**
- [ ] Write `test/print-plot/Sparkline.test.tsx`: (a) renders `[data-testid="sparkline"]`; (b) renders one `[data-role="spark-line"]` path with a non-empty `d`; (c) renders `[data-role="spark-baseline"]`; (d) phase prop drives stroke — with `phase="Lamellar"` the line stroke `getAttribute("stroke")` `includes("oklch")` (a real phase color), and with `phase={null}` the stroke resolves to the ink-faint token (`includes("--color-ink-faint")` or the oklch value — assert it differs from the Lamellar stroke); (e) a degenerate trace (`{q:[0.1], I:[5]}`) still renders the baseline and does not throw.
- [ ] Run `npm test -- Sparkline` → fails (no file).
- [ ] Build `Sparkline.tsx` per spec.
- [ ] Write `Sparkline.stories.tsx`: stories `Indexed` (a `realTraces[65]` slice + `phase="Ia3d"`), `Unindexed` (`phase={null}`), `Lamellar`, and a `Row` story placing several in a flex row to show the filmstrip feel.
- [ ] `npm test -- Sparkline` + `npm run lint:design` + `npx tsc --noEmit -p tsconfig.build.json` → green.
- [ ] Flip `Sparkline` ⬜→✅ in `docs/greenfield-component-ledger.md` (file `plot/Sparkline.tsx`).
- [ ] Commit: `git add src/print/plot/Sparkline.tsx src/print/plot/Sparkline.stories.tsx test/print-plot/Sparkline.test.tsx docs/greenfield-component-ledger.md`

---

## Task 2: `Thumbnail` (Layer-2 composite)

Mini detector-image thumbnail with frame-number, representative dot, and reject overlay. Used in the sample-table exposures cell, the focus DetectorPanel filmstrip, and the loupe strip (**3 sites**).

**Files:**
- Create: `src/print/components/Thumbnail.tsx`
- Create: `src/print/components/Thumbnail.stories.tsx`
- Test: `test/print-components/Thumbnail.test.tsx`

**Spec (from `sample-table.html` `.thumb`):**
- Root `<button data-testid="thumbnail" data-state={…}>`: `relative`, fixed square (default 62px; loupe variant 70px), `rounded-sm` (mockup uses 3px — acceptable as `rounded-sm`=5px is the single control radius; if 3px fidelity matters, that's a refactor-on-contact note, not a one-off literal), `overflow-hidden`, `bg-frame-edge`, `border border-frame-edge`, `flex-shrink-0`, `p-0`, `cursor-pointer`. Hover ring + selected ring via tokens (see states).
- Inner detector image: compose the existing **`DetectorImage`** primitive (from `src/print/detector/`) filling the button (`w-full h-full`, `display:block`). Pass it the thumbnail source the consumer supplies. (Check `DetectorImage`'s props — it renders the warm-LUT detector SVG; if it doesn't accept a plain raster thumb source, use its image/href prop or fall back to an `<img>`/`<svg image>` element inside the button. Resolve which at build time by reading `src/print/detector/DetectorImage.tsx`; prefer composing the primitive.)
- Frame-number `<span data-role="thumb-fno">`: bottom-left (`absolute left-1 bottom-0.5`), mono, ~`text-[8.5px]`-class — there is NO role this small; use `text-xs` (10.5) as the nearest role and note the divergence, OR add an `xxs` micro-role in `styles.css` as a reviewed refactor-on-contact if 8.5px fidelity is required. Color: the mockup uses a light mono over the dark frame — reuse `text-frame-tag` (the existing token for mono labels over frame-edge). `opacity-80`.
- Representative marker `<span data-role="thumb-rep">`: top-right (`absolute top-1 right-1`), shown only when `representative`. Compose **`Dot` tone="accent"** with a plate ring. Mockup is 9px + 1.5px plate border; `Dot size="sm"` is 8px — close; the plate ring is `ring`/`border` via token. If exact 9px+ring needed, extend `Dot` with a `ring` affordance (refactor-on-contact) rather than inlining.
- Reject overlay: compose **`RejectOverlay`** primitive (absolute inset-0), shown only when `rejected`; ALSO dim the image (`opacity-[0.32]` is appearance — instead apply a token/util: use `opacity-30` nearest util, acceptable as opacity is layout-ish; confirm lint passes) and dim the frame-no.
- **States** (drive via `data-state` and boolean props, assert via data-attrs): `normal`, `representative` (rep dot shown), `rejected` (RejectOverlay shown + image dimmed), `selected` (accent ring: `ring-2 ring-accent`-style via token util). Hover (non-selected): `hover:` ring in `hair-strong`. These are togglable independently (a thumb can be representative AND selected).

**Props interface:**
```ts
export interface ThumbnailProps {
  /** Detector thumbnail source (raster href or the data DetectorImage needs). */
  src: string;
  /** Frame number label, e.g. "65". */
  frameNo?: string | number;
  representative?: boolean;
  rejected?: boolean;
  selected?: boolean;
  /** Loupe variant → 70px instead of 62px. */
  size?: "sheet" | "loupe";
  onClick?: () => void;
  title?: string;
  /** PLACEMENT ONLY. */
  className?: string;
}
```

**TDD steps:**
- [ ] Write `test/print-components/Thumbnail.test.tsx`: (a) renders `[data-testid="thumbnail"]` as a `<button>`; (b) shows `[data-role="thumb-fno"]` with the frameNo text; (c) `representative` → `[data-role="thumb-rep"]` present, absent otherwise; (d) `rejected` → `[data-testid="reject-overlay"]` present, absent otherwise; (e) `selected` reflected via `data-state` containing `selected` (or a dedicated `data-selected`); (f) click fires `onClick`; (g) `size="loupe"` reflected via a `data-size` attr.
- [ ] Run `npm test -- Thumbnail` → fails.
- [ ] Read `src/print/detector/DetectorImage.tsx` to settle how to pass the thumb source; build `Thumbnail.tsx` composing `DetectorImage` + `Dot` + `RejectOverlay` + text.
- [ ] Write `Thumbnail.stories.tsx` using `fixtures/thumbs/65.png` etc.: stories `Normal`, `Representative`, `Rejected`, `Selected`, `RepresentativeAndSelected`, `Loupe`.
- [ ] `npm test -- Thumbnail` + `npm run lint:design` + tsc → green. (If lint flags a needed appearance, do the refactor-on-contact in `ui/`/`styles.css` and commit that file too.)
- [ ] Flip `Thumbnail` ⬜→✅ in the ledger (file `components/Thumbnail.tsx`).
- [ ] Commit only the named files (+ any exempt-layer primitive/role file touched by a refactor-on-contact).

---

## Task 3: `ThumbnailGallery` (Layer-2 composite)

The filmstrip wrapper holding N `Thumbnail`s — the exposures cell in a table row AND the loupe strip.

**Files:**
- Create: `src/print/components/ThumbnailGallery.tsx`
- Create: `src/print/components/ThumbnailGallery.stories.tsx`
- Test: `test/print-components/ThumbnailGallery.test.tsx`

**Spec (from `sample-table.html` `.strip` and `.loupe-strip`):**
- `<div data-testid="thumbnail-gallery">`: flex row. Sheet variant: `gap-[7px]` → use `gap-2` (8px) nearest util, or keep 7px as a placement literal (gaps are placement, allowed). `flex-nowrap overflow-x-auto`. Loupe variant: `gap-2 justify-center`, `mt-3`.
- Renders one `Thumbnail` per exposure; passes through `representative`/`rejected`/`selected` per item; forwards `onSelect(exposureId)`.

**Props interface:**
```ts
export interface GalleryExposure {
  id: number;
  src: string;
  frameNo?: string | number;
  representative?: boolean;
  rejected?: boolean;
}
export interface ThumbnailGalleryProps {
  exposures: GalleryExposure[];
  /** Currently-selected exposure id (drives Thumbnail `selected`). */
  selectedId?: number;
  onSelect?: (id: number) => void;
  variant?: "sheet" | "loupe";
  className?: string;
}
```

**TDD steps:**
- [ ] Write `test/print-components/ThumbnailGallery.test.tsx`: (a) renders `[data-testid="thumbnail-gallery"]`; (b) one `[data-testid="thumbnail"]` per exposure; (c) clicking a thumb fires `onSelect` with its id; (d) `selectedId` marks exactly that thumb selected; (e) `variant="loupe"` forwards `size="loupe"` to children (assert via child `data-size`).
- [ ] Run → fails.
- [ ] Build `ThumbnailGallery.tsx` composing `Thumbnail`.
- [ ] Write `ThumbnailGallery.stories.tsx`: `Sheet` (5 fixture exposures, one representative, one rejected), `Loupe`, `Selected`.
- [ ] Tests + lint:design + tsc green.
- [ ] Flip `ThumbnailGallery` ⬜→✅ in the ledger.
- [ ] Commit named files.

---

## Task 4: Sample-table cell composites — `SpecCell`, `KeptCell`, `StatusCell`

Three small, cohesive Layer-2 composites for the sample-table row. (`TagsCell` is NOT a new component — it is the existing `TagList` primitive used directly in the row.) One implementer builds all three; reviewed together.

**Files:**
- Create: `src/print/components/SpecCell.tsx`, `src/print/components/KeptCell.tsx`, `src/print/components/StatusCell.tsx`
- Create: `src/print/components/SampleCells.stories.tsx` (one stories file covering all three)
- Test: `test/print-components/SampleCells.test.tsx` (one test file covering all three)

**`SpecCell` spec (`sample-table.html` `.spec`):** flex row, `gap-[11px]`→`gap-2.5`(10px) or keep 11px placement. Compose **`CheckCircle`** (`checked={screened}`, `label={screened ? "Screened" : "Not screened"}`, `self-start mt-0.5`) + a text column: `<span data-role="spec-name">` (sample name; nearest role `text-title` is 15px/600 — mockup is 13.5px/600; use `text-title` and note divergence OR `text-meta`/`text-body`; pick `text-body` + `font-semibold` if 13px reads closer — implementer's fidelity call against the mockup) and `<span data-role="spec-id">` mono id (`text-data text-ink-faint`, `text-xs`, `mt-[3px]`).
```ts
export interface SpecCellProps { name: string; sampleId: string; screened?: boolean; className?: string; }
```

**`KeptCell` spec (`.kept`):** mono. `<span data-role="kept-count">{kept}</span>` (large mono, `text-data` + size — mockup 16px; nearest role is small, so use `text-data` and bump weight, or add a `text-data-lg` role as refactor-on-contact) + `<span data-role="kept-total"> / {total}</span>` (`text-ink-faint`, smaller). When `dropped > 0`, `<span data-role="kept-dropped">{dropped} dropped</span>` block in `text-accent`, sans, `text-xs font-semibold`, `mt-[3px]`.
```ts
export interface KeptCellProps { kept: number; total: number; dropped?: number; className?: string; }
```

**`StatusCell` spec (`.status`):** if `phase` set → render **`PhaseChip`** `phase={phase}` (its tint variant matches the mockup's `.pchip` faint-tint exactly). Else → unset placeholder `<span data-role="status-unset">`: compose **`Dot` tone="muted" size="xs"`** (the hollow 6px ring, `mr-1.5`) + text "Not indexed" (`text-ink-faint`, `text-xs`/`text-caption`).
```ts
export interface StatusCellProps { phase?: string | null; className?: string; }
```

**TDD steps:**
- [ ] Write `test/print-components/SampleCells.test.tsx`:
  - SpecCell: renders name + id text; `screened` → `[data-testid="check-circle"]` `data-checked="true"`, unset → no `data-checked`.
  - KeptCell: renders `kept`/`total`; `dropped>0` → `[data-role="kept-dropped"]` present with count, `dropped=0`/undefined → absent.
  - StatusCell: `phase="Pn3m"` → `[data-testid]` of PhaseChip present with "Pn3m" text; `phase={null}` → `[data-role="status-unset"]` present with "Not indexed".
- [ ] Run → fails.
- [ ] Build the three cell files.
- [ ] Write `SampleCells.stories.tsx`: a story per cell covering each state (SpecCell screened/unscreened; KeptCell with & without dropped; StatusCell indexed/unset).
- [ ] Tests + lint:design + tsc green.
- [ ] Flip `SpecCell`/`KeptCell`/`StatusCell` (the `KeptCell/TagsCell/StatusCell/SpecCell` ledger row) ⬜→✅, noting `TagsCell = TagList (primitive, no new file)`.
- [ ] Commit named files.

---

## Task 5: `SampleTableRow` (Layer-3 — the user-named "table rows")

One full contact-sheet row. Composes everything above. This is the deliverable that makes the Samples page legible.

**Files:**
- Create: `src/print/components/SampleTableRow.tsx`
- Create: `src/print/components/SampleTableRow.stories.tsx`
- Test: `test/print-components/SampleTableRow.test.tsx`

**Spec (`sample-table.html` `.srow` + `.COLS`):**
- Root `<div data-testid="sample-table-row" data-screened={…}>`. Bottom hairline `border-b border-hair` (last row: caller drops it). Hover bg `hover:bg-paper-sunk` (mockup uses a subtle warm tint — `bg-paper-sunk` is the nearest token; confirm acceptable). Unscreened rows: `data-screened="false"` + a faint resting tint (`bg-paper-sunk` or a dedicated token — implementer's call; must be token-based).
- Grid `.COLS`: `grid grid-cols-[244px_minmax(360px,1fr)_78px_168px_150px]` (Tailwind arbitrary grid-template is a LAYOUT utility, not appearance — allowed). Each cell `flex items-center px-4 py-[13px] min-h-[92px]` (padding/min-height are placement — allowed).
- Column order + content:
  1. **Sample** → `SpecCell` (name, sampleId, screened).
  2. **Exposures** → `ThumbnailGallery variant="sheet"` (exposures, selectedId, onSelect).
  3. **Kept** → `KeptCell` (kept, total, dropped).
  4. **Tags** → `TagList` (the existing primitive; pass `tags`, `editable`, `onAdd`/`onRemove`). The hover-reveal of the add invite is already built into `TagList` (`group/tags`) — the row just needs the `group/tags` hover context, which `TagList` owns.
  5. **Status** → `StatusCell` (phase).

**Props interface:**
```ts
import type { Tag } from "../ui/tag";
export interface SampleTableRowProps {
  name: string;
  sampleId: string;
  screened?: boolean;
  exposures: GalleryExposure[];        // from ThumbnailGallery
  selectedExposureId?: number;
  onSelectExposure?: (id: number) => void;
  kept: number;
  total: number;
  dropped?: number;
  tags: Tag[];
  onAddTag?: (t: Tag) => void;
  onRemoveTag?: (t: Tag) => void;
  phase?: string | null;
  className?: string;
}
```

**TDD steps:**
- [ ] Write `test/print-components/SampleTableRow.test.tsx`: (a) renders `[data-testid="sample-table-row"]`; (b) contains the SpecCell name + id; (c) contains a `[data-testid="thumbnail-gallery"]` with one thumb per exposure; (d) contains the KeptCell counts; (e) `phase` set → PhaseChip, unset → status-unset placeholder; (f) `data-screened` reflects the `screened` prop; (g) selecting a thumb bubbles `onSelectExposure`.
- [ ] Run → fails.
- [ ] Build `SampleTableRow.tsx`.
- [ ] Write `SampleTableRow.stories.tsx`: `Indexed` (full row, phase set, some dropped, real fixture thumbs), `Unindexed` (phase null, unscreened), `AllKept` (no dropped), and a `Stack` story rendering 3–4 rows under a faux header to show alignment against the `grid-cols` track.
- [ ] Tests + lint:design + tsc green.
- [ ] Flip `SampleTableRow` ⬜→✅ in the ledger (Layer-3 table).
- [ ] Commit named files.

---

## Batch verification (after all 5 tasks)

- [ ] `npm test -- print-components print-plot` → all green.
- [ ] `npm run lint:design` → clean (proves placement-only held across the batch).
- [ ] `npx tsc --noEmit -p tsconfig.build.json` → clean.
- [ ] `npm run build` → exit 0.
- [ ] `npm run build-storybook` → exit 0.
- [ ] Manual fidelity (optional): `npm run storybook`, compare the new `components/*` + `plot/Sparkline` stories against `sample-table.html` / `series-scoping.html`.
- [ ] Ledger coverage summary updated (tier-1 composites count, Layer-1 renderer count, the one Layer-3 `SampleTableRow` ✅).

## Out of scope for this batch (roadmapped, not built here)

`CullBar` (floating batch-reject bar — its own small task), `SheetTable` (the grid wrapper that stacks `SampleTableRow`s + header + `CullBar` — Layer-3, next batch), `BigFrame`/`LoupeSidePanel` (loupe panels — consume `ThumbnailGallery` but are a separate surface), and the Series tier-1 composites (`PhaseBlock`/`FolioHeader`/`AutoGroup`/`MemberRow`/`ReadingRow`). The `CustomIndexModal` remains gated on the `CustomPreview` Layer-1 renderer.
