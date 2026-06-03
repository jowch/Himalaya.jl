# Figure Export (`useFigureExport` + `ExportButton`) — design

**Date:** 2026-06-03 · **Branch:** `worktree-greenfield-ui-rebuild` (unmerged) · **Status:** approved design, pre-plan

## Context

Greenfield "The Print" rebuild. Layers 0–3 are complete; the deferred component
list named two "modals" with no mockup — **ExportSheet** and **ConflictModal**.
A 2026-06-03 brainstorm collapsed both:

- **ExportSheet is not a sheet.** Export is *content-WYSIWYG, style-transformed*
  via a split button. No preview, no options dialog. Buildable now (the engine
  already exists in `src/lib/figure-export/`).
- **ConflictModal is cancelled** — there is no conflict UI. Multiplayer becomes
  edit-tracking → undo/redo → versioning, designed during Layer 4.

This spec covers **the export control only**. The conflict decision is recorded
here as a non-goal and cross-referenced in `docs/redesign-notes.md`,
`docs/future-feature-ideas.md`, `docs/event-log.md`, and `CLAUDE.md`.

## Decisions

### Export model: content-WYSIWYG, style-transformed

The exported figure preserves **exactly the composed content** — which traces,
which peaks, the offsets, the order — but **sheds The Print's decorative in-app
styling** for clean scientific presentation (white ground, Arial, heavier
strokes). That clean idiom already exists as the engine's `CLEAN_SCIENTIFIC`
preset (`presets.ts`). There is **no export-time tuning, no preview, no options
dialog**. Copy short-circuits straight to the clipboard.

### Interaction: split button

Primary action **Copy** (one click → clean-render PNG → clipboard). A chevron
reveals a two-item download menu: **Download as PNG** / **Download as SVG**.
This mirrors the legacy `FigureExportControls` ergonomics, rationalized so the
dropdown is purely the two formats. It is **not** a modal.

### Non-goal: conflict resolution (cancelled, not deferred)

No conflict modal. The legacy `SeriesCommitConflictModal` / `ConflictModalShell`
pattern — server-vs-draft panels, Discard/Overwrite, driven by a `409` +
`expected_content_hash` — is **cancelled for the rebuild**. The product stance
is that commit conflicts should not surface as friction. Multiplayer is
last-write-wins, with **edit-tracking for undo/redo + versioning** as the
positive replacement, designed during Layer 4 (it touches the commit contract,
not just UI). This extends the existing deliberate LWW grain
(`exposures.selected`, `packages/HimalayaUI/src/AGENTS.md`). See
`docs/event-log.md` §"Conflict resolution" and `docs/future-feature-ideas.md`
§"Multi-user / review".

## Architecture

Two units, with the logic/presentation split the Print layer mandates
(components are presentational; side effects live in a hook or the page).

### `useFigureExport(spec, filenameStem, ariaContext)` — the logic hook

**File:** `src/print/hooks/useFigureExport.ts` (new `hooks/` dir under `print/`).

The **only** consumer of `src/lib/figure-export/*`. Encapsulates the engine
wiring once so every call site stays trivial. Verified engine API it composes:
`buildExportPng(spec, scale=2): Promise<Blob>`, `buildExportSvg(spec):
SVGSVGElement`, `canExportPng(): boolean`, `canCopyPngToClipboard(): boolean`,
`copyPngToClipboard(blob): Promise<void>`, `downloadBlob(blob, filename): void`,
`buildFilename(stem, "png"|"svg"): string`, `ExportSpec` (`types.ts`).

```ts
function useFigureExport(
  spec: () => ExportSpec,   // thunk — evaluated at click time for fresh state
  filenameStem: string,     // slug stem; hook appends date + extension
  ariaContext: string,      // e.g. "this figure", "the LL37 series"
): {
  onCopy: () => void;
  onDownloadPng: () => void;
  onDownloadSvg: () => void;
  copyDisabled: boolean;
  pngDisabled: boolean;
  pending: boolean;
}
```

- **onCopy** — `buildExportPng(spec())` → `copyPngToClipboard` → success toast;
  `catch` → error toast ("Couldn't copy figure. Try Download instead.").
- **onDownloadPng** — `buildExportPng(spec())` →
  `downloadBlob(blob, buildFilename(stem, "png"))`; `catch` → error toast.
- **onDownloadSvg** — `buildExportSvg(spec())` → `XMLSerializer` →
  `Blob([xml], { type: "image/svg+xml" })` →
  `downloadBlob(blob, buildFilename(stem, "svg"))`; `catch` → error toast.
- **pending** — true while a render is in flight; every handler bails while
  pending (the double-invocation guard).
- **copyDisabled** = `!canCopyPngToClipboard() || !canExportPng() || pending`.
- **pngDisabled** = `!canExportPng() || pending`. SVG is never capability-gated
  (it doesn't go through `OffscreenCanvas`); only `pending` blocks it.
- Toasts via `showToast` from `src/lib/toast` (the Print `ToastContainer`
  already wires it through `setToastImpl`).

### `ExportButton` — the presentational split button

**File:** `src/print/components/ExportButton.tsx`.

Composes `Button` + `Menu` primitives. Owns **only** its menu-open UI state and
Escape / outside-pointerdown dismissal — the exact precedent `Field` sets. Every
side effect is a prop.

```
┌───────────────┬───┐
│     Copy      │ ▾ │   ▾ opens Menu: "Download as PNG" / "Download as SVG"
└───────────────┴───┘
```

**Props:** `{ onCopy, onDownloadPng, onDownloadSvg, copyDisabled?, pngDisabled?,
pending?, disabled?, ariaContext, className? }`.

- Primary **Copy** `Button` → `onCopy`; disabled when `disabled || copyDisabled`.
- Chevron `IconButton` → toggles the `Menu`; carries `aria-haspopup="menu"` /
  `aria-expanded`.
- `Menu` options:
  `[{ value: "png", label: "Download as PNG", disabled: disabled || pngDisabled },
    { value: "svg", label: "Download as SVG", disabled: disabled || pending }]`,
  `onSelect` routes to `onDownloadPng` / `onDownloadSvg` then closes.
- Split-group border + divider use token classes (`border-hair`, `w-px bg-hair`,
  `rounded`) — design-guard-clean, **no new primitive**.
- `data-testid`: `export-button` (root), `export-copy`, `export-menu-trigger`,
  `export-menu` (`Menu` `aria-label="Download formats"`).
- `ariaContext` fills aria-labels: `Copy {ariaContext} to clipboard`,
  `Download {ariaContext} as PNG/SVG`.

A page wires the two together:

```tsx
const ex = useFigureExport(spec, stem, "the LL37 series");
<ExportButton {...ex} ariaContext="the LL37 series" disabled={dataNotReady} />
```

### Styling responsibility (why this build is tiny)

The "scientific presentation" render is the **engine + per-figure adapter's**
job: the `CLEAN_SCIENTIFIC` preset and the `ExportSpec` the adapter builds,
realized at **Layer-4 wire time**. `useFigureExport` and `ExportButton` are
**style-agnostic** — they render whatever spec is handed them. This build adds
**no** engine or adapter changes.

## Data flow

page `spec` thunk → `useFigureExport` → `figure-export` engine →
clipboard / download; toasts via `lib/toast` → `ToastContainer`.

## Error & capability handling

- Clipboard unsupported (no HTTPS / no Clipboard API) → Copy disabled, `title`
  explains.
- PNG renderer unsupported → Copy **and** Download-PNG disabled; SVG stays
  available.
- Render or copy failure → error toast (distinct copy vs download messages).
- Double invocation guarded by `pending` — handlers bail while a render is in
  flight.

## Testing

- **Hook test** (`test/print-components/useFigureExport.test.tsx`): mock
  `figure-export/*` + `lib/toast`; assert each handler's call chain, the
  success/error toast on copy + download, and `copyDisabled` / `pngDisabled`
  across capability combinations (clipboard on/off, PNG on/off, pending).
- **Component test** (`test/print-components/ExportButton.test.tsx`): split
  structure, chevron opens `Menu` (two items), primary click → `onCopy`,
  PNG/SVG select → respective handler + menu closes, disabled states (copy and
  PNG), Escape + outside-pointerdown close the menu. Class-agnostic
  (`data-testid` / role / text only).
- **Stories** (`ExportButton.stories.tsx`): spy handlers — default, menu-open,
  `copyDisabled`, `pending`, fully `disabled`. No engine, no canvas.
- **Gate:** `npm test` (print suites) green; `npm run lint:design` clean;
  `npx tsc --noEmit -p tsconfig.build.json` clean; `npm run build-storybook`
  exit 0.

## Out of scope (this build)

- **Call-site placement** (focus/series plate heads, the builder rail-foot
  "Copy as PNG", loupe) — Layer-4 page wiring.
- **The clean-render adapters** per figure — they exist in `lib/figure-export`;
  wired at Layer 4.
- **Conflict resolution** — cancelled (see Decisions); the undo/redo/versioning
  replacement is its own Layer-4 spec.

## Ledger

The `ExportSheet` entry is replaced by **`ExportButton` + `useFigureExport`**.
`ConflictModal` is marked **⛔ out of scope**. Recorded in the
*Out-of-scope & deferred registry* of
`packages/HimalayaUI/frontend/docs/greenfield-component-ledger.md`.
