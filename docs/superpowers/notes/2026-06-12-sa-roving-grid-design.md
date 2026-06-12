# SA-ROVING — Samples contact-sheet roving-tabindex grid (APG data grid)

> Design note for the impeccable polish loop's highest-complexity build. Owner decision
> 2026-06-12 (Jonathan): **full APG roving-tabindex grid**, not the skip-link minimum.
> Closes the last Samples P1 (`loop-backlog.md` → SA-ROVING). Companion to that row.

## Problem

`/samples` renders 139 rows × ~7 focusables/row = **~1105 page tab stops** (measured). A
keyboard-first user (persona Maya) Tabs through hundreds of stops to reach a lower row, the
footer, or the floating CullBar. Everything is *operable* (WCAG 2.1.1 / 2.4.3 pass — the
strict-AA remedy is a skip-link), but it is not *tractable*. The owner chose the flagship
remedy: turn the table into a real APG data grid with **one** tab stop and arrow-key
navigation.

## Current structure (verified in source 2026-06-12)

`SheetTable.tsx` → `role="table"` wrapper holding two rowgroups:
- **head** (`sheet-head`): one `role="row"` of `role="columnheader"` cells. Post-SA-SORT the
  4 data headers (Sample / Exposures / Kept / Status) are `ColumnSortButton`s carrying
  `aria-sort`; Tags is a static `Kicker` columnheader; the checkbox header is an empty cell.
- **body** (`sheet-rows`): one `SampleTableRow` per sample. Each row = a `role="row"` whose
  inner grid (`role="presentation"`) holds **6 columns**:

  | col | cell | focusables |
  |---|---|---|
  | 0 | checkbox | 1 (`Checkbox`) — only when `onCheck` set |
  | 1 | Sample (`SpecCell`) | 1 (name `<button>`, when `onOpenLoupe`) — **sticky** |
  | 2 | Exposures (`ThumbnailGallery`) | **N** thumbnails (click=select, dblclick=loupe) |
  | 3 | Kept (`KeptCell`) | **0** (static text) |
  | 4 | Tags (`TagList`) | **multi** (chips + remove ✕ + add) |
  | 5 | Status (`StatusCell`) | 1 (Focus-door `<button>`, when `onOpenFocus`) |

Load-bearing invariants that MUST survive: the shared `sampleTableCols()` grid template
(alignment), the sticky checkbox + Sample columns (`CHECKBOX_TRACK_PX`), the SA-SORT
`aria-sort` headers + sort-button behavior, the one shared `overflow-x-auto` scroller, and
every existing pointer flow (single-click select, ⇧-click range, dbl-click loupe, cull/keep).

## Keyboard model (APG "Data Grid")

**Grid = header row (index 0) + 139 data rows (1..N), 6 columns.** Exactly one cell is in the
tab order at a time (`tabindex=0`); all others `-1`. This removes the 4 sort-button tab stops
too — the header participates in the grid (the owner's "full APG" choice).

**Navigation mode** (default). Key handler on the grid container; `preventDefault` only on
keys it consumes so unrelated keys (⌘K, X/K cull) still bubble:

| key | effect |
|---|---|
| ↓ / ↑ | row ± 1 (same column), clamped |
| → / ← | col ± 1 (same row), clamped, **no wrap** |
| Home / End | first / last column in the row |
| Ctrl/⌘ + Home / End | grid corner (0,0) / (N, last) |
| PageDown / PageUp | row ± 10, clamped |

Landing on a cell focuses: the cell's **single widget** if it has exactly one (checkbox,
Sample button, Status door, a header sort button); otherwise the **gridcell element itself**
(Kept = 0 widgets; Exposures / Tags = multi → see interaction mode). `Enter`/`Space` on a
single-widget cell just activates the native control.

**Interaction mode** — only the two multi-widget cells (Exposures gallery, Tags). From the
focused gridcell, **Enter or F2** enters: focus moves to the first (or selected) inner widget;
**← / →** now rove the inner widgets (thumbnails / tag controls); **Tab** may also step through
them; **Escape** exits back to navigation mode (focus returns to the gridcell, arrow keys
resume cross-cell movement). While in interaction mode the grid's cross-cell arrow handler is
suppressed. The thumbnail gallery keeps its own `overflow-x-auto`; the active thumb is scrolled
into view on focus.

**Out of the grid:** one `Tab` from the active cell leaves the grid entirely and reaches the
footer / CullBar / ComposeBar (they are after the grid in DOM and remain normal tab stops).
⌘K and the X/K/Esc cull keys stay global (unaffected — the grid handler ignores them).

## Architecture

Three pieces, smallest-blast-radius first:

1. **`src/lib/grid/rovingGrid.ts` (pure, no DOM — fully unit-tested).**
   ```ts
   export interface GridCoord { row: number; col: number }
   export interface GridDims { rows: number; cols: number }
   /** Next coord for a navigation key, clamped (no wrap). null = not a nav key
    *  (let it bubble). pageSize default 10. ctrl = ⌘/Ctrl held (corner jumps). */
   export function nextGridCoord(
     c: GridCoord, key: string, dims: GridDims,
     opts?: { pageSize?: number; ctrl?: boolean },
   ): GridCoord | null
   ```
   This is the brain. ~40 lines, ~20 table-driven tests (each key, both clamp edges, ctrl
   corners, page jumps past the ends, non-nav key → null). Land it ALONE first (commit 1):
   zero behavior change, trivially green.

2. **Roving context + the grid wiring in `SheetTable` (commit 2).**
   - `SheetTable` becomes the grid owner: `role="grid"` (was `table`), holds the active
     `GridCoord` state, a cell registry (`Map<"r,c", HTMLElement>` populated via a
     `registerCell(row,col,el)` callback), an interaction-mode flag, and the keydown handler.
     A `useLayoutEffect` focuses the active cell's registered element **only after a
     user-driven coord change** — guard with a ref flag set by the key/click handler and
     cleared in the effect, exactly like `RepresentativeBox.tsx` so SSE / foreign re-renders
     never yank focus.
   - Provide a `RovingGridContext` (default = **inert**: `tabIndexFor` returns `undefined`,
     `register` is a no-op) so `SampleTableRow`'s other consumer (its story) and any future
     non-grid use render unchanged. The grid is "on" only when `SheetTable` supplies the
     provider (gate on `onSort != null`-style opt-in, or a new `roving` prop the page sets).
   - Cells read `tabIndex` + `ref`(register) from context by their `(row,col)`. Header cells
     register at row 0. `role="cell"` → `role="gridcell"` on body cells (grid requires
     gridcell); `columnheader` stays.
   - Row index: the **page** passes `rowIndex` to each `SampleTableRow` (it already `.map`s
     with an index available); the row passes `rowIndex` + a fixed col number to each cell's
     context lookup. Keeps SheetTable presentational; the page owns row identity.

3. **Interaction mode for the gallery + tags (same commit 2, or a tight follow-up).**
   `ThumbnailGallery` gains an optional roving mode (its thumbnails already are buttons); when
   the gallery cell is active+entered, one thumb is `tabindex=0` and ←/→ rove. `TagList`
   similarly. Both expose an `onExitInteraction` (Escape) that returns focus to the gridcell.

## Tests

- **Unit (`test/lib/grid/rovingGrid.test.ts`):** the reducer, table-driven (above).
- **Unit (`test/print-components/SheetTable.test.tsx` + `SampleTableRow.test.tsx`):** exactly
  one `tabindex=0` in the grid at rest; arrow keydown moves it (assert which element is
  focused / which coord is active via a `data-active-cell` attr — semantic, never class
  strings); `role="grid"` + `gridcell`; `aria-sort` headers still present; Enter enters
  interaction mode on the gallery cell, Escape exits.
- **e2e (`e2e/samples.spec.ts` or a new spec):** real-Chrome Tab reaches the grid in **one**
  stop; ArrowDown/Right move focus; Tab from inside the grid lands on the footer/CullBar, not
  the next cell. (jsdom focus is unreliable for the multi-listener path — the e2e is the
  honest gate, per the jsdom-dispatch false-green memo.)
- **Render-verify at :5182** in real Chrome: keyboard-drive the grid, confirm one tab stop,
  arrow nav, interaction mode in/out of the gallery, no SSE focus-yank, sort still works,
  sticky columns still freeze, 0 console errors.

## Risk & sequencing

This is the app's most load-bearing table and a brand-new interaction model. De-risk by:
commit 1 = pure reducer (proven in isolation) → commit 2 = wiring built on the proven reducer,
reviewed hard (frontend-reviewer + queue-reviewer for the SSE-focus interplay) and
render-verified before it lands. Preserve every pointer flow and the SA-SORT / sticky / scroller
invariants. If interaction mode proves too large for one clean commit, land nav-mode +
single/zero-widget cells first (gallery/tags reachable as a whole cell, thumbs still
individually tabbable as a transient) then interaction mode — but never ship an intermediate
where a gallery thumb is *less* reachable than today.
