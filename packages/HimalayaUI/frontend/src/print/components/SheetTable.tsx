import React from "react";
import type { ReactNode } from "react";
import { Card, Kicker, ColumnSortButton } from "../ui";
import type { ColumnSortDir } from "../ui";
import { sampleTableCols, CHECKBOX_TRACK_PX } from "./SampleTableRow";
import {
  useRovingGrid,
  RovingGridProvider,
  type RovingGridValue,
} from "../../lib/grid/useRovingGrid";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** The page-owned sort the header reflects. `key` is one of the four sortable
 *  column keys (or null = ingest order); `dir` only matters when key is set. */
export interface SheetSort {
  key: string | null;
  dir: ColumnSortDir;
}

export interface SheetTableProps {
  /** SampleTableRow elements, already mapped/sorted by the page. */
  children?: ReactNode;
  /** Rendered (instead of the rows region) when there are zero children. */
  empty?: ReactNode;
  /** When true, prepends a blank checkbox header cell + a grid track
   *  CHECKBOX_TRACK_PX wide so the header stays aligned with rows that render
   *  their own checkbox column. */
  checkboxColumn?: boolean;
  /** Current sort, reflected into the columnheaders' `aria-sort` + the active
   *  column's caret. Only consulted when `onSort` is also provided. */
  sort?: SheetSort;
  /** When provided, the four data headers (Sample / Exposures / Kept / Status)
   *  render as sortable buttons that call this with the clicked column key; the
   *  page owns the toggle cycle. When ABSENT, the static kicker headers render
   *  (other consumers stay byte-compatible). Tags is never sortable. */
  onSort?: (key: SortColumnKey) => void;
  /** When true, the table becomes a WAI-ARIA roving-tabindex DATA GRID
   *  (`role="grid"`, body cells `role="gridcell"`, ONE tab stop, arrow-key
   *  navigation). The page also passes `rowIndex={i+1}` to each SampleTableRow
   *  so cells know their (row,col). When absent the static (non-roving)
   *  `role="table"` path renders byte-identically — other consumers unaffected. */
  roving?: boolean;
  /** Number of data rows — required when `roving` is on so the grid dims are
   *  correct (header row 0 + dataRows). Defaults to the slotted child count. */
  dataRowCount?: number;
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

/** The four sortable column keys, in header order (Tags is excluded). */
export type SortColumnKey = "sample" | "exposures" | "kept" | "status";

/** The Samples grid is always 6 columns: checkbox(0) Sample(1) Exposures(2)
 *  Kept(3) Tags(4) Status(5). The page always sets `checkboxColumn`, so the grid
 *  is 6-wide whenever roving is on. */
const GRID_COLS = 6;
/** Multi-widget cells (Enter/F2 → interaction mode): Exposures(2), Tags(4). */
const INTERACTION_COLS = new Set([2, 4]);

/**
 * The contact-sheet samples table: an elevated Card plate containing an
 * aligned column-header row (paper-sunk, hairline-strong bottom border) and
 * a slotted rows region. Column alignment is guaranteed because BOTH this
 * header and every SampleTableRow child use the SAME exported
 * `SAMPLE_TABLE_COLS` constant (not a re-derived arbitrary) — that shared
 * identity is what proves the header and rows are actually aligned.
 *
 * Narrow viewports (WCAG 1.4.10): the grid's intrinsic min-width (~1018px +
 * the CHECKBOX_TRACK_PX checkbox track) exceeds a 1024px viewport, so the
 * header rowgroup
 * AND the rows live together in ONE shared `overflow-x-auto` scroller
 * (`data-testid="sheet-scroll"`) — a narrow viewport SCROLLS rather than
 * clipping columns (Carbon/Polaris). The alignment invariant therefore also
 * depends on the shared scroller: header and rows must scroll as one. The
 * identity columns (checkbox + Sample) are sticky — see SampleTableRow.
 *
 * The inner wrapper is `min-w-min` (min-width: min-content), NOT min-w-max:
 * the grid's MIN-content width is exactly the sum of the fixed track minimums
 * (~1054px with the checkbox track) because the gallery cell is
 * overflow-x-auto/min-w-0 and the text cells clamp/truncate. max-content
 * would instead let the fr tracks expand to fit a fully-unwrapped
 * ThumbnailGallery (flex-nowrap) — a many-exposure corpus blew the scroller
 * out to ~3400px. With min-w-min the gallery keeps scrolling INTERNALLY and
 * the sheet scroller caps at the intrinsic track minimum. At wide viewports
 * the block wrapper's auto width fills the Card (min-width only binds when
 * the container is narrower) — zero wide-viewport change.
 *
 * ARIA tree: `role="table"` (or `role="grid"` when `roving`) lives on the inner
 * `min-w-min` wrapper (NOT the Card), so the scroll container sits BETWEEN the
 * visual plate and the table role tree. That keeps the owned-element chain
 * contiguous (table > rowgroup > row > columnheader/cell) with no role-bearing
 * wrapper in between. The empty slot renders OUTSIDE the table role tree (a
 * rowgroup may only own rows) and outside the scroller so it spans the Card
 * width.
 *
 * Children-slotting keeps SheetTable presentational: the PAGE maps its sample
 * data → SampleTableRow children (owning per-row handlers + selection state);
 * SheetTable only provides the plate chrome and the aligned header.
 *
 * Roving data grid (SA-ROVING): when `roving`, SheetTable owns a `useRovingGrid`
 * instance, wraps its whole subtree in `RovingGridProvider`, makes its OWN
 * header cells (row 0) roving cells, and puts the container keydown handler +
 * grid ref on the `role="grid"` element. Each SampleTableRow reads its
 * `tabIndex`/`ref`/`requestActivate` from the context by `(rowIndex, col)`. When
 * `roving` is absent the context default is INERT and nothing roves.
 *
 * Design-guard contract: `src/print/components/**` is NOT guard-exempt.
 * Only named-role classes, token color classes, and named layout utilities
 * appear here. `className` is PLACEMENT-ONLY.
 */
export function SheetTable({
  children,
  empty,
  checkboxColumn,
  sort,
  onSort,
  roving,
  dataRowCount,
  className,
}: SheetTableProps): JSX.Element {
  const childCount = React.Children.count(children);
  const isEmpty = childCount === 0;
  const showEmpty = isEmpty && empty != null;
  // The sticky offset for the Sample column derives from the SAME constant as
  // the checkbox grid track, so track width and frozen offset cannot diverge.
  const sampleLeft = checkboxColumn ? CHECKBOX_TRACK_PX : 0;

  // Roving grid owner (SA-ROVING). Hook is always called (rules of hooks); the
  // provider is only WIRED when `roving` is on — otherwise the children read the
  // INERT context default and render unchanged. rows = header(1) + data rows.
  const rovingDataRows = dataRowCount ?? childCount;
  const grid = useRovingGrid({
    rows: 1 + rovingDataRows,
    cols: GRID_COLS,
    interactionCols: INTERACTION_COLS,
  });
  const rv = roving ? grid : null;

  // aria-sort for a sortable column: "ascending"/"descending" on the active
  // column, "none" on the other sortable headers. Only set when sorting is on.
  const ariaSortFor = (key: SortColumnKey): "ascending" | "descending" | "none" =>
    sort?.key === key ? (sort.dir === "asc" ? "ascending" : "descending") : "none";
  const activeDir = (key: SortColumnKey): ColumnSortDir | undefined =>
    sort?.key === key ? sort.dir : undefined;

  // Header-cell roving wiring helpers (row 0). For a SORT-BUTTON header the
  // tabindex + register go on the button (its single widget); for a ZERO-widget
  // header (the checkbox col 0, the Tags label col 4) they go on the
  // columnheader div itself.
  const headTabIndex = (col: number): number | undefined =>
    rv ? rv.tabIndexFor(0, col) : undefined;
  const headRef =
    (col: number) =>
    (el: HTMLElement | null): void => {
      rv?.registerCellEl(0, col, el);
    };
  // Pointer parity for header cells: mousedown makes the header cell the active
  // roving cell WITHOUT forcing focus (the native click already focuses the sort
  // button). NOT onFocus — an onFocus latch traps Tab inside the grid (BUG D).
  const headActivatePointer = (col: number): (() => void) | undefined =>
    rv ? () => rv.requestActivate(0, col, { focus: false }) : undefined;

  const content = (
    <Card
      elevated
      data-testid="sheet-table"
      className={cx("overflow-hidden", className)}
    >
      {/* ONE shared horizontal scroller for header + rows (alignment invariant).
          tabIndex=-1 opts this overflow-x-auto scroller OUT of Chrome's automatic
          keyboard Tab-stop for scrollable containers — at narrow viewports it
          overflows and would otherwise become a Tab stop, trapping focus before
          the grid is exited (SA-ROVING BUG D). Horizontal scroll stays reachable
          via the grid's arrow-key cell navigation (cells scroll into view). */}
      <div data-testid="sheet-scroll" className="overflow-x-auto" tabIndex={-1}>
        <div
          role={roving ? "grid" : "table"}
          aria-label="Samples"
          className="min-w-min"
          {...(rv
            ? {
                ref: rv.gridRef,
                onKeyDown: rv.onContainerKeyDown,
              }
            : {})}
        >
          {/* Header: paper-sunk tint + hairline-strong bottom border.
              Header kickers are tone="soft" (NOT faint): ink-faint is 2.92:1
              on the sunk band — a WCAG 1.4.3 failure for these informational
              column labels. */}
          <div data-testid="sheet-head" role="rowgroup" className="bg-paper-sunk border-b border-hair-strong">
            <div
              role="row"
              className="grid"
              style={{ gridTemplateColumns: sampleTableCols(checkboxColumn ?? false) }}
            >
              {checkboxColumn && (
                <div
                  className="sticky left-0 z-10 bg-paper-sunk px-1 py-2.5 flex items-center justify-center"
                  role="columnheader"
                  aria-label="Select"
                  data-sticky="true"
                  {...(rv
                    ? {
                        tabIndex: headTabIndex(0),
                        ref: headRef(0),
                        onMouseDown: headActivatePointer(0),
                      }
                    : {})}
                />
              )}
              <div
                role="columnheader"
                data-sticky="true"
                {...(onSort ? { "aria-sort": ariaSortFor("sample") } : {})}
                className="sticky z-10 bg-paper-sunk border-r border-hair-strong"
                style={{ left: sampleLeft }}
                {...(rv && !onSort
                  ? {
                      tabIndex: headTabIndex(1),
                      ref: headRef(1),
                      onMouseDown: headActivatePointer(1),
                    }
                  : {})}
              >
                {onSort ? (
                  <ColumnSortButton
                    label="Sample"
                    active={activeDir("sample")}
                    onClick={() => onSort("sample")}
                    {...(rv
                      ? {
                          tabIndex: headTabIndex(1),
                          ref: headRef(1),
                          onMouseDown: headActivatePointer(1),
                        }
                      : {})}
                  />
                ) : (
                  <Kicker tone="soft" className="px-4 py-2.5">Sample</Kicker>
                )}
              </div>
              {onSort ? (
                <>
                  <div role="columnheader" aria-sort={ariaSortFor("exposures")}>
                    <ColumnSortButton
                      label="Exposures"
                      active={activeDir("exposures")}
                      onClick={() => onSort("exposures")}
                      {...(rv
                        ? {
                            tabIndex: headTabIndex(2),
                            ref: headRef(2),
                            onMouseDown: headActivatePointer(2),
                          }
                        : {})}
                    />
                  </div>
                  <div role="columnheader" aria-sort={ariaSortFor("kept")}>
                    <ColumnSortButton
                      label="Kept"
                      active={activeDir("kept")}
                      onClick={() => onSort("kept")}
                      {...(rv
                        ? {
                            tabIndex: headTabIndex(3),
                            ref: headRef(3),
                            onMouseDown: headActivatePointer(3),
                          }
                        : {})}
                    />
                  </div>
                  {/* Tags is multi-valued → never sortable: a plain label, no aria-sort.
                      Zero-widget header → the columnheader div itself is the roving cell. */}
                  <Kicker
                    tone="soft"
                    role="columnheader"
                    className="px-4 py-2.5"
                    {...(rv
                      ? {
                          tabIndex: headTabIndex(4),
                          ref: headRef(4),
                          onMouseDown: headActivatePointer(4),
                        }
                      : {})}
                  >
                    Tags
                  </Kicker>
                  <div role="columnheader" aria-sort={ariaSortFor("status")}>
                    <ColumnSortButton
                      label="Status"
                      active={activeDir("status")}
                      onClick={() => onSort("status")}
                      {...(rv
                        ? {
                            tabIndex: headTabIndex(5),
                            ref: headRef(5),
                            onMouseDown: headActivatePointer(5),
                          }
                        : {})}
                    />
                  </div>
                </>
              ) : (
                <>
                  <Kicker tone="soft" role="columnheader" className="px-4 py-2.5">Exposures</Kicker>
                  <Kicker tone="soft" role="columnheader" className="px-4 py-2.5">Kept</Kicker>
                  <Kicker tone="soft" role="columnheader" className="px-4 py-2.5">Tags</Kicker>
                  <Kicker tone="soft" role="columnheader" className="px-4 py-2.5">Status</Kicker>
                </>
              )}
            </div>
          </div>

          {/* Rows region */}
          {!showEmpty && (
            <div data-testid="sheet-rows" role="rowgroup">
              {children}
            </div>
          )}
        </div>
      </div>

      {showEmpty && <div data-testid="sheet-empty">{empty}</div>}
    </Card>
  );

  // Only wrap in the provider when roving is on; otherwise the INERT default
  // context keeps every consumer (incl. SampleTableRow's story) unchanged.
  return rv ? (
    <RovingGridProvider value={rv as RovingGridValue}>{content}</RovingGridProvider>
  ) : (
    content
  );
}
