import React from "react";
import type { ReactNode } from "react";
import { Card, Kicker, ColumnSortButton } from "../ui";
import type { ColumnSortDir } from "../ui";
import { sampleTableCols, CHECKBOX_TRACK_PX } from "./SampleTableRow";

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
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

/** The four sortable column keys, in header order (Tags is excluded). */
export type SortColumnKey = "sample" | "exposures" | "kept" | "status";

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
 * ARIA tree: `role="table"` lives on the inner `min-w-min` wrapper (NOT the
 * Card), so the scroll container sits BETWEEN the visual plate and the table
 * role tree. That keeps the owned-element chain contiguous
 * (table > rowgroup > row > columnheader/cell) with no role-bearing wrapper
 * in between. The empty slot renders OUTSIDE the table role tree (a rowgroup
 * may only own rows) and outside the scroller so it spans the Card width.
 *
 * Children-slotting keeps SheetTable presentational: the PAGE maps its sample
 * data → SampleTableRow children (owning per-row handlers + selection state);
 * SheetTable only provides the plate chrome and the aligned header.
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
  className,
}: SheetTableProps): JSX.Element {
  const isEmpty = React.Children.count(children) === 0;
  const showEmpty = isEmpty && empty != null;
  // The sticky offset for the Sample column derives from the SAME constant as
  // the checkbox grid track, so track width and frozen offset cannot diverge.
  const sampleLeft = checkboxColumn ? CHECKBOX_TRACK_PX : 0;

  // aria-sort for a sortable column: "ascending"/"descending" on the active
  // column, "none" on the other sortable headers. Only set when sorting is on.
  const ariaSortFor = (key: SortColumnKey): "ascending" | "descending" | "none" =>
    sort?.key === key ? (sort.dir === "asc" ? "ascending" : "descending") : "none";
  const activeDir = (key: SortColumnKey): ColumnSortDir | undefined =>
    sort?.key === key ? sort.dir : undefined;

  return (
    <Card
      elevated
      data-testid="sheet-table"
      className={cx("overflow-hidden", className)}
    >
      {/* ONE shared horizontal scroller for header + rows (alignment invariant). */}
      <div data-testid="sheet-scroll" className="overflow-x-auto">
        <div role="table" aria-label="Samples" className="min-w-min">
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
                />
              )}
              <div
                role="columnheader"
                data-sticky="true"
                {...(onSort ? { "aria-sort": ariaSortFor("sample") } : {})}
                className="sticky z-10 bg-paper-sunk border-r border-hair-strong"
                style={{ left: sampleLeft }}
              >
                {onSort ? (
                  <ColumnSortButton
                    label="Sample"
                    active={activeDir("sample")}
                    onClick={() => onSort("sample")}
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
                    />
                  </div>
                  <div role="columnheader" aria-sort={ariaSortFor("kept")}>
                    <ColumnSortButton
                      label="Kept"
                      active={activeDir("kept")}
                      onClick={() => onSort("kept")}
                    />
                  </div>
                  {/* Tags is multi-valued → never sortable: a plain label, no aria-sort. */}
                  <Kicker tone="soft" role="columnheader" className="px-4 py-2.5">Tags</Kicker>
                  <div role="columnheader" aria-sort={ariaSortFor("status")}>
                    <ColumnSortButton
                      label="Status"
                      active={activeDir("status")}
                      onClick={() => onSort("status")}
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
}
