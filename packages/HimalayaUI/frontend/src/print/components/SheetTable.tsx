import React from "react";
import type { ReactNode } from "react";
import { Card, Kicker } from "../ui";
import { SAMPLE_TABLE_COLS } from "./SampleTableRow";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface SheetTableProps {
  /** SampleTableRow elements, already mapped/sorted by the page. */
  children?: ReactNode;
  /** Rendered (instead of the rows region) when there are zero children. */
  empty?: ReactNode;
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}

/**
 * The contact-sheet samples table: an elevated Card plate containing an
 * aligned column-header row (paper-sunk, hairline-strong bottom border) and
 * a slotted rows region. Column alignment is guaranteed because BOTH this
 * header and every SampleTableRow child use the SAME exported
 * `SAMPLE_TABLE_COLS` constant (not a re-derived arbitrary) — that shared
 * identity is what proves the header and rows are actually aligned.
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
  className,
}: SheetTableProps): JSX.Element {
  const isEmpty = React.Children.count(children) === 0;

  return (
    <Card elevated data-testid="sheet-table" className={cx("overflow-hidden", className)}>
      {/* Header: paper-sunk tint + hairline-strong bottom border */}
      <div data-testid="sheet-head" className="bg-paper-sunk border-b border-hair-strong">
        <div className="grid" style={{ gridTemplateColumns: SAMPLE_TABLE_COLS }}>
          <Kicker tone="faint" className="px-4 py-2.5">Sample</Kicker>
          <Kicker tone="faint" className="px-4 py-2.5">Exposures</Kicker>
          <Kicker tone="faint" className="px-4 py-2.5">Kept</Kicker>
          <Kicker tone="faint" className="px-4 py-2.5">Tags</Kicker>
          <Kicker tone="faint" className="px-4 py-2.5">Status</Kicker>
        </div>
      </div>

      {/* Rows region */}
      <div data-testid="sheet-rows">
        {isEmpty && empty != null ? (
          <div data-testid="sheet-empty">{empty}</div>
        ) : (
          children
        )}
      </div>
    </Card>
  );
}
