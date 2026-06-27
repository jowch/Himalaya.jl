export interface KeptCellProps {
  kept: number;
  total: number;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/** The samples-table "Kept" cell: a plain `kept / total` count, where kept =
 *  not-dropped. Status is binary (dropped vs kept) so there is nothing to
 *  restore here — un-dropping happens on the frame, via the Keep verb. */
export function KeptCell({ kept, total, className }: KeptCellProps): JSX.Element {
  return (
    <div
      data-testid="kept-cell"
      className={`font-mono overflow-hidden${className ? ` ${className}` : ""}`}
    >
      <div className="whitespace-nowrap">
        <span data-role="kept-count" className="text-data text-ink font-semibold">
          {kept}
        </span>
        <span data-role="kept-total" className="text-data text-ink-soft">
          {" / "}{total}
        </span>
      </div>
    </div>
  );
}
