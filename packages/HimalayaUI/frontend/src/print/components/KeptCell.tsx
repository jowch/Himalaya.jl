export interface KeptCellProps {
  kept: number;
  total: number;
  dropped?: number;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function KeptCell({
  kept,
  total,
  dropped,
  className,
}: KeptCellProps): JSX.Element {
  const showDropped = typeof dropped === "number" && dropped > 0;

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
      {showDropped && (
        // WCAG 1.4.3: this ~10.5px informational caption read 4.40:1 in terracotta
        // (text-accent) on paper-sunk — under the 4.5:1 AA-normal floor. ink-soft
        // measures 6.00:1 on paper-sunk (the F-CONTRAST precedent: small
        // informational text takes a darker ink token, not a hue). The word
        // "dropped" carries the verdict; the sibling kept-total uses the same tone.
        <span
          data-role="kept-dropped"
          className="block whitespace-nowrap text-ink-soft text-xs font-semibold mt-0.5"
        >
          {dropped} dropped
        </span>
      )}
    </div>
  );
}
