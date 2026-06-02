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
      className={`font-mono${className ? ` ${className}` : ""}`}
    >
      <div>
        <span data-role="kept-count" className="text-data text-ink font-semibold">
          {kept}
        </span>
        <span data-role="kept-total" className="text-data text-ink-faint">
          {" / "}{total}
        </span>
      </div>
      {showDropped && (
        <span
          data-role="kept-dropped"
          className="block text-accent text-xs font-semibold mt-0.5"
        >
          {dropped} dropped
        </span>
      )}
    </div>
  );
}
