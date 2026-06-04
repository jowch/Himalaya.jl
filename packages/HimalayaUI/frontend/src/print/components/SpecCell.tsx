import { Chip } from "../ui/Chip";

export interface SpecCellProps {
  name: string;
  sampleId: string;
  screened?: boolean;
  /** When provided, promotes the name to a keyboard-accessible button that
   *  opens the loupe view for this sample. When absent, the name renders as
   *  a static span (other consumers unaffected). */
  onOpenLoupe?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/** The samples-table identity cell: sample name + id, plus a quiet "Needs
 *  review" chip on UN-screened rows. We flag the EXCEPTION (the actionable,
 *  not-yet-screened state) rather than marking every screened row — a leading
 *  checkbox-in-a-circle reads as row SELECTION (Carbon/Polaris). The chip sits
 *  inline on the id line, right of the id: the id may `truncate` to yield room
 *  while the chip holds its size (`flex-shrink-0`), so the column still scans. */
export function SpecCell({
  name,
  sampleId,
  screened,
  onOpenLoupe,
  className,
}: SpecCellProps): JSX.Element {
  return (
    <div
      data-testid="spec-cell"
      className={`min-w-0${className ? ` ${className}` : ""}`}
    >
      {onOpenLoupe ? (
        <button
          type="button"
          data-role="spec-name"
          onClick={onOpenLoupe}
          className="text-body font-semibold block leading-tight line-clamp-2 text-left hover:text-print-accent"
        >
          {name}
        </button>
      ) : (
        <span data-role="spec-name" className="text-body font-semibold block leading-tight line-clamp-2">
          {name}
        </span>
      )}
      <div className="flex items-center gap-2 mt-0.5 min-w-0">
        <span
          data-role="spec-id"
          className="text-data text-ink-faint truncate min-w-0"
        >
          {sampleId}
        </span>
        {!screened && (
          <Chip variant="static" size="sm" testId="needs-review" className="flex-shrink-0">
            Needs review
          </Chip>
        )}
      </div>
    </div>
  );
}
