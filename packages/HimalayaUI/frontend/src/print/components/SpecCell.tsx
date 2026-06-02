import { Chip } from "../ui/Chip";

export interface SpecCellProps {
  name: string;
  sampleId: string;
  screened?: boolean;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/** The samples-table identity cell: sample name + id, plus a quiet "Needs
 *  review" chip on UN-screened rows. We flag the EXCEPTION (the actionable,
 *  not-yet-screened state) rather than marking every screened row — a leading
 *  checkbox-in-a-circle reads as row SELECTION, and a mark floated right of the
 *  subtitle breaks vertical scannability (Carbon/Polaris). The chip sits in a
 *  fixed slot under the id, left-aligned, so the column still scans. */
export function SpecCell({
  name,
  sampleId,
  screened,
  className,
}: SpecCellProps): JSX.Element {
  return (
    <div
      data-testid="spec-cell"
      className={`min-w-0${className ? ` ${className}` : ""}`}
    >
      <span data-role="spec-name" className="text-body font-semibold block leading-tight line-clamp-2">
        {name}
      </span>
      <span
        data-role="spec-id"
        className="text-data text-ink-faint block truncate mt-0.5"
      >
        {sampleId}
      </span>
      {!screened && (
        <Chip variant="static" size="sm" testId="needs-review" className="mt-1.5">
          Needs review
        </Chip>
      )}
    </div>
  );
}
