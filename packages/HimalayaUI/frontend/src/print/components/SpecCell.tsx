import { forwardRef } from "react";
import { Chip } from "../ui/Chip";

export interface SpecCellProps {
  name: string;
  sampleId: string;
  /** Slot coordinate shown as a chip in the identity cell (e.g. "slot 5").
   *  When undefined, no slot chip renders. */
  slotIndex?: number;
  /** When provided, promotes the name to a keyboard-accessible button that
   *  opens the loupe view for this sample. When absent, the name renders as
   *  a static span (other consumers unaffected). */
  onOpenLoupe?: () => void;
  /** Override tabIndex on the name button. Absent → the button's natural tab order. */
  tabIndex?: number;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/** The samples-table identity cell: sample name + id, plus an optional "slot N"
 *  chip when a slot coordinate is provided. The chip sits inline on the id line,
 *  right of the id: the id may `truncate` to yield room while the chip holds its
 *  size (`flex-shrink-0`), so the column still scans. */
/** Forwards its ref to the focusable name button when `onOpenLoupe` is provided. */
export const SpecCell = forwardRef<HTMLButtonElement, SpecCellProps>(function SpecCell(
  { name, sampleId, slotIndex, onOpenLoupe, tabIndex, className },
  ref,
): JSX.Element {
  return (
    <div
      data-testid="spec-cell"
      className={`min-w-0${className ? ` ${className}` : ""}`}
    >
      {onOpenLoupe ? (
        <button
          ref={ref}
          type="button"
          data-role="spec-name"
          onClick={onOpenLoupe}
          {...(tabIndex !== undefined ? { tabIndex } : {})}
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
          className="text-data text-ink-soft truncate min-w-0"
        >
          {sampleId}
        </span>
        {slotIndex !== undefined && (
          <Chip variant="static" size="sm" testId="slot-chip" className="flex-shrink-0">slot {slotIndex}</Chip>
        )}
      </div>
    </div>
  );
});
