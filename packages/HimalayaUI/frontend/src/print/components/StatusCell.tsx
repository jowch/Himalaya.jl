import { Dot } from "../ui/Dot";
import { PhaseChip } from "../ui/PhaseChip";

export interface StatusCellProps {
  phase?: string | null;
  /** The sample has zero exposures: there is nothing to index, so the cell reads
   *  a terminal "No exposures" status instead of the "Not indexed" invitation
   *  (which would imply an action the empty sample can't take). SA-ZEROEXP. */
  noExposures?: boolean;
  /** This cell sits inside the Index door button (the sample is indexable).
   *  SA-RESCORE3 F9: an unindexed door then reads as an INVITATION ("Index" + a
   *  → chevron) rather than a dead "Not indexed" status, so the door visibly
   *  affords its click instead of looking like a terminal label. An indexed door
   *  keeps the PhaseChip (already an obvious thing to open). */
  door?: boolean;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function StatusCell({ phase, noExposures = false, door = false, className }: StatusCellProps): JSX.Element {
  const hasPhase = typeof phase === "string" && phase.length > 0;

  return (
    <div
      data-testid="status-cell"
      className={className ?? ""}
    >
      {noExposures ? (
        <span data-role="status-no-exposures" className="text-ink-soft text-xs whitespace-nowrap">
          <Dot tone="muted" size="xs" className="mr-1.5 align-middle" />
          No exposures
        </span>
      ) : hasPhase ? (
        <PhaseChip phase={phase as string} />
      ) : door ? (
        <span data-role="status-index-invite" className="text-accent text-xs font-semibold whitespace-nowrap">
          Index
          <span aria-hidden="true" className="ml-1">&rarr;</span>
        </span>
      ) : (
        <span data-role="status-unset" className="text-ink-soft text-xs whitespace-nowrap">
          <Dot tone="muted" size="xs" className="mr-1.5 align-middle" />
          Not indexed
        </span>
      )}
    </div>
  );
}
