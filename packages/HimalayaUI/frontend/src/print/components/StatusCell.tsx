import { Dot } from "../ui/Dot";
import { PhaseChip } from "../ui/PhaseChip";

export interface StatusCellProps {
  phase?: string | null;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function StatusCell({ phase, className }: StatusCellProps): JSX.Element {
  const hasPhase = typeof phase === "string" && phase.length > 0;

  return (
    <div
      data-testid="status-cell"
      className={className ?? ""}
    >
      {hasPhase ? (
        <PhaseChip phase={phase as string} />
      ) : (
        <span data-role="status-unset" className="text-ink-faint text-xs">
          <Dot tone="muted" size="xs" className="mr-1.5 align-middle" />
          Not indexed
        </span>
      )}
    </div>
  );
}
