import type { ReactNode } from "react";
import {
  Card,
  CheckCircle,
  PhaseChip,
  BonnetBadge,
  ScoreBar,
} from "../ui";

export interface CandidateRowProps {
  /** Candidate phase name (drives the PhaseChip + ScoreBar color). */
  phase: string;
  /** Fit score 0..1, or null for a custom/unscored candidate (renders "custom", no bar). */
  score: number | null;
  /** One-line rationale (the mockup .c-cover note). */
  why: ReactNode;
  /** Show the Gauss–Bonnet flag. */
  bonnet?: boolean;
  /** In the assignment cart. */
  selected?: boolean;
  /** Toggle add/remove. */
  onToggle?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export interface CandidateListProps {
  children: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function CandidateRow({
  phase,
  score,
  why,
  bonnet,
  selected = false,
  onToggle,
  className,
}: CandidateRowProps): JSX.Element {
  return (
    <Card
      as="button"
      selected={selected}
      aria-pressed={selected}
      aria-label={`${phase}${selected ? ", in assignment" : ""}`}
      className={`w-full flex items-center gap-3 px-3 py-2.5 text-left${className ? ` ${className}` : ""}`}
      {...(onToggle ? { onClick: onToggle } : {})}
    >
      {/* selection mark */}
      <CheckCircle checked={selected} />

      {/* body */}
      <div className="flex-1 min-w-0">
        <div className="flex items-center gap-1.5">
          <PhaseChip phase={phase} variant="tint" size="md" />
          {bonnet && <BonnetBadge />}
        </div>
        <div className="text-caption text-ink-faint mt-0.5">{why}</div>
      </div>

      {/* score block */}
      <div className="flex flex-col items-end flex-shrink-0">
        <span className="text-data-strong text-ink font-mono font-bold">
          {score == null ? "custom" : score.toFixed(2)}
        </span>
        {score != null && (
          <ScoreBar value={score} phase={phase} size="compact" className="mt-1" />
        )}
      </div>
    </Card>
  );
}

export function CandidateList({ children, className }: CandidateListProps): JSX.Element {
  return (
    <div className={`flex flex-col gap-2${className ? ` ${className}` : ""}`}>
      {children}
    </div>
  );
}
