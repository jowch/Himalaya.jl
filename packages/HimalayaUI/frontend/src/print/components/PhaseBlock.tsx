import type { ReactNode } from "react";
import { PhaseLabel, ScoreBar, IconButton } from "../ui";

export interface PhaseBlockProps {
  phase: string;
  score: number;        // 0..1
  meta: ReactNode;      // e.g. "a = 14.2 nm · 5 reflections"
  onRemove?: () => void;
  /** Optional series name shown below the score bar. */
  series?: ReactNode;
  /** PLACEMENT-ONLY. Appended last. */
  className?: string;
}

export function PhaseBlock({
  phase,
  score,
  meta,
  onRemove,
  series,
  className,
}: PhaseBlockProps): JSX.Element {
  return (
    <div
      data-testid="phase-block"
      className={`py-3${className ? ` ${className}` : ""}`}
    >
      <div className="flex items-baseline gap-2">
        {/* Name is 22px serif in the mockup; the named type scale has no 22px
            step, so we snap to the nearest serif role (text-headline, 19px). */}
        <PhaseLabel phase={phase} className="text-headline flex-1 min-w-0 truncate">
          {phase}
        </PhaseLabel>
        <PhaseLabel phase={phase} className="text-data-strong font-bold shrink-0">
          {score.toFixed(2)}
        </PhaseLabel>
        <IconButton
          label={`Remove ${phase}`}
          dismiss
          tone="ghost"
          className="shrink-0 hover:text-print-accent"
          {...(onRemove ? { onClick: onRemove } : {})}
        />
      </div>
      <div className="text-data text-ink-soft mt-1">{meta}</div>
      <ScoreBar value={score} phase={phase} size="bar" className="mt-2" />
      {series != null && (
        <div className="text-caption text-ink-faint mt-1.5">
          series <b className="text-ink-soft font-semibold font-mono">{series}</b>
        </div>
      )}
    </div>
  );
}
