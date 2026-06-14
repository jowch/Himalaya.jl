import { phaseColor } from "../../phases";

export type ScoreBarSize = "bar" | "compact";

interface ScoreBarProps {
  value: number;
  phase: string;
  size?: ScoreBarSize;
  className?: string;
}

const sizeClass: Record<ScoreBarSize, string> = {
  bar: "h-1 w-full",
  // FO-SCOREBAR-SLIVER: 46px was too narrow for the ranked candidate comparison
  // (0.28 vs 0.32 differed by ~1px). 88px + a touch taller makes the gap legible.
  compact: "h-[5px] w-[88px]",
};

export function ScoreBar({
  value,
  phase,
  size = "bar",
  className = "",
}: ScoreBarProps): JSX.Element {
  const pct = `${Math.round(Math.min(1, Math.max(0, value)) * 100)}%`;
  return (
    <div className={`overflow-hidden rounded-full bg-hair ${sizeClass[size]} ${className}`}>
      <i
        data-score-bar
        data-phase={phase}
        className="block h-full"
        style={{ width: pct, background: phaseColor(phase) }}
      />
    </div>
  );
}
