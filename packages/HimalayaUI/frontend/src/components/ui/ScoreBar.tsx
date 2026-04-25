interface ScoreBarProps {
  score: number;
  color: string;
}

export function ScoreBar({ score, color }: ScoreBarProps): JSX.Element {
  const pct = `${Math.round(Math.min(1, Math.max(0, score)) * 100)}%`;
  return (
    <div className="h-0.5 w-full bg-bg-hover rounded-full overflow-hidden mt-1">
      <div
        data-score-bar
        className="h-full rounded-full"
        style={{ width: pct, background: color }}
      />
    </div>
  );
}
