function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

interface ProgressBarProps {
  value: number;
  total: number;
  label?: string;
  className?: string;
}

/** Accent capsule track (mockup `.pbar`): a 4px `rounded-full bg-hair` rail whose
 *  `bg-accent` fill width = value/total %. Progress is carried by the fill WIDTH
 *  (and the ARIA value), never by hue alone — the second channel is position, so it
 *  reads in grayscale (A). Terracotta `accent` is the one rationed progress mark (B).
 *  Non-interactive indicator, so C/D are N/A, but it exposes `role="progressbar"` +
 *  `aria-valuenow/min/max` (+ optional `aria-label`) for a11y. H — flat (rail tonal
 *  step + fill, no shadow/gradient). */
export function ProgressBar({
  value,
  total,
  label,
  className = "",
}: ProgressBarProps): JSX.Element {
  const pct = total > 0 ? Math.min(100, (value / total) * 100) : 0;
  return (
    <div
      data-testid="progressbar"
      role="progressbar"
      aria-valuenow={value}
      aria-valuemin={0}
      aria-valuemax={total}
      aria-label={label}
      className={cx("h-1 rounded-full bg-hair overflow-hidden", className)}
    >
      <span className="block h-full bg-accent" style={{ width: `${pct}%` }} />
    </div>
  );
}
