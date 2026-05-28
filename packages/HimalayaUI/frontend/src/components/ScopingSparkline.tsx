import type { Trace } from "../api";
import { phaseColor } from "../phases";
import { sparklinePath, SPARK_W, SPARK_H } from "../lib/plot/sparkline";

interface Props {
  trace: Trace | undefined;
  phase: string | null;
}

/**
 * The per-row mini-trace on the scoping worksheet (series-scoping.html `.spark`).
 * Stroked + 10%-fill in the row's dominant phase colour; unindexed → neutral.
 * An empty frame stands in while the trace is loading or the sample is
 * unindexed, so rows keep a stable 76×28 rhythm.
 */
export function ScopingSparkline({ trace, phase }: Props): JSX.Element {
  const d = trace ? sparklinePath(trace) : null;
  const color = phase ? phaseColor(phase) : "var(--color-ink-faint)";
  if (!d) {
    return (
      <span
        data-testid="scoping-sparkline-empty"
        className="h-7 w-[76px] shrink-0 overflow-hidden rounded-sm border border-hair bg-paper-sunk"
        aria-hidden
      />
    );
  }
  const y0 = SPARK_H - 4;
  return (
    <span className="h-7 w-[76px] shrink-0 overflow-hidden rounded-sm border border-hair bg-paper-sunk">
      <svg viewBox={`0 0 ${SPARK_W} ${SPARK_H}`} className="block h-full w-full" aria-hidden>
        <line x1={4} y1={y0} x2={SPARK_W - 4} y2={y0} stroke="var(--color-hair)" strokeWidth={1} />
        <path d={`${d} L${SPARK_W - 4} ${y0} L4 ${y0} Z`} fill={color} opacity={0.1} />
        <path d={d} fill="none" stroke={color} strokeWidth={1.4} strokeLinejoin="round" />
      </svg>
    </span>
  );
}
