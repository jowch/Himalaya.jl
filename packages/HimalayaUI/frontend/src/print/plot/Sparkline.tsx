import { makeAxis } from "./projection";
import { phaseColor } from "../../phases";

export interface SparklineProps {
  /** Measured trace (q ascending, I parallel). */
  trace: { q: number[]; I: number[] };
  /** Dominant phase for hue; null/undefined → unindexed (ink-faint). */
  phase?: string | null;
  /** PLACEMENT ONLY. */
  className?: string;
}

/** Degenerate guard: returns true when the trace has too few points, bad maxI,
 *  or a degenerate q range to render lines. */
function isDegenerate(trace: { q: number[]; I: number[] }): boolean {
  if (trace.q.length <= 1) return true;
  const positiveQ = trace.q.filter((v) => v > 0);
  if (positiveQ.length < 2) return true;
  // maxI must be taken over the positive-q positions only — a beamstop/detector
  // artefact at near-zero or negative q (filtered out of the x-mapping) would
  // otherwise squash the visible curve.
  const positiveQI = trace.I.filter((_, idx) => trace.q[idx]! > 0);
  if (positiveQI.length === 0) return true;
  const maxI = Math.max(...positiveQI);
  if (maxI <= 0) return true;
  const qMin = Math.min(...positiveQ);
  const qMax = Math.max(...positiveQ);
  if (qMin === qMax) return true;
  return false;
}

export function Sparkline({ trace, phase, className }: SparklineProps): JSX.Element {
  const color = phase ? phaseColor(phase) : "var(--color-ink-faint)";

  const containerClass = [
    "inline-block flex-shrink-0 overflow-hidden rounded-sm bg-paper-sunk border border-hair",
    className,
  ]
    .filter(Boolean)
    .join(" ");

  const degen = isDegenerate(trace);

  let lineD = "";
  let areaD = "";

  if (!degen) {
    const positiveQ = trace.q.filter((v) => v > 0);
    const qMin = Math.min(...positiveQ);
    const qMax = Math.max(...positiveQ);
    const maxI = Math.max(...trace.I.filter((_, idx) => trace.q[idx]! > 0));

    // x: log-q mapping over [4, 72]
    const x = makeAxis([qMin, qMax], [4, 72], "log");

    // y: linear over [0, maxI] → pixels [24, 5] (higher I = smaller y)
    const yScale = (i: number) => 24 - ((i / maxI) * (24 - 5));

    // Build point arrays, filtering to valid (q > 0)
    const points: Array<[number, number]> = [];
    for (let i = 0; i < trace.q.length; i++) {
      const q = trace.q[i]!;
      const intensity = trace.I[i]!;
      if (q <= 0) continue;
      const px = x.to(q);
      const py = yScale(intensity);
      if (!Number.isFinite(px) || !Number.isFinite(py)) continue;
      points.push([px, py]);
    }

    if (points.length >= 2) {
      // Line path: M x0,y0 L x1,y1 ...
      lineD = points
        .map(([px, py], idx) => `${idx === 0 ? "M" : "L"}${px.toFixed(2)},${py.toFixed(2)}`)
        .join(" ");

      // Area path: start at (4,24), trace across the points, end at (72,24), close
      const firstPx = points[0]![0].toFixed(2);
      const lastPx = points[points.length - 1]![0].toFixed(2);
      areaD =
        `M4,24 ` +
        (firstPx !== "4.00" ? `L${firstPx},24 ` : "") +
        points
          .map(([px, py]) => `L${px.toFixed(2)},${py.toFixed(2)}`)
          .join(" ") +
        (lastPx !== "72.00" ? ` L${lastPx},24` : "") +
        ` L72,24 Z`;
    }
  }

  const canDraw = !degen && lineD.length > 0;

  return (
    <span
      data-testid="sparkline"
      style={{ width: 76, height: 28 }}
      className={containerClass}
    >
      <svg
        viewBox="0 0 76 28"
        width="100%"
        height="100%"
        aria-hidden
        style={{ display: "block" }}
      >
        {/* 1. Baseline */}
        <line
          data-role="spark-baseline"
          x1={4}
          y1={24}
          x2={72}
          y2={24}
          stroke="var(--color-hair)"
          strokeWidth={1}
        />

        {/* 2. Filled area (only when non-degenerate) */}
        {canDraw && (
          <path
            data-role="spark-area"
            d={areaD}
            fill={color}
            fillOpacity={0.1}
          />
        )}

        {/* 3. Trace stroke (only when non-degenerate) */}
        {canDraw && (
          <path
            data-role="spark-line"
            d={lineD}
            fill="none"
            stroke={color}
            strokeWidth={1.4}
            strokeLinejoin="round"
          />
        )}
      </svg>
    </span>
  );
}
