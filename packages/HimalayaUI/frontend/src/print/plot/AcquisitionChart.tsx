import type { JSX } from "react";
import type { AcqSession } from "../components/AcquisitionTimeline";

const H = 132;          // chart height (px)
const PAD_TOP = 8;
const PAD_BOTTOM = 22;  // room for the session label row
const BAR_W = 14;
const BAR_GAP = 4;      // between bars in a cluster
const CLUSTER_GAP = 22; // between session clusters
const PLOT_H = H - PAD_TOP - PAD_BOTTOM;

export interface AcquisitionChartProps {
  sessions: AcqSession[];
  className?: string;
}

/** Categorical bar cluster (per session, per load) on the trace-plot token
 *  language. Appearance-exempt render layer -- authors fill/stroke directly. */
export function AcquisitionChart({ sessions, className }: AcquisitionChartProps): JSX.Element {
  const maxCount = Math.max(1, ...sessions.flatMap((s) => s.loadFrameCounts));
  const yOf = (count: number): number => PLOT_H - (count / maxCount) * PLOT_H; // top edge of a bar

  // Walk left to right, accumulating x and recording per-cluster centers for labels.
  let x = 0;
  const clusters = sessions.map((s) => {
    const bars = s.loadFrameCounts.map((count) => {
      const bx = x;
      x += BAR_W + BAR_GAP;
      return { x: bx, count };
    });
    if (bars.length === 0) x += BAR_W; // empty session still occupies a slot
    const start = bars[0]?.x ?? x;
    const end = (bars[bars.length - 1]?.x ?? x) + BAR_W;
    const center = (start + end) / 2;
    x += CLUSTER_GAP;
    return { label: s.label, bars, center };
  });
  const width = Math.max(0, x - CLUSTER_GAP);

  return (
    <svg
      className={className}
      width="100%"
      viewBox={`0 0 ${Math.max(width, 1)} ${H}`}
      role="img"
      aria-label="Acquisition timeline -- exposures per load, grouped by session"
    >
      {/* baseline hairline in the plot's frame-edge token */}
      <line x1={0} y1={PAD_TOP + PLOT_H} x2={width} y2={PAD_TOP + PLOT_H}
            stroke="var(--color-hair)" strokeWidth={1} />
      {clusters.map((c) => (
        <g key={c.label}>
          {c.bars.map((b, i) => (
            <rect
              key={i}
              data-testid="acq-bar"
              x={b.x}
              y={PAD_TOP + yOf(b.count)}
              width={BAR_W}
              height={PLOT_H - yOf(b.count)}
              rx={2}
              fill="var(--color-accent)"
            />
          ))}
          <text
            x={c.center}
            y={H - 6}
            textAnchor="middle"
            fontSize={10}
            fontWeight={700}
            fill="var(--color-ink-faint)"
            style={{ textTransform: "uppercase", letterSpacing: "0.04em" }}
          >
            {c.label}
          </text>
        </g>
      ))}
    </svg>
  );
}
