import { useId } from "react";
import type { Peak, Trace } from "../../api";
import { TracePlot } from "../../print/plot/TracePlot";
import { toTraceModel } from "../../print/plot/traceModel";
import type { PeakRef } from "../../lib/peakRef";
import type { Projection } from "../../print/plot/projection";

export interface CombLine {
  id: number;
  qs: readonly number[];
  color: string;
  emphasis: boolean;
}

export interface CompareLine {
  id: number;
  name: string;
  q: readonly number[];
  I: readonly number[];
  dash: string;
  emphasis: boolean;
}

export const COMPARE_DASHES = ["6 4", "2 3", "10 3 2 3"] as const;

interface Props {
  trace: Trace | undefined;
  peaks: readonly Peak[];
  phase: string | null;
  combs: readonly CombLine[];
  compares?: readonly CompareLine[];
  xDomain?: [number, number] | null;
  logY?: boolean;
  onXDomain?: (d: [number, number] | null) => void;
  onAddPeak?: (q: number) => void;
  onClickPeak?: (peak: PeakRef) => void;
  height?: number;
}

function comparePath(c: CompareLine, projection: Projection): string {
  let d = "";
  let pen = false;
  c.q.forEach((q, i) => {
    const v = c.I[i]!;
    const x = projection.x.to(q);
    const y = projection.y.to(v);
    if ((!(v > 0) && projection.y.type === "log") || !Number.isFinite(x) || !Number.isFinite(y)) {
      pen = false;
      return;
    }
    d += `${pen ? "L" : "M"}${x.toFixed(1)},${y.toFixed(1)}`;
    pen = true;
  });
  return d;
}

export function WorkspacePlot({
  trace,
  peaks,
  phase,
  combs,
  compares = [],
  xDomain = null,
  logY = true,
  onXDomain,
  onAddPeak,
  onClickPeak,
  height = 360,
}: Props): JSX.Element {
  const clipId = `workspace-clip-${useId().replace(/:/g, "")}`;
  if (!trace || trace.q.length < 2) {
    return <div className="rounded-sm bg-paper-sunk" style={{ height }} data-testid="workspace-trace-empty" />;
  }
  return (
    <TracePlot
      trace={toTraceModel(trace, [...peaks], phase)}
      height={height}
      xType="linear"
      yType={logY ? "log" : "linear"}
      xDomain={xDomain}
      axes
      paperColor="var(--color-paper)"
      figureLabel="Integration trace: intensity vs q"
      data-testid="workspace-trace"
      {...(onXDomain
        ? {
            interaction: {
              onXDomain,
              doubleClickResets: false,
              ...(onAddPeak ? { onAddPeak } : {}),
              ...(onClickPeak ? { onClickPeak } : {}),
            },
          }
        : {})}
      overlay={({ projection, dims }) => (
        <>
          <defs>
            <clipPath id={clipId}>
              <rect x={0} y={0} width={dims.plotWidth} height={dims.plotHeight} />
            </clipPath>
          </defs>
          <g data-testid="workspace-compares" clipPath={`url(#${clipId})`} pointerEvents="none">
            {compares.map((c) => (
              <path
                key={c.id}
                data-compare={c.id}
                d={comparePath(c, projection)}
                fill="none"
                stroke="var(--color-ink-soft)"
                strokeOpacity={c.emphasis ? 0.9 : 0.6}
                strokeWidth={c.emphasis ? 1.6 : 1.1}
                strokeDasharray={c.dash}
              />
            ))}
          </g>
          <g data-testid="workspace-combs">
            {combs.map((c) =>
              c.qs.map((q) => {
                const x = projection.x.to(q);
                if (!Number.isFinite(x) || x < 0 || x > dims.plotWidth) return null;
                return (
                  <line
                    key={`${c.id}-${q}`}
                    data-comb={c.id}
                    x1={x}
                    x2={x}
                    y1={0}
                    y2={dims.plotHeight}
                    stroke={c.color}
                    strokeOpacity={c.emphasis ? 0.85 : 0.45}
                    strokeDasharray={c.emphasis ? undefined : "3 3"}
                  />
                );
              }),
            )}
          </g>
        </>
      )}
    />
  );
}
