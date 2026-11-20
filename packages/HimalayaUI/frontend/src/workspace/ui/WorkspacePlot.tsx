import type { Peak, Trace } from "../../api";
import { TracePlot } from "../../print/plot/TracePlot";
import { toTraceModel } from "../../print/pages/focusAdapters";

export interface CombLine {
  id: number;
  qs: readonly number[];
  color: string;
  emphasis: boolean;
}

interface Props {
  trace: Trace | undefined;
  peaks: readonly Peak[];
  phase: string | null;
  combs: readonly CombLine[];
  height?: number;
}

export function WorkspacePlot({ trace, peaks, phase, combs, height = 360 }: Props): JSX.Element {
  if (!trace || trace.q.length < 2) {
    return <div className="rounded-sm bg-paper-sunk" style={{ height }} data-testid="workspace-trace-empty" />;
  }
  return (
    <TracePlot
      trace={toTraceModel(trace, [...peaks], phase)}
      height={height}
      xType="linear"
      axes
      paperColor="var(--color-paper)"
      figureLabel="Integration trace: intensity vs q"
      data-testid="workspace-trace"
      overlay={({ projection, dims }) => (
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
      )}
    />
  );
}
