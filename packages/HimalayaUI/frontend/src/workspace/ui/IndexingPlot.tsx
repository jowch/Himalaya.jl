import { useId } from "react";
import { scaleLinear } from "d3-scale";
import { PlotFrame } from "../../print/plot/PlotFrame";
import { makeAxis } from "../../print/plot/projection";
import type { IndexFit } from "../indexing";

export interface IndexSeries {
  key: string;
  phase: string;
  color: string;
  fit: IndexFit;
  lattice: number | null;
  param?: string;
  r2: number | null;
}

const MARGINS = { top: 8, right: 8, bottom: 22, left: 34 };
const HEIGHT = 190;

export function IndexingPlot({ series }: { series: readonly IndexSeries[] }): JSX.Element {
  const clipId = `workspace-index-clip-${useId().replace(/:/g, "")}`;
  const empty = series.length === 0;
  const xmax = empty ? 4 : Math.max(1, ...series.flatMap((s) => s.fit.points.map((p) => p.ratio))) * 1.08;
  const ymax = empty ? 0.3 : Math.max(0, ...series.flatMap((s) => s.fit.points.map((p) => Math.max(p.predicted, p.observed ?? 0)))) * 1.08 || 1;
  return (
    <div className="flex flex-col" data-testid="workspace-indexing">
      <PlotFrame
        height={HEIGHT}
        margins={MARGINS}
        defaultWidth={272}
        svgRole="img"
        svgLabel="Indexing plot: peak q against reflection ratio"
        render={(dims) => {
          const x = makeAxis([0, xmax], [0, dims.plotWidth], "linear");
          const y = makeAxis([0, ymax], [dims.plotHeight, 0], "linear");
          const xt = scaleLinear().domain([0, xmax]).ticks(4);
          const yt = scaleLinear().domain([0, ymax]).ticks(4);
          return (
            <>
              <defs>
                <clipPath id={clipId}>
                  <rect x={0} y={0} width={dims.plotWidth} height={dims.plotHeight} />
                </clipPath>
              </defs>
              <path d={`M0 0V${dims.plotHeight}H${dims.plotWidth}`} fill="none" stroke="var(--color-hair)" />
              {xt.map((t) => (
                <text key={`x${t}`} x={x.to(t)} y={dims.plotHeight + 12} textAnchor="middle" className="fill-ink-faint text-caption">
                  {t}
                </text>
              ))}
              {yt.map((t) => (
                <text key={`y${t}`} x={-4} y={y.to(t) + 3} textAnchor="end" className="fill-ink-faint text-caption">
                  {t.toFixed(2)}
                </text>
              ))}
              <text x={dims.plotWidth / 2} y={dims.plotHeight + MARGINS.bottom - 1} textAnchor="middle" className="fill-ink-faint text-caption">
                q / q₁
              </text>
              <g clipPath={`url(#${clipId})`}>
                {series.map((s) => (
                  <g key={s.key} data-index-series={s.key}>
                    <line x1={x.to(0)} y1={y.to(0)} x2={x.to(xmax)} y2={y.to(s.fit.slope * xmax)} stroke={s.color} strokeWidth={1.2} />
                    {s.fit.points.map((p, i) => (
                      <circle key={`p${i}`} cx={x.to(p.ratio)} cy={y.to(p.predicted)} r={3} fill="none" stroke={s.color} />
                    ))}
                    {s.fit.points.map((p, i) =>
                      p.observed === null ? null : (
                        <circle key={`o${i}`} data-observed cx={x.to(p.ratio)} cy={y.to(p.observed)} r={3} fill={s.color} />
                      ),
                    )}
                  </g>
                ))}
              </g>
            </>
          );
        }}
      />
      <div className="mt-1.5 min-h-[3.2em] text-xs text-ink-faint">
        {series.map((s) => (
          <p key={s.key} className="m-0 flex flex-wrap gap-x-2" data-testid="workspace-indexing-fit">
            <span style={{ color: s.color }}>{s.phase}</span>{" "}
            {s.lattice !== null && <span>{`${s.param ?? "a"} = ${s.lattice.toFixed(1)} Å`}</span>}{" "}
            {s.r2 !== null && <span>{`R² ${s.r2.toFixed(4)}`}</span>}{" "}
            <span>{`${s.fit.landed}/${s.fit.points.length} on peaks`}</span>
          </p>
        ))}
      </div>
    </div>
  );
}
