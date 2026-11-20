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
  r2: number | null;
}

const MARGINS = { top: 8, right: 8, bottom: 26, left: 50 };

export function IndexingPlot({ series, empty }: { series: readonly IndexSeries[]; empty: string }): JSX.Element {
  const clipId = `workspace-index-clip-${useId().replace(/:/g, "")}`;
  if (series.length === 0) return <p className="text-caption">{empty}</p>;
  const xmax = Math.max(1, ...series.flatMap((s) => s.fit.points.map((p) => p.ratio))) * 1.08;
  const ymax = Math.max(0, ...series.flatMap((s) => s.fit.points.map((p) => Math.max(p.predicted, p.observed ?? 0)))) * 1.08 || 1;
  return (
    <div className="flex flex-col gap-1" data-testid="workspace-indexing">
      <PlotFrame
        height={180}
        margins={MARGINS}
        defaultWidth={268}
        svgRole="img"
        svgLabel="Indexing plot: peak q against reflection ratio"
        render={(dims) => {
          const x = makeAxis([0, xmax], [0, dims.plotWidth], "linear");
          const y = makeAxis([0, ymax], [dims.plotHeight, 0], "linear");
          const xt = scaleLinear().domain([0, xmax]).ticks(5).filter((t) => t > 0);
          const ys = scaleLinear().domain([0, ymax]);
          const yt = ys.ticks(4).filter((t) => t > 0);
          const yfmt = ys.tickFormat(4);
          return (
            <>
              <defs>
                <clipPath id={clipId}>
                  <rect x={0} y={0} width={dims.plotWidth} height={dims.plotHeight} />
                </clipPath>
              </defs>
              <line x1={0} x2={0} y1={0} y2={dims.plotHeight} stroke="var(--color-hair-strong)" />
              <line x1={0} x2={dims.plotWidth} y1={dims.plotHeight} y2={dims.plotHeight} stroke="var(--color-hair-strong)" />
              {xt.map((t) => (
                <text key={`x${t}`} x={x.to(t)} y={dims.plotHeight + 13} textAnchor="middle" className="fill-ink-faint text-caption">
                  {t}
                </text>
              ))}
              {yt.map((t) => (
                <text key={`y${t}`} x={-5} y={y.to(t) + 4} textAnchor="end" className="fill-ink-faint text-caption">
                  {yfmt(t)}
                </text>
              ))}
              <text x={dims.plotWidth} y={dims.plotHeight + 24} textAnchor="end" className="fill-ink-faint text-caption">
                q / q₁
              </text>
              <text transform={`translate(-38 ${dims.plotHeight / 2}) rotate(-90)`} textAnchor="middle" className="fill-ink-faint text-caption">
                q (Å⁻¹)
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
      {series.map((s) => (
        <p key={s.key} className="text-caption" data-testid="workspace-indexing-fit">
          <span className="font-semibold" style={{ color: s.color }}>
            {s.phase}
          </span>{" "}
          {s.lattice !== null ? `${s.lattice.toFixed(1)} Å · ` : ""}
          {s.r2 !== null ? `R² ${s.r2.toFixed(4)} · ` : ""}
          {s.fit.landed} of {s.fit.points.length} on peaks
        </p>
      ))}
    </div>
  );
}
