import { useId } from "react";
import { stripRange } from "../compare";
import { PlotFrame } from "../../print/plot/PlotFrame";
import { makeAxis } from "../../print/plot/projection";
import { viewAxes, viewIntensity, type PlotView } from "../plotView";

interface Props {
  q: readonly number[];
  d: readonly number[];
  label: string;
  xDomain: [number, number] | null;
  xExtent: [number, number];
  peakQs: readonly number[];
  view?: PlotView;
  height?: number;
}

const MARGINS = { top: 18, right: 24, bottom: 6, left: 60 };

export function SubtractionStrip({ q, d: raw, label, xDomain, xExtent, peakQs, view = "log", height = 110 }: Props): JSX.Element {
  const dom = xDomain ?? xExtent;
  const d = raw.map((v, i) => viewIntensity(view, q[i]!, v));
  const { xType } = viewAxes(view);
  const clipId = `workspace-subtraction-clip-${useId().replace(/:/g, "")}`;
  const [lo, hi] = stripRange(q, d, dom, peakQs);
  const pad = (hi - lo) * 0.1;
  return (
    <PlotFrame
      height={height}
      margins={MARGINS}
      data-testid="workspace-subtraction"
      svgRole="img"
      svgLabel={`Difference: current trace ${label}`}
      render={(dims) => {
        const x = makeAxis(dom, [0, dims.plotWidth], xType);
        const y = makeAxis([lo - pad, hi + pad], [dims.plotHeight, 0], "linear");
        let path = "";
        q.forEach((qq, i) => {
          const px = x.to(qq);
          const py = y.to(d[i]!);
          if (!Number.isFinite(px) || !Number.isFinite(py)) return;
          path += `${path ? "L" : "M"}${px.toFixed(1)},${py.toFixed(1)}`;
        });
        return (
          <>
            <defs>
              <clipPath id={clipId}>
                <rect x={0} y={0} width={dims.plotWidth} height={dims.plotHeight} />
              </clipPath>
            </defs>
            <text x={0} y={-6} className="fill-ink-faint text-caption">{label}</text>
            <line x1={0} x2={dims.plotWidth} y1={y.to(0)} y2={y.to(0)} stroke="var(--color-hair-strong)" strokeDasharray="3 3" />
            <path data-role="difference-line" d={path} fill="none" stroke="var(--color-accent)" strokeWidth={1.4} clipPath={`url(#${clipId})`} />
          </>
        );
      }}
    />
  );
}
