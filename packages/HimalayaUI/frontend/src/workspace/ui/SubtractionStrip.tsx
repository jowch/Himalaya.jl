import { useId } from "react";
import { stripRange } from "../compare";
import { PlotFrame } from "../../print/plot/PlotFrame";
import { makeAxis } from "../../print/plot/projection";
import { viewAxes, viewIntensity, type PlotView } from "../plotView";
import { PLOT_MARGINS } from "../plotScale";
import { linePath } from "./WorkspacePlot";

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

const MARGINS = { top: 6, right: PLOT_MARGINS.right, bottom: 8, left: PLOT_MARGINS.left };

export function SubtractionStrip({ q, d: raw, label, xDomain, xExtent, peakQs, view = "log", height = 96 }: Props): JSX.Element {
  const dom = xDomain ?? xExtent;
  const d = raw.map((v, i) => viewIntensity(view, q[i]!, v));
  const { xType } = viewAxes(view);
  const clipId = `workspace-subtraction-clip-${useId().replace(/:/g, "")}`;
  const [rlo, rhi] = stripRange(q, d, dom, peakQs);
  const lo = Math.min(rlo, 0);
  const hi = Math.max(rhi, 0);
  const pad = (hi - lo) * 0.1 || 1;
  return (
    <div className="border-t border-dashed border-hair">
      <PlotFrame
        height={height}
        margins={MARGINS}
        data-testid="workspace-subtraction"
        svgRole="img"
        svgLabel={`Difference: current trace ${label}`}
        render={(dims) => {
          const x = makeAxis(dom, [0, dims.plotWidth], xType);
          const y = makeAxis([lo - pad, hi + pad], [dims.plotHeight, 0], "linear");
          return (
            <>
              <defs>
                <clipPath id={clipId}>
                  <rect x={0} y={-MARGINS.top} width={dims.plotWidth} height={height} />
                </clipPath>
              </defs>
              <line x1={0} x2={dims.plotWidth} y1={y.to(0)} y2={y.to(0)} stroke="var(--color-hair-strong)" strokeDasharray="3 3" />
              <path data-role="difference-line" d={linePath(q, d, x, y)} fill="none" stroke="var(--color-accent)" strokeWidth={1.4} clipPath={`url(#${clipId})`} />
              <text x={6} y={8} className="fill-ink-faint text-caption" stroke="var(--color-plate)" strokeWidth={4} paintOrder="stroke">{label}</text>
            </>
          );
        }}
      />
    </div>
  );
}
