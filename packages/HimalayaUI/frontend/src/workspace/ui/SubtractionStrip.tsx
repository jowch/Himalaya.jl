import { PlotFrame } from "../../print/plot/PlotFrame";
import { makeAxis } from "../../print/plot/projection";

interface Props {
  q: readonly number[];
  d: readonly number[];
  label: string;
  xDomain: [number, number] | null;
  xExtent: [number, number];
  height?: number;
}

const MARGINS = { top: 18, right: 24, bottom: 6, left: 60 };

export function SubtractionStrip({ q, d, label, xDomain, xExtent, height = 110 }: Props): JSX.Element {
  const dom = xDomain ?? xExtent;
  let lo = 0;
  let hi = 0;
  q.forEach((x, i) => {
    if (x < dom[0] || x > dom[1]) return;
    lo = Math.min(lo, d[i]!);
    hi = Math.max(hi, d[i]!);
  });
  const pad = (hi - lo) * 0.1 || 1;
  return (
    <PlotFrame
      height={height}
      margins={MARGINS}
      data-testid="workspace-subtraction"
      svgRole="img"
      svgLabel={`Difference: current trace ${label}`}
      render={(dims) => {
        const x = makeAxis(dom, [0, dims.plotWidth], "linear");
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
              <clipPath id="workspace-subtraction-clip">
                <rect x={0} y={0} width={dims.plotWidth} height={dims.plotHeight} />
              </clipPath>
            </defs>
            <text x={0} y={-6} className="fill-ink-faint text-caption">{label}</text>
            <line x1={0} x2={dims.plotWidth} y1={y.to(0)} y2={y.to(0)} stroke="var(--color-hair-strong)" strokeDasharray="3 3" />
            <path d={path} fill="none" stroke="var(--color-accent)" strokeWidth={1.4} clipPath="url(#workspace-subtraction-clip)" />
          </>
        );
      }}
    />
  );
}
