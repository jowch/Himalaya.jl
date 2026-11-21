import { useCallback, useId, useRef, useState } from "react";
import type { Peak, Trace } from "../../api";
import { peakKey, type PeakRef } from "../../lib/peakRef";
import { PlotFrame, type PlotDims } from "../../print/plot/PlotFrame";
import { makeAxis, positiveExtent, type Axis1D, type ScaleType } from "../../print/plot/projection";
import { hitTestPeaks, zoomXDomain, PEAK_HIT_PX } from "../../print/plot/interaction";
import { viewAxes, viewIntensity, viewLabel, type PlotView } from "../plotView";
import { fitYDomain, PLOT_MARGINS, ticksFor } from "../plotScale";

export interface CombLine {
  id: number;
  qs: readonly number[];
  color: string;
}

export interface CompareLine {
  id: number;
  name: string;
  q: readonly number[];
  I: readonly number[];
  color: string;
  emphasis: boolean;
}

export const COMPARE_DASH = "5 4";
export const COMPARE_COLORS = [1, 2, 3, 4, 5].map((n) => `var(--color-compare-${n})`);
export const PLOT_HEIGHT = 330;
const Q_LABEL_ROOM = 70;

interface Props {
  trace: Trace | undefined;
  peaks: readonly Peak[];
  combs: readonly CombLine[];
  compares?: readonly CompareLine[];
  xDomain?: [number, number] | null;
  view?: PlotView;
  onXDomain?: (d: [number, number] | null) => void;
  onAddPeak?: (q: number) => void;
  onClickPeak?: (peak: PeakRef) => void;
  height?: number;
  onRetry?: () => void;
  retrying?: boolean;
  emptyMessage?: string;
}

interface Geometry {
  x: Axis1D;
  dims: PlotDims;
  xDomain: [number, number];
  xExtent: [number, number];
  xType: ScaleType;
}

export function linePath(q: readonly number[], I: readonly number[], x: Axis1D, y: Axis1D): string {
  let d = "";
  let pen = false;
  q.forEach((qq, i) => {
    const v = I[i]!;
    const px = x.to(qq);
    const py = y.to(v);
    if ((y.type === "log" && !(v > 0)) || (x.type === "log" && !(qq > 0)) || !Number.isFinite(px) || !Number.isFinite(py)) {
      pen = false;
      return;
    }
    d += `${pen ? "L" : "M"}${px.toFixed(1)},${py.toFixed(1)}`;
    pen = true;
  });
  return d;
}

function peakPoints(px: number): string {
  const pts: [number, number][] = [[px - 4, -8], [px + 4, -8], [px, -1]];
  return pts.map(([a, b]) => `${a.toFixed(1)},${b}`).join(" ");
}

function peakFill(p: Peak): string {
  if (p.excluded) return "none";
  return p.source === "manual" ? "var(--color-accent)" : "var(--color-ink)";
}

export function WorkspacePlot({
  trace,
  peaks,
  combs,
  compares = [],
  xDomain = null,
  view = "log",
  onXDomain,
  onAddPeak,
  onClickPeak,
  height = PLOT_HEIGHT,
  onRetry,
  retrying = false,
  emptyMessage,
}: Props): JSX.Element {
  const clipId = `workspace-clip-${useId().replace(/:/g, "")}`;
  const geometry = useRef<Geometry | null>(null);
  const [hoverQ, setHoverQ] = useState<number | null>(null);
  const { xType, yType } = viewAxes(view);

  const plotX = (px: number): number | null => {
    const g = geometry.current;
    if (!g) return null;
    const p = px - g.dims.margins.left;
    return p >= 0 && p <= g.dims.plotWidth ? p : null;
  };

  const onWheelPx = useCallback(
    (deltaY: number, px: number) => {
      const g = geometry.current;
      if (!g || !onXDomain) return;
      const cursorQ = g.x.invert(px - g.dims.margins.left);
      if (!Number.isFinite(cursorQ)) return;
      const next = zoomXDomain({ cursorQ, deltaY, current: g.xDomain, extent: g.xExtent, type: g.xType });
      if (next) onXDomain(next);
    },
    [onXDomain],
  );

  if (!trace && onRetry) {
    return (
      <div className="flex flex-col items-center justify-center gap-2 text-base text-ink-soft" style={{ height }} role="alert" data-testid="workspace-trace-error">
        <span>Couldn't load this sample</span>
        <button type="button" onClick={onRetry} disabled={retrying} className="rounded-sm border border-hair bg-paper px-2 py-0.5 text-sm text-ink hover:border-hair-strong disabled:text-ink-faint">
          {retrying ? "Retrying…" : "Retry"}
        </button>
      </div>
    );
  }

  const empty = trace ? (trace.q.length < 2 ? "This frame's trace has no points." : undefined) : emptyMessage ?? "";
  if (empty !== undefined || !trace) {
    return (
      <div role={empty ? "status" : undefined} className="flex items-center justify-center rounded-sm bg-paper-sunk px-6 text-center text-base text-ink-soft" style={{ height }} data-testid="workspace-trace-empty">
        {empty}
      </div>
    );
  }

  const xExtent = positiveExtent(trace.q);
  const shownX = xDomain ?? xExtent;
  const I = trace.I.map((v, i) => viewIntensity(view, trace.q[i]!, v));
  const shownCompares = compares.map((c) => ({ ...c, I: c.I.map((v, i) => viewIntensity(view, c.q[i]!, v)) }));
  const yDomain = fitYDomain({ q: trace.q, I }, shownCompares, shownX, yType);
  const yLabel = viewLabel(view);

  const lowX = Math.min(...shownX);
  const highX = Math.max(...shownX);
  const shownPeaks = peaks.filter((p) => p.q >= lowX && p.q <= highX);

  const onClickPx = (px: number, py: number): void => {
    const g = geometry.current;
    const p = plotX(px);
    const top = g ? g.dims.margins.top : 0;
    if (!g || p === null || py < top - 12 || py > top + g.dims.plotHeight) return;
    const q = g.x.invert(p);
    if (!Number.isFinite(q) || q <= 0) return;
    const hit = hitTestPeaks([...shownPeaks], p, (qq) => g.x.to(qq), PEAK_HIT_PX);
    if (hit && onClickPeak) onClickPeak({ id: hit.id, source: hit.source });
    else if (!hit && onAddPeak && py >= top) onAddPeak(q);
  };

  const onPointerMovePx = (px: number, py: number): void => {
    const g = geometry.current;
    const p = plotX(px);
    if (!g || p === null || py < 0 || py > g.dims.margins.top + g.dims.plotHeight) {
      setHoverQ(null);
      return;
    }
    const q = g.x.invert(p);
    setHoverQ(Number.isFinite(q) && q > 0 ? q : null);
  };

  return (
    <PlotFrame
      height={height}
      margins={PLOT_MARGINS}
      data-testid="workspace-trace"
      svgRole={onClickPeak ? "group" : "img"}
      svgLabel="Integration trace: intensity against q"
      onPointerMovePx={onPointerMovePx}
      onPointerLeave={() => setHoverQ(null)}
      {...(onXDomain ? { onWheelPx, onClickPx } : {})}
      render={(dims) => {
        const x = makeAxis(shownX, [0, dims.plotWidth], xType);
        const y = makeAxis(yDomain, [dims.plotHeight, 0], yType);
        geometry.current = { x, dims, xDomain: shownX, xExtent, xType };
        const xTicks = ticksFor(shownX, xType, Math.max(3, Math.floor(dims.plotWidth / 90)), false);
        const yTicks = ticksFor(yDomain, yType, 5, true);
        const inX = (q: number): boolean => q >= lowX && q <= highX;
        return (
          <>
            <defs>
              <clipPath id={clipId}>
                <rect x={0} y={-4} width={dims.plotWidth} height={dims.plotHeight + 4} />
              </clipPath>
            </defs>
            <g data-role="grid" stroke="var(--color-hair)" strokeOpacity={0.6} pointerEvents="none">
              {xTicks.map((t) => (
                <line key={`x${t.value}`} x1={x.to(t.value)} x2={x.to(t.value)} y1={0} y2={dims.plotHeight} />
              ))}
              {yTicks.map((t) => (
                <line key={`y${t.value}`} x1={0} x2={dims.plotWidth} y1={y.to(t.value)} y2={y.to(t.value)} />
              ))}
            </g>
            <g data-role="axis-bottom" className="fill-ink-faint text-caption" textAnchor="middle" aria-hidden="true">
              {xTicks.filter((t) => x.to(t.value) < dims.plotWidth - Q_LABEL_ROOM).map((t) => (
                <text key={t.value} x={x.to(t.value)} y={dims.plotHeight + 15}>
                  {t.label}
                </text>
              ))}
              <text x={dims.plotWidth} y={dims.plotHeight + 15} textAnchor="end">
                q (Å⁻¹)
              </text>
            </g>
            <g data-role="axis-left" className="fill-ink-faint text-caption" textAnchor="end" aria-hidden="true">
              {yTicks.map((t) => (
                <text key={t.value} x={-6} y={y.to(t.value) + 4}>
                  {t.label}
                </text>
              ))}
            </g>
            <g clipPath={`url(#${clipId})`} pointerEvents="none">
              <g data-testid="workspace-compares">
                {shownCompares.map((c) => (
                  <path
                    key={c.id}
                    data-compare={c.id}
                    d={linePath(c.q, c.I, x, y)}
                    fill="none"
                    stroke={c.color}
                    strokeWidth={c.emphasis ? 1.6 : 1.1}
                    strokeDasharray={COMPARE_DASH}
                  />
                ))}
              </g>
              <g data-role="trace-line">
                <path d={linePath(trace.q, I, x, y)} fill="none" stroke="var(--color-ink)" strokeWidth={1.8} strokeLinejoin="round" />
              </g>
              <g data-testid="workspace-combs">
                {combs.map((c) =>
                  c.qs.filter(inX).map((q) => (
                    <line
                      key={`${c.id}-${q}`}
                      data-comb={c.id}
                      x1={x.to(q)}
                      x2={x.to(q)}
                      y1={10}
                      y2={dims.plotHeight}
                      stroke={c.color}
                      strokeDasharray="2 3"
                    />
                  )),
                )}
              </g>
            </g>
            <g data-role="axis-title" className="fill-ink-faint text-caption" stroke="var(--color-plate)" strokeWidth={4} paintOrder="stroke" aria-hidden="true">
              {yLabel && (
                <text x={6} y={12}>
                  {yLabel}
                </text>
              )}
            </g>
            <g data-role="peak-glyph">
              {shownPeaks.map((p) => (
                <polygon
                  key={peakKey(p)}
                  data-peak-key={peakKey(p)}
                  data-excluded={p.excluded || undefined}
                  points={peakPoints(x.to(p.q))}
                  fill={peakFill(p)}
                  stroke={p.excluded ? "var(--color-ink-faint)" : "none"}
                  pointerEvents="all"
                  {...(onClickPeak
                    ? {
                        role: "button",
                        tabIndex: 0,
                        "aria-label": `${p.source === "manual" ? "Manual" : "Auto"} peak at q = ${p.q.toFixed(4)} Å⁻¹${p.excluded ? ", off" : ""}`,
                        className: "cursor-pointer outline-none focus-visible:outline focus-visible:outline-2 focus-visible:outline-accent",
                        onKeyDown: (ev: React.KeyboardEvent) => {
                          if (ev.key !== "Enter" && ev.key !== " ") return;
                          ev.preventDefault();
                          ev.stopPropagation();
                          onClickPeak({ id: p.id, source: p.source });
                        },
                      }
                    : {})}
                />
              ))}
            </g>
            {hoverQ !== null && (
              <text
                data-role="q-readout"
                x={6}
                y={dims.plotHeight - 6}
                className="fill-ink-soft text-caption"
                stroke="var(--color-plate)"
                strokeWidth={4}
                paintOrder="stroke"
                pointerEvents="none"
              >
                q {hoverQ.toFixed(4)} Å⁻¹ · d {((2 * Math.PI) / hoverQ).toFixed(1)} Å
              </text>
            )}
          </>
        );
      }}
    />
  );
}
