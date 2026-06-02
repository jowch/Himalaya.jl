import { useCallback, useRef, useState } from "react";
import { PlotFrame, type Margins, type PlotDims } from "./PlotFrame";
import {
  makeProjection,
  positiveExtent,
  type Projection,
  type ScaleType,
} from "./projection";
import { Axis } from "./Axis";
import { TraceLine } from "./marks/TraceLine";
import { PlotPeaks, type PlotPeak } from "./marks/PlotPeaks";
import { PlotLabels } from "./marks/PlotLabels";
import { hitTestPeaks, zoomXDomain, PEAK_HIT_PX } from "./interaction";
import { phaseColor } from "../../phases";
import { type Trace } from "../../api";
import { formatAxis } from "../../lib/plot/formatAxis";

export interface TraceModel {
  trace: Trace;
  peaks: PlotPeak[];
  phase: string | null;
}

export interface PlotContext {
  projection: Projection;
  dims: PlotDims;
  hitTest: (peaks: PlotPeak[], px: number, tol?: number) => PlotPeak | null;
}

export interface TracePlotInteraction {
  onXDomain: (d: [number, number] | null) => void;
  onAddPeak?: (q: number) => void;
  onClickPeak?: (peakId: number, altKey: boolean) => void;
  onReset?: () => void;
  hitTolerancePx?: number;
}

export interface TracePlotProps {
  /** 1 = hero/mini; >1 overlays in shared scales. */
  traces: TraceModel[];
  /** Visible x window; null = full data extent. */
  xDomain?: [number, number] | null;
  xType?: ScaleType;
  yType?: ScaleType;
  axes?: boolean;
  xLabel?: string;
  yLabel?: string;
  interaction?: TracePlotInteraction | false;
  overlay?: (ctx: PlotContext) => React.ReactNode;
  height?: number;
  width?: number;
  /** Paper colour for peak-glyph halos. */
  paperColor?: string;
  className?: string;
  "data-testid"?: string;
  /** Gate which annotation layers render. Axes/grid stay governed by `axes`. */
  show?: { peaks?: boolean; labels?: boolean; band?: boolean };
}

const UNINDEXED_COLOR = "var(--color-ink-faint)";

interface PlotState {
  projection: Projection;
  dims: PlotDims;
  curXDomain: [number, number];
  xExtent: [number, number];
  xType: ScaleType;
  interaction: TracePlotInteraction | null;
  peaks: PlotPeak[];
}

export function TracePlot(props: TracePlotProps): JSX.Element {
  const {
    traces,
    xDomain = null,
    xType = "log",
    yType = "log",
    axes = true,
    xLabel = "q (Å⁻¹)",
    yLabel = "intensity (a.u.)",
    interaction = false,
    overlay,
    height = 320,
    width,
    paperColor,
    className,
    "data-testid": testid,
    show,
  } = props;

  const layers = { peaks: true, labels: false, band: true, ...(show ?? {}) };

  const xExtent = positiveExtent(traces.flatMap((t) => t.trace.q));
  const yExtent = positiveExtent(traces.flatMap((t) => t.trace.I));
  const curXDomain = xDomain ?? xExtent;

  const margins: Margins = axes
    ? { top: 20, right: 24, bottom: 48, left: 60 }
    : { top: 4, right: 4, bottom: 4, left: 4 };

  const stateRef = useRef<PlotState | null>(null);

  const [hoverId, setHoverId] = useState<number | null>(null);

  const handlePointerMovePx = useCallback((px: number, py: number) => {
    const s = stateRef.current;
    if (!s || !s.interaction) return;
    const plotPx = px - s.dims.margins.left;
    const plotPy = py - s.dims.margins.top;
    if (plotPx < 0 || plotPx > s.dims.plotWidth || plotPy < 0 || plotPy > s.dims.plotHeight) {
      setHoverId(null);
      return;
    }
    const tol = s.interaction.hitTolerancePx ?? PEAK_HIT_PX;
    const hit = hitTestPeaks(s.peaks, plotPx, (q) => s.projection.x.to(q), tol);
    setHoverId(hit ? hit.id : null);
  }, []);

  const handlePointerLeave = useCallback(() => setHoverId(null), []);

  const handleWheelPx = useCallback((deltaY: number, px: number) => {
    const s = stateRef.current;
    if (!s || !s.interaction) return;
    const cursorQ = s.projection.x.invert(px - s.dims.margins.left);
    if (!Number.isFinite(cursorQ)) return;
    const next = zoomXDomain({
      cursorQ,
      deltaY,
      current: s.curXDomain,
      extent: s.xExtent,
      type: s.xType,
    });
    if (next) s.interaction.onXDomain(next);
  }, []);

  const handleClickPx = useCallback(
    (px: number, py: number, altKey: boolean) => {
      const s = stateRef.current;
      if (!s || !s.interaction) return;
      const plotPx = px - s.dims.margins.left;
      const plotPy = py - s.dims.margins.top;
      // Ignore clicks landing in the axis margins, not the plot body — parity
      // with PlotSurface's interior guard (no spurious peak from the gutters).
      if (
        plotPx < 0 ||
        plotPx > s.dims.plotWidth ||
        plotPy < 0 ||
        plotPy > s.dims.plotHeight
      ) {
        return;
      }
      const q = s.projection.x.invert(plotPx);
      // q must be a positive, finite value (linear invert can go negative in
      // the gutters; log invert stays positive). Matches PlotSurface:252.
      if (!Number.isFinite(q) || q <= 0) return;
      const tol = s.interaction.hitTolerancePx ?? PEAK_HIT_PX;
      const hit = hitTestPeaks(
        s.peaks,
        plotPx,
        (qq) => s.projection.x.to(qq),
        tol,
      );
      if (hit && s.interaction.onClickPeak) {
        s.interaction.onClickPeak(hit.id, altKey);
      } else if (s.interaction.onAddPeak) {
        s.interaction.onAddPeak(q);
      }
    },
    [],
  );

  const handleDblClick = useCallback(() => {
    const s = stateRef.current;
    if (!s || !s.interaction) return;
    if (s.interaction.onReset) s.interaction.onReset();
    else s.interaction.onXDomain(null);
  }, []);

  return (
    <PlotFrame
      height={height}
      margins={margins}
      {...(width !== undefined ? { width } : {})}
      {...(className ? { className } : {})}
      {...(testid ? { "data-testid": testid } : {})}
      {...(interaction
        ? {
            onWheelPx: handleWheelPx,
            onClickPx: handleClickPx,
            onDoubleClickPx: handleDblClick,
            onPointerMovePx: handlePointerMovePx,
            onPointerLeave: handlePointerLeave,
          }
        : {})}
      render={(dims) => {
        const projection = makeProjection({
          xDomain: curXDomain,
          yDomain: yExtent,
          plotWidth: dims.plotWidth,
          plotHeight: dims.plotHeight,
          xType,
          yType,
        });
        const allPeaks = traces.flatMap((t) => t.peaks);
        stateRef.current = {
          projection,
          dims,
          curXDomain,
          xExtent,
          xType,
          interaction: interaction || null,
          peaks: allPeaks,
        };
        const ctx: PlotContext = {
          projection,
          dims,
          hitTest: (ps, px, tol) =>
            hitTestPeaks(ps, px, (q) => projection.x.to(q), tol ?? PEAK_HIT_PX),
        };
        return (
          <>
            {axes ? (
              <>
                <Axis
                  axis={projection.x}
                  orientation="bottom"
                  plotWidth={dims.plotWidth}
                  plotHeight={dims.plotHeight}
                  label={xLabel}
                />
                <Axis
                  axis={projection.y}
                  orientation="left"
                  plotWidth={dims.plotWidth}
                  plotHeight={dims.plotHeight}
                  label={yLabel}
                />
              </>
            ) : null}
            {traces.map((t, i) => (
              <TraceLine key={`line-${i}`} trace={t.trace} projection={projection} band={layers.band} />
            ))}
            {layers.peaks ? traces.map((t, i) => {
              // C4: inject hot flag for the hovered peak (hover does NOT dim others)
              const peaksWithHover = hoverId == null
                ? t.peaks
                : t.peaks.map((p) => (p.id === hoverId ? { ...p, hot: true } : p));
              return (
                <PlotPeaks
                  key={`peaks-${i}`}
                  peaks={peaksWithHover}
                  projection={projection}
                  color={t.phase ? phaseColor(t.phase) : UNINDEXED_COLOR}
                  baselineI={yExtent[0]}
                  {...(paperColor ? { paperColor } : {})}
                  {...(interaction ? { onPeakFocus: setHoverId } : {})}
                />
              );
            }) : null}
            {layers.labels ? traces.map((t, i) => (
              <PlotLabels
                key={`labels-${i}`}
                peaks={t.peaks}
                projection={projection}
                color={t.phase ? phaseColor(t.phase) : UNINDEXED_COLOR}
                baselineI={yExtent[0]}
              />
            )) : null}
            {overlay ? overlay(ctx) : null}
            {/* C5: q-readout chip — anchored to axis bottom, on top of all other layers.
                Motion is instant (no CSS transition), which is inherently reduced-motion-safe.
                A 90ms opacity fade on enter/leave is intentionally deferred. */}
            {(() => {
              const hovered = hoverId == null
                ? null
                : traces.flatMap((t) => t.peaks).find((p) => p.id === hoverId) ?? null;
              if (!hovered) return null;
              const qx = projection.x.to(hovered.q);
              const baseY = dims.plotHeight;
              const w = 46, h = 16;
              return (
                <g data-role="q-readout" transform={`translate(${qx},${baseY})`}>
                  <rect x={-w / 2} y={6} width={w} height={h} rx={3}
                    fill="var(--color-plate)" stroke="var(--color-hair-strong)" strokeWidth={1} />
                  <text x={0} y={6 + h / 2} dy="0.32em" textAnchor="middle"
                    style={{ fontFamily: "var(--font-mono)", fontSize: 10.5, fill: "var(--color-ink)" }}>
                    {formatAxis(hovered.q)}
                  </text>
                </g>
              );
            })()}
          </>
        );
      }}
    />
  );
}
