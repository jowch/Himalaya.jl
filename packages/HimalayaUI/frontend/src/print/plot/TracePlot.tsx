import { useCallback, useRef } from "react";
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
import { hitTestPeaks, zoomXDomain, PEAK_HIT_PX } from "./interaction";
import { phaseColor } from "../../phases";
import { type Trace } from "../../api";

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
            {layers.peaks ? traces.map((t, i) => (
              <PlotPeaks
                key={`peaks-${i}`}
                peaks={t.peaks}
                projection={projection}
                color={t.phase ? phaseColor(t.phase) : UNINDEXED_COLOR}
                baselineI={yExtent[0]}
                {...(paperColor ? { paperColor } : {})}
              />
            )) : null}
            {overlay ? overlay(ctx) : null}
          </>
        );
      }}
    />
  );
}
