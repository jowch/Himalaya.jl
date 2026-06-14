import { useCallback, useEffect, useId, useRef, useState } from "react";
import { PlotFrame, type Margins, type PlotDims } from "./PlotFrame";
import {
  makeProjection,
  positiveExtent,
  type Projection,
  type ScaleType,
} from "./projection";
import { Axis } from "./Axis";
import { TraceLine } from "./marks/TraceLine";
import { PlotPeaks, type PlotPeak, type PeakFocusRequest } from "./marks/PlotPeaks";
import { PlotLabels } from "./marks/PlotLabels";
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
  onClickPeak?: (peakId: number) => void;
  onReset?: () => void;
  hitTolerancePx?: number;
}

export interface TracePlotProps {
  /** The single trace to render (hero / mini / one waterfall row). */
  trace: TraceModel;
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
  /** When non-empty, peaks/labels NOT in this set fade to neutral gray. Hot peaks are exempt. */
  highlightPeakIds?: ReadonlySet<number>;
  /** Multiply the y-domain top by (1 + yHeadroom) so peaks keep headroom below
   *  the ceiling — used by the stacked waterfall. Default 0 (no change). */
  yHeadroom?: number;
  /** Emitted when the USER hovers a peak (internal hover only — frame hit-test
   *  or glyph focus). The cross-panel q-link source. */
  onHoverQ?: (q: number | undefined) => void;
  /** Incoming hot q from another panel: light the peak matching this q (within
   *  the hit tolerance). The cross-panel q-link sink. */
  hoveredQ?: number;
  /** One-shot keyboard-focus re-anchor after a destructive peak edit (WCAG
   *  2.4.3). Forwarded to PlotPeaks. */
  focusRequest?: PeakFocusRequest;
  /** Called when a focusRequest has no surviving mark to land on. Forwarded to
   *  PlotPeaks; the consumer parks focus elsewhere (e.g. "+ Peak" button). */
  onFocusFallback?: () => void;
  /** Accessible name for the trace figure svg (WCAG 1.1.1). Without it the svg
   *  surfaces as a nameless img in the a11y tree. Forwarded to PlotFrame as
   *  role="img" + aria-label. */
  figureLabel?: string;
}

// BU-TRACE-DIM: unassigned (form-factor / unindexed) traces read at FULL visual
// strength at rest — a present neutral (ink-soft), NOT ink-faint, which looked
// like a stuck hover-dim next to the phase-colored traces. De-emphasis is
// reserved for an explicit hover/focus (the waterfall row's opacity gate), and
// clears completely on mouse-out.
const UNINDEXED_COLOR = "var(--color-ink-soft)";

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
    trace,
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
    highlightPeakIds,
    yHeadroom = 0,
    onHoverQ,
    hoveredQ,
    focusRequest,
    onFocusFallback,
    figureLabel,
  } = props;

  const layers = { peaks: true, labels: false, band: true, ...(show ?? {}) };

  const xExtent = positiveExtent(trace.trace.q);
  const rawYExtent = positiveExtent(trace.trace.I);
  const yExtent: [number, number] = [rawYExtent[0], rawYExtent[1] * (1 + yHeadroom)];
  const curXDomain = xDomain ?? xExtent;

  const margins: Margins = axes
    ? { top: 20, right: 24, bottom: 48, left: 60 }
    : { top: 4, right: 4, bottom: 4, left: 4 };

  const stateRef = useRef<PlotState | null>(null);

  const [hoverId, setHoverId] = useState<number | null>(null);
  // Unique clip-path id (multiple TracePlots coexist — e.g. the waterfall stack).
  // useId() returns colon-wrapped ids; strip them so url(#…) resolves.
  const clipId = `trace-clip-${useId().replace(/:/g, "")}`;

  // Outward q-link: emit the hovered peak's q whenever the INTERNAL hover changes
  // (covers both the frame hit-test and the per-glyph focus, which both write
  // hoverId). Keyed on hoverId — NOT hoveredQ — so an incoming external hover can
  // never round-trip back out and form a feedback loop.
  useEffect(() => {
    if (!onHoverQ) return;
    const p =
      hoverId == null ? undefined : trace.peaks.find((pk) => pk.id === hoverId);
    onHoverQ(p?.q);
  }, [hoverId, onHoverQ, trace.peaks]);

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
    (px: number, py: number) => {
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
        s.interaction.onClickPeak(hit.id);
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
      {...(figureLabel ? { svgRole: "img", svgLabel: figureLabel } : {})}
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
        stateRef.current = {
          projection,
          dims,
          curXDomain,
          xExtent,
          xType,
          interaction: interaction || null,
          peaks: trace.peaks,
        };
        const ctx: PlotContext = {
          projection,
          dims,
          hitTest: (ps, px, tol) =>
            hitTestPeaks(ps, px, (q) => projection.x.to(q), tol ?? PEAK_HIT_PX),
        };
        // Incoming q-link: resolve hoveredQ to the nearest peak within the same
        // pixel tolerance the hit-test uses. Internal hover stays authoritative.
        const externalHotId =
          hoveredQ == null
            ? null
            : (() => {
                let best: number | null = null;
                let bestPx = Infinity;
                const tol =
                  (interaction && interaction.hitTolerancePx) || PEAK_HIT_PX;
                const targetPx = projection.x.to(hoveredQ);
                for (const p of trace.peaks) {
                  const d = Math.abs(projection.x.to(p.q) - targetPx);
                  if (d <= tol && d < bestPx) {
                    bestPx = d;
                    best = p.id;
                  }
                }
                return best;
              })();
        const effectiveHoverId = hoverId ?? externalHotId;
        return (
          <>
            {/* Clip the trace + peaks + labels to the plot rect so a curve / glyphs
                whose q falls outside the visible window — common after a zoom — stop
                at the axes instead of drawing over the spines and tick labels. The
                top is left open to the SVG edge (y = -margins.top) so peak labels
                keep their designed headroom above the curve; only left / right /
                bottom (the three axis edges) bound the annotations. */}
            <defs>
              <clipPath id={clipId}>
                <rect
                  x={0}
                  y={-dims.margins.top}
                  width={dims.plotWidth}
                  height={dims.plotHeight + dims.margins.top}
                />
              </clipPath>
            </defs>
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
            <g clipPath={`url(#${clipId})`}>
              <TraceLine
                trace={trace.trace}
                projection={projection}
                band={layers.band}
                color={trace.phase ? phaseColor(trace.phase) : UNINDEXED_COLOR}
              />
              {layers.peaks ? (() => {
              const peaksWithHover = effectiveHoverId == null
                ? trace.peaks
                : trace.peaks.map((p) => (p.id === effectiveHoverId ? { ...p, hot: true } : p));
              return (
                <PlotPeaks
                  peaks={peaksWithHover}
                  projection={projection}
                  color={trace.phase ? phaseColor(trace.phase) : UNINDEXED_COLOR}
                  baselineI={yExtent[0]}
                  {...(paperColor ? { paperColor } : {})}
                  {...(interaction ? { onPeakFocus: setHoverId } : {})}
                  {
                    // Keyboard peak editing rides the same gate as pointer
                    // editing: onClickPeak is only threaded while the consumer
                    // is armed (TracePlate strips it otherwise), so disarmed
                    // marks stay roleless and unfocusable.
                    ...(interaction && interaction.onClickPeak
                      ? { onPeakActivate: interaction.onClickPeak }
                      : {})
                  }
                  {...(highlightPeakIds ? { highlightPeakIds } : {})}
                  {...(focusRequest ? { focusRequest } : {})}
                  {...(onFocusFallback ? { onFocusFallback } : {})}
                />
              );
            })() : null}
              {layers.labels ? (
                <PlotLabels
                  peaks={trace.peaks}
                  projection={projection}
                  color={trace.phase ? phaseColor(trace.phase) : UNINDEXED_COLOR}
                  baselineI={yExtent[0]}
                  {...(highlightPeakIds ? { highlightPeakIds } : {})}
                />
              ) : null}
            </g>
            {overlay ? overlay(ctx) : null}
            {/* C5: q-readout chip — anchored to axis bottom, on top of all other layers.
                Motion is instant (no CSS transition), which is inherently reduced-motion-safe.
                A 90ms opacity fade on enter/leave is intentionally deferred. */}
            {(() => {
              // Gate the readout on the peaks layer: with peaks hidden there is
              // no glyph to anchor it to, so a floating chip would be dishonest.
              const hovered = effectiveHoverId == null || !layers.peaks
                ? null
                : trace.peaks.find((p) => p.id === effectiveHoverId) ?? null;
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
                    {hovered.q.toFixed(3)}
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
