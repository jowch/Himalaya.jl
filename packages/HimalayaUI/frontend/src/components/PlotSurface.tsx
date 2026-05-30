/**
 * PlotSurface — the shared Observable-Plot host (Plan C plot spine).
 *
 * Owns the single `Plot.plot({...})` instance, a shared log/linear x-scale,
 * the gesture layer (wheel-zoom-about-cursor, click→hit-test/add, dblclick
 * reset), rAF-friendly resize (ResizeObserver → key bump; clientWidth read in
 * render), and exposes an `overlay(ctx)` render slot + a `hitTest` closure so
 * consumers (TraceViewer, MultiTracePlot, and Plans D/E) attach interaction
 * WITHOUT re-implementing the scale plumbing.
 *
 * Margins + xType stay per-consumer props — they genuinely differ across the
 * Focus hero (50/14/36/40) and Compare waterfall (50/14/8/32).
 *
 * NOTE: MillerPlot is intentionally NOT a consumer — it's a q-vs-ratio scatter
 * with regression overlays, not a q-axis plot, and does not fit this contract.
 */
import { useCallback, useEffect, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import { invertQ as invertQHelper, applyQ as applyQHelper } from "../lib/plot/invertQ";
import { formatAxis } from "../lib/plot/formatAxis";
import { prettifyUnits } from "../lib/units";
import type { Peak } from "../api";

/** Click within this many pixels of a peak to act on it. Preserved from
 *  TraceViewer/MultiTracePlot (#180 q-link feel). */
export const PEAK_HIT_PX = 10;

type Scale =
  | { invert?: (v: number) => number; apply?: (v: number) => number }
  | undefined;

export interface PlotMargins {
  left: number;
  right: number;
  top: number;
  bottom: number;
}

export interface PlotOverlayContext {
  /** q → px (null when the scale isn't ready). */
  applyQ: (q: number) => number | null;
  /** px → q (null when the scale isn't ready). */
  invertQ: (px: number) => number | null;
  /** data-value → px on the y axis. */
  applyY: (v: number) => number | null;
  width: number;
  height: number;
  margins: PlotMargins;
  /** Nearest clickable (id>=0) peak to a pixel x, or null. Skips optimistic. */
  hitTest: (peaks: Peak[], px: number, tolerancePx?: number) => Peak | null;
}

export interface PlotSurfaceProps {
  /** Plot marks for the data layer (trace line, band, member rows, heatmap). */
  marks: Plot.Markish[];
  xType: "log" | "linear";
  xDomain: [number, number] | null;
  yDomain?: [number, number] | null;
  /** Y scale type. Focus hero is log-I; Compare bands are linear pixel-space. */
  yType?: "log" | "linear";
  /** Reverse the y domain (Compare bands: small-y at top). */
  yReversed?: boolean;
  /** Hide the y axis (Compare: bands ARE the visualization). */
  hideYAxis?: boolean;
  yLabel?: string;
  margins: PlotMargins;
  qUnits?: string;
  onXDomain: (d: [number, number] | null) => void;
  /** Full extent used to clamp wheel-zoom. Defaults to xDomain. */
  xExtent?: [number, number] | null;
  /** Peaks for click hit-testing. Omit on read-only surfaces. */
  peaks?: Peak[];
  onAddPeak?: (q: number) => void;
  onClickPeak?: (peakId: number, altKey: boolean) => void;
  onReset?: () => void;
  /** Render the interactive SVG overlay given live scales. */
  overlay?: (ctx: PlotOverlayContext) => React.ReactNode | void;
  hitTolerancePx?: number;
  /** Forwarded to the root for E2E / placement. */
  className?: string;
  "data-testid"?: string;
}

/** Nearest clickable (id>=0) peak to a pixel x. Skips optimistic placeholders
 *  (id<0): routing a click to one would hit a 404 on the not-yet-created row.
 *  Shared hit-test for the surface + the overlay context. */
export function hitTestPeaks(
  peaks: Peak[],
  clickX: number,
  toPx: (q: number) => number,
  tolerancePx: number,
): Peak | null {
  let best: Peak | null = null;
  let bestDist = tolerancePx;
  for (const p of peaks) {
    if (p.id < 0) continue;
    const px = toPx(p.q);
    if (!Number.isFinite(px)) continue;
    const d = Math.abs(px - clickX);
    if (d <= bestDist) {
      best = p;
      bestDist = d;
    }
  }
  return best;
}

export function PlotSurface(props: PlotSurfaceProps): JSX.Element {
  const {
    marks, xType, xDomain, yDomain = null,
    yType = "log", yReversed = false, hideYAxis = false, yLabel = "I (a.u.)",
    margins, qUnits, onXDomain, xExtent = null,
    peaks = [], onAddPeak, onClickPeak, onReset, overlay,
    hitTolerancePx = PEAK_HIT_PX,
    className, "data-testid": dataTestId = "plot-surface",
  } = props;

  const hostRef = useRef<HTMLDivElement>(null);
  const plotContainer = useRef<HTMLDivElement>(null);
  const plotElRef = useRef<HTMLElement | SVGElement | null>(null);

  const [_resizeKey, setResizeKey] = useState(0);
  // The overlay subtree reads scales off the live plot element, which only
  // exists AFTER the first effect runs. `ready` flips false→true exactly once
  // so we re-render to paint the overlay; thereafter `_resizeKey` re-renders
  // the overlay against fresh scales (it's also a renderPlot effect trigger).
  const [ready, setReady] = useState(false);

  useEffect(() => {
    const el = plotContainer.current;
    if (!el) return;
    const obs = new ResizeObserver(() => setResizeKey((k) => k + 1));
    obs.observe(el);
    return () => obs.disconnect();
  }, []);

  // Imperative render-and-bind. useCallback with true deps so the effect below
  // depends on [renderPlot, _resizeKey] alone — handlers never close over stale
  // onXDomain/onAddPeak (CLAUDE.md "imperative renderers in effects").
  const renderPlot = useCallback((): (() => void) | undefined => {
    const container = plotContainer.current;
    if (!container) return;

    const width = container.clientWidth || 400;
    const height = container.clientHeight || 300;

    const el = Plot.plot({
      width,
      height,
      marginLeft: margins.left,
      marginRight: margins.right,
      marginTop: margins.top,
      marginBottom: margins.bottom,
      style: {
        fontFamily: "var(--font-sans)",
        color: "var(--color-ink-soft)",
        background: "transparent",
        overflow: "visible",
      },
      x: {
        type: xType,
        label: `q (${prettifyUnits(qUnits ?? "A-1")})`,
        tickFormat: (d: number) => formatAxis(d),
        ...(xDomain ? { domain: xDomain } : {}),
      },
      y: {
        type: yType,
        label: yLabel,
        tickFormat: (d: number) => formatAxis(d),
        ...(hideYAxis ? { axis: null } : {}),
        ...(yReversed && yDomain
          ? { domain: [yDomain[1], yDomain[0]] }
          : yDomain
            ? { domain: yDomain }
            : {}),
      },
      marks,
    });

    container.replaceChildren(el);
    plotElRef.current = el as unknown as HTMLElement;

    // ── click: hit-test → onClickPeak, else onAddPeak(invertQ) ─────────────
    function handleClick(ev: Event): void {
      const me = ev as MouseEvent;
      const xScale: Scale = (
        plotElRef.current as unknown as { scale: (n: string) => Scale }
      )?.scale("x");
      if (!xScale?.invert || !xScale.apply) return;
      const rect = container!.getBoundingClientRect();
      const clickX = me.clientX - rect.left;
      const clickY = me.clientY - rect.top;
      if (!insideInterior(clickX, clickY, rect.width, rect.height, margins)) return;

      if (onClickPeak) {
        const hit = hitTestPeaks(peaks, clickX, xScale.apply, hitTolerancePx);
        if (hit) {
          onClickPeak(hit.id, me.altKey);
          return;
        }
      }
      if (onAddPeak) {
        const q = xScale.invert(clickX);
        if (Number.isFinite(q) && q > 0) onAddPeak(q);
      }
    }
    (el as unknown as EventTarget).addEventListener("click", handleClick);

    // ── wheel: zoom x-domain about the cursor (log-aware) ──────────────────
    function handleWheel(evRaw: Event): void {
      const ev = evRaw as WheelEvent;
      ev.preventDefault();
      const rect = container!.getBoundingClientRect();
      const cursorQ = invertQHelper(plotElRef.current, ev.clientX - rect.left);
      if (cursorQ === null) return;
      const ext = xExtent ?? xDomain;
      const curMin = xDomain ? xDomain[0] : ext ? ext[0] : cursorQ * 0.5;
      const curMax = xDomain ? xDomain[1] : ext ? ext[1] : cursorQ * 2;
      const factor = Math.exp(ev.deltaY * 0.001);
      const q0 = ext ? ext[0] : curMin;
      const qN = ext ? ext[1] : curMax;
      let newMin: number;
      let newMax: number;
      if (xType === "log") {
        const logMin = Math.log(curMin);
        const logMax = Math.log(curMax);
        const logCur = Math.log(Math.max(cursorQ, 1e-6));
        newMin = Math.max(q0, Math.exp(logCur - (logCur - logMin) * factor));
        newMax = Math.min(qN, Math.exp(logCur + (logMax - logCur) * factor));
      } else {
        newMin = Math.max(q0, cursorQ - (cursorQ - curMin) * factor);
        newMax = Math.min(qN, cursorQ + (curMax - cursorQ) * factor);
      }
      if (newMax - newMin < (qN - q0) * 1e-4) return;
      onXDomain([newMin, newMax]);
    }
    (el as unknown as EventTarget).addEventListener("wheel", handleWheel, {
      passive: false,
    } as AddEventListenerOptions);

    // ── dblclick: reset domain (delegates to onReset if provided) ──────────
    function handleDblClick(): void {
      if (onReset) onReset();
      else onXDomain(null);
    }
    (el as unknown as EventTarget).addEventListener("dblclick", handleDblClick);

    // Paint the overlay once the scales exist (one-shot; subsequent rebuilds
    // are picked up via the `_resizeKey` re-render below).
    setReady(true);

    return () => {
      (el as unknown as EventTarget).removeEventListener("click", handleClick);
      (el as unknown as EventTarget).removeEventListener("wheel", handleWheel);
      (el as unknown as EventTarget).removeEventListener("dblclick", handleDblClick);
      container.replaceChildren();
      plotElRef.current = null;
    };
  }, [
    marks, xType, xDomain, yDomain, yType, yReversed, hideYAxis, yLabel,
    margins, qUnits, onXDomain, xExtent, peaks, onAddPeak, onClickPeak,
    onReset, hitTolerancePx,
  ]);

  useEffect(() => {
    return renderPlot();
  }, [renderPlot, _resizeKey]);

  // Build the overlay context from the live plot element + scales.
  const overlayContent = ((): React.ReactNode => {
    if (!overlay) return null;
    void _resizeKey; // re-derive context on resize
    void ready;      // and once the plot element first exists
    const plotEl = plotElRef.current;
    if (!plotEl) return null;
    const yScale: Scale = (
      plotEl as unknown as { scale: (n: string) => Scale }
    ).scale("y");
    const container = plotContainer.current;
    const width = container?.clientWidth ?? 0;
    const height = container?.clientHeight ?? 0;
    const ctx: PlotOverlayContext = {
      applyQ: (q) => applyQHelper(plotEl, q),
      invertQ: (px) => invertQHelper(plotEl, px),
      applyY: (v) => (yScale?.apply ? yScale.apply(v) : null),
      width,
      height,
      margins,
      hitTest: (ps, px, tol) =>
        hitTestPeaks(ps, px, (q) => applyQHelper(plotEl, q) ?? NaN, tol ?? hitTolerancePx),
    };
    return (overlay(ctx) ?? null) as React.ReactNode;
  })();

  return (
    <div
      ref={hostRef}
      className={["w-full h-full relative", className].filter(Boolean).join(" ")}
      data-testid={dataTestId}
    >
      <div ref={plotContainer} className="w-full h-full" />
      {overlay && (
        <svg
          className="absolute inset-0 pointer-events-none"
          aria-hidden="true"
          data-role="plot-surface-overlay"
        >
          {overlayContent}
        </svg>
      )}
    </div>
  );
}

function insideInterior(
  x: number,
  y: number,
  w: number,
  h: number,
  m: PlotMargins,
): boolean {
  return x >= m.left && x <= w - m.right && y >= m.top && y <= h - m.bottom;
}
