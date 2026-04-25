import { useCallback, useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import type { Trace, Peak, IndexEntry } from "../api";
import { phaseColor } from "../phases";

export interface TraceViewerProps {
  trace: Trace;
  peaks: Peak[];
  activeGroupIndices: IndexEntry[];
  hoveredIndex: IndexEntry | undefined;
  onAddPeak: (q: number) => void;
  onRemovePeak: (peakId: number) => void;
  onTogglePeakExclusion: (peakId: number, excluded: boolean) => void;
  /** Visible q-range. null = auto (full trace). */
  xDomain: [number, number] | null;
  onXDomain: (d: [number, number] | null) => void;
}

// ── constants ──────────────────────────────────────────────────────────────

/** Click within this many pixels of a peak triangle to act on it. */
const PEAK_HIT_PX = 10;

/** Pixel offset of a peak triangle above the trace. ~half of the previous log-space lift. */
const PEAK_OFFSET_PX = 7;

/** Plot margins (kept in sync with the effect so overlay can clip / position). */
const MARGIN_LEFT   = 50;
const MARGIN_RIGHT  = 14;
const MARGIN_TOP    = 22;
const MARGIN_BOTTOM = 40;

/** Triangle marker geometry. */
const TRIANGLE_HALF_W = 4;
const TRIANGLE_H      = 7;

/** When a predicted-q tick has no nearby detected peak, terminate this many pixels above the trace. */
const TICK_END_OFFSET_PX = 14;

/** Pixel tolerance to consider a predicted q "matched" by a detected peak. */
const TICK_MATCH_PX = 12;

// ── helpers ────────────────────────────────────────────────────────────────

export function findNearestPeak(peaks: Peak[], qClick: number, tolerance: number): Peak | null {
  let best: Peak | null = null;
  let bestDist = Infinity;
  for (const p of peaks) {
    const d = Math.abs(p.q - qClick);
    if (d < bestDist) { best = p; bestDist = d; }
  }
  return best && bestDist <= tolerance ? best : null;
}

/**
 * Plain-decimal axis label formatter for SAXS log scales. Avoids Plot's
 * default SI-suffix formatter (which renders 0.04 as "40m"). Uses scientific
 * notation only at the extremes where decimals would be unreadable.
 */
function formatAxis(d: number): string {
  const ad = Math.abs(d);
  if (ad === 0) return "0";
  if (ad < 1e-3 || ad >= 1e4) return d.toExponential(0).replace("e+", "e").replace("e-0", "e-").replace("e0", "e");
  if (ad < 0.01) return d.toFixed(3);
  if (ad < 1)    return d.toFixed(2);
  if (ad < 100)  return d.toFixed(ad < 10 ? 1 : 0);
  return d.toFixed(0);
}

function interpolateI(q: number, trace: Trace): number {
  const qs = trace.q;
  let nearest = 0;
  for (let i = 1; i < qs.length; i++) {
    if (Math.abs(qs[i]! - q) < Math.abs(qs[nearest]! - q)) nearest = i;
  }
  return trace.I[nearest]!;
}

type Scale = { invert?: (v: number) => number; apply?: (v: number) => number } | undefined;

interface IndexTick { q: number; color: string; indexId: number }
function indexTicks(indices: IndexEntry[]): IndexTick[] {
  return indices.flatMap((ix) =>
    ix.predicted_q.map((q) => ({ q, color: phaseColor(ix.phase), indexId: ix.id })));
}

// ── component ──────────────────────────────────────────────────────────────

export function TraceViewer({
  trace, peaks, activeGroupIndices, hoveredIndex,
  onAddPeak, onRemovePeak, onTogglePeakExclusion,
  xDomain, onXDomain,
}: TraceViewerProps): JSX.Element {
  const hostRef       = useRef<HTMLDivElement>(null);
  const plotContainer = useRef<HTMLDivElement>(null);
  const overlayRef    = useRef<SVGSVGElement>(null);
  const plotElRef     = useRef<HTMLElement | SVGElement | null>(null);

  // Re-render trigger for the overlay (needed because our deps include peaks
  // and indices; effects close over those values).
  useEffect(() => {
    const host      = hostRef.current;
    const container = plotContainer.current;
    if (!host || !container) return;

    const bandData = trace.q.map((q, i) => ({
      q, I: trace.I[i]!,
      lo: Math.max(1e-12, trace.I[i]! - trace.sigma[i]!),
      hi: trace.I[i]! + trace.sigma[i]!,
    }));

    const el = Plot.plot({
      width:  container.clientWidth  || 400,
      height: container.clientHeight || 300,
      marginLeft: MARGIN_LEFT, marginRight: MARGIN_RIGHT,
      marginTop: MARGIN_TOP,  marginBottom: MARGIN_BOTTOM,
      style: {
        fontFamily: "var(--font-sans)",
        color: "var(--color-fg-muted)",
        background: "transparent",
        overflow: "visible",
      },
      x: {
        type: "log",
        label: "q (Å⁻¹)",
        // Plain decimal tick labels — Plot's default SI-suffix formatter
        // renders 0.040 as "40m" which is unhelpful for SAXS q values.
        tickFormat: (d: number) => formatAxis(d),
        ...(xDomain ? { domain: xDomain } : {}),
      },
      y: {
        type: "log",
        label: "I (a.u.)",
        tickFormat: (d: number) => formatAxis(d),
      },
      // The Plot only renders the data trace + sigma band. Everything else
      // (peak triangles, predicted-q lines, cursor) lives in the overlay so
      // we have full control over geometry, fade-on-hover, and cursor follow.
      marks: [
        Plot.areaY(bandData, {
          x: "q", y1: "lo", y2: "hi",
          fill: "var(--color-accent)", fillOpacity: 0.12,
        }),
        Plot.line(bandData, {
          x: "q", y: "I",
          stroke: "var(--color-fg)", strokeWidth: 1,
        }),
      ],
    });

    container.replaceChildren(el);
    plotElRef.current = el as unknown as HTMLElement;

    // ── click: add / remove / toggle-exclude based on what's near the cursor ─
    function handleClick(ev: Event): void {
      const me = ev as MouseEvent;
      const xScale: Scale = (plotElRef.current as unknown as { scale: (n: string) => Scale })?.scale("x");
      if (!xScale?.invert || !xScale.apply) return;
      const rect = container!.getBoundingClientRect();
      const clickX = me.clientX - rect.left;
      const clickY = me.clientY - rect.top;
      if (!insideInterior(clickX, clickY, rect.width, rect.height)) return;

      // Pixel-proximity: walk visible peaks and check which the user hit.
      let bestPeak: Peak | null = null;
      let bestDist = PEAK_HIT_PX;
      for (const p of peaks) {
        const px = xScale.apply!(p.q);
        if (!Number.isFinite(px)) continue;
        const d = Math.abs(px - clickX);
        if (d <= bestDist) { bestPeak = p; bestDist = d; }
      }

      if (bestPeak) {
        if (bestPeak.source === "manual") onRemovePeak(bestPeak.id);
        else                              onTogglePeakExclusion(bestPeak.id, !bestPeak.excluded);
        return;
      }

      // Empty area → add a manual peak at the exact clicked q (no snap).
      const q = xScale.invert(clickX);
      if (Number.isFinite(q) && q > 0) onAddPeak(q);
    }
    (el as unknown as EventTarget).addEventListener("click", handleClick);

    // ── wheel: zoom x-domain around cursor ───────────────────────────────
    function handleWheel(evRaw: Event): void {
      const ev = evRaw as WheelEvent;
      ev.preventDefault();
      const xScale: Scale = (plotElRef.current as unknown as { scale: (n: string) => Scale })?.scale("x");
      if (!xScale?.invert) return;
      const rect = container!.getBoundingClientRect();
      const cursorQ = xScale.invert(ev.clientX - rect.left);
      const curMin  = xDomain ? xDomain[0] : trace.q[0]!;
      const curMax  = xDomain ? xDomain[1] : trace.q[trace.q.length - 1]!;
      const factor  = Math.exp(ev.deltaY * 0.001);
      const logMin  = Math.log(curMin);
      const logMax  = Math.log(curMax);
      const logCur  = Math.log(Math.max(cursorQ, 1e-6));
      const newLogMin = logCur - (logCur - logMin) * factor;
      const newLogMax = logCur + (logMax - logCur) * factor;
      const q0 = trace.q[0]!;
      const qN = trace.q[trace.q.length - 1]!;
      const newMin = Math.max(q0, Math.exp(newLogMin));
      const newMax = Math.min(qN, Math.exp(newLogMax));
      if (newMax - newMin < (qN - q0) * 1e-4) return;
      onXDomain([newMin, newMax]);
    }
    (el as unknown as EventTarget).addEventListener("wheel", handleWheel, { passive: false } as AddEventListenerOptions);

    // ── dblclick: reset x-domain ────────────────────────────────────────
    function handleDblClick(): void { onXDomain(null); }
    (el as unknown as EventTarget).addEventListener("dblclick", handleDblClick);

    renderOverlay();

    return () => {
      (el as unknown as EventTarget).removeEventListener("click", handleClick);
      (el as unknown as EventTarget).removeEventListener("wheel", handleWheel);
      (el as unknown as EventTarget).removeEventListener("dblclick", handleDblClick);
      container.replaceChildren();
      plotElRef.current = null;
    };
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [trace, peaks, activeGroupIndices, hoveredIndex, xDomain, onAddPeak, onRemovePeak, onTogglePeakExclusion, onXDomain]);

  // ── overlay renderer (peaks + predicted-q lines + cursor) ───────────────
  const renderOverlay = useCallback((): void => {
    const host    = hostRef.current;
    const overlay = overlayRef.current;
    if (!host || !overlay) return;
    const plotEl  = plotElRef.current as unknown as SVGElement | null;
    if (!plotEl) return;
    const xScale: Scale = (plotEl as unknown as { scale: (n: string) => Scale }).scale("x");
    const yScale: Scale = (plotEl as unknown as { scale: (n: string) => Scale }).scale("y");
    if (!xScale?.apply || !yScale?.apply) return;

    const bbox = (plotEl as Element).getBoundingClientRect();
    overlay.setAttribute("width",  String(bbox.width));
    overlay.setAttribute("height", String(bbox.height));

    // Wipe and redraw the peak + tick layers.
    const peakRoot = overlay.querySelector<SVGGElement>("[data-role=peak-root]")!;
    const tickRoot = overlay.querySelector<SVGGElement>("[data-role=tick-root]")!;
    while (peakRoot.firstChild) peakRoot.removeChild(peakRoot.firstChild);
    while (tickRoot.firstChild) tickRoot.removeChild(tickRoot.firstChild);

    const dimOthers = hoveredIndex !== undefined;
    const insideX = (px: number): boolean =>
      bbox.width === 0 || (px >= MARGIN_LEFT && px <= bbox.width - MARGIN_RIGHT);

    // ── 1. Peak triangles. Compute pixel positions; cache for tick lookup ──
    interface PeakDraw { peak: Peak; px: number; py: number }
    const peakDraws: PeakDraw[] = [];
    for (const p of peaks) {
      const px = xScale.apply!(p.q);
      if (!Number.isFinite(px) || !insideX(px)) continue;
      const I = interpolateI(p.q, trace);
      const py = yScale.apply!(I) - PEAK_OFFSET_PX;
      if (!Number.isFinite(py)) continue;
      peakDraws.push({ peak: p, px, py });
    }

    for (const { peak, px, py } of peakDraws) {
      const isAuto    = peak.source === "auto";
      // Bright/neon for "active workflow": auto = ice blue, manual = magenta.
      const baseColor = isAuto ? "var(--color-accent)" : "var(--color-peak-manual)";
      // Auto peaks: filled triangle. Excluded auto peaks: same color but ~30% opacity.
      // Manual peaks: filled magenta triangle (always full opacity when not faded).
      let fill: string;
      let opacity: number;
      const excludedAuto = isAuto && peak.excluded;
      if (excludedAuto) {
        // Excluded auto peaks keep their identity (ice blue, ghosted) — they
        // are user curation, not "context to fade away."
        fill = baseColor;
        opacity = 0.3;
      } else if (dimOthers) {
        // Hovering an index: peaks are not the focus → desaturate to gray
        // rather than just thin the alpha. Removing the color signal entirely
        // is what makes the hovered phase pop.
        fill = "var(--color-fg-dim)";
        opacity = 0.5;
      } else {
        fill = baseColor;
        opacity = isAuto ? 0.95 : 1;
      }

      const tri = document.createElementNS("http://www.w3.org/2000/svg", "polygon");
      tri.setAttribute("points",
        `${px - TRIANGLE_HALF_W},${py - TRIANGLE_H} ` +
        `${px + TRIANGLE_HALF_W},${py - TRIANGLE_H} ` +
        `${px},${py}`);
      tri.setAttribute("fill", fill);
      tri.setAttribute("fill-opacity", String(opacity));
      tri.setAttribute("stroke", "var(--color-bg)");
      tri.setAttribute("stroke-width", "0.75");
      peakRoot.appendChild(tri);
    }

    // ── 2. Predicted-q vlines. Adaptive endpoint: terminate above a matched
    //       peak triangle if any, else just above the trace. Defaults are
    //       deliberately quiet (thin, semi-transparent) so they don't
    //       compete with the data; on hover the active phase pops to the
    //       previous emphasised stroke and the others gray out entirely.
    function drawTickLine(t: IndexTick, opts: { strong: boolean; faded: boolean }): void {
      const px = xScale!.apply!(t.q);
      if (!Number.isFinite(px) || !insideX(px)) return;

      // Find a matching peak triangle within TICK_MATCH_PX.
      let matchedPy: number | null = null;
      for (const draw of peakDraws) {
        if (Math.abs(draw.px - px) <= TICK_MATCH_PX) {
          matchedPy = draw.py - TRIANGLE_H - 7;
          break;
        }
      }
      const traceY = yScale!.apply!(interpolateI(t.q, trace)) - TICK_END_OFFSET_PX;
      const y2     = matchedPy ?? traceY;

      // Faded → neutral gray (color signal removed entirely).
      const stroke = opts.faded ? "var(--color-fg-dim)" : t.color;
      const strokeWidth   = opts.strong ? "1.5" : "1";
      const strokeOpacity = opts.faded ? "0.3" : (opts.strong ? "1" : "0.35");

      const line = document.createElementNS("http://www.w3.org/2000/svg", "line");
      line.setAttribute("x1", String(px));
      line.setAttribute("x2", String(px));
      line.setAttribute("y1", String(MARGIN_TOP - 2));
      line.setAttribute("y2", String(y2));
      line.setAttribute("stroke", stroke);
      line.setAttribute("stroke-width", strokeWidth);
      line.setAttribute("stroke-opacity", strokeOpacity);
      line.setAttribute("stroke-linecap", "round");
      tickRoot.appendChild(line);
    }

    for (const t of indexTicks(activeGroupIndices)) {
      drawTickLine(t, { strong: false, faded: dimOthers });
    }
    if (hoveredIndex) {
      for (const t of indexTicks([hoveredIndex])) {
        drawTickLine(t, { strong: true, faded: false });
      }
    }
  }, [peaks, trace, hoveredIndex, activeGroupIndices, xDomain]);

  // Re-render overlay whenever anything that affects it changes.
  useEffect(() => {
    renderOverlay();
  }, [renderOverlay]);

  // ── cursor crosshair (separate, doesn't rebuild the plot) ──────────────
  useEffect(() => {
    const host    = hostRef.current;
    const overlay = overlayRef.current;
    if (!host || !overlay) return;

    let rafId = 0;
    function drawCursor(mx: number, my: number): void {
      const plotEl = plotElRef.current as unknown as SVGElement | null;
      if (!plotEl) return;
      const xScale: Scale = (plotEl as unknown as { scale: (n: string) => Scale }).scale("x");
      const yScale: Scale = (plotEl as unknown as { scale: (n: string) => Scale }).scale("y");
      if (!xScale?.invert || !yScale?.apply) return;

      const bbox = (plotEl as Element).getBoundingClientRect();
      const relX = mx - bbox.left;
      const relY = my - bbox.top;

      const line = overlay!.querySelector<SVGLineElement>("[data-role=cursor-line]")!;
      const dot  = overlay!.querySelector<SVGCircleElement>("[data-role=cursor-dot]")!;

      if (!insideInterior(relX, relY, bbox.width, bbox.height)) {
        line.setAttribute("opacity", "0");
        dot.setAttribute("opacity", "0");
        return;
      }

      const q  = xScale.invert(relX);
      const Iv = interpolateI(q, trace);
      const py = yScale.apply(Iv);

      line.setAttribute("x1", String(relX));
      line.setAttribute("x2", String(relX));
      line.setAttribute("y1", String(MARGIN_TOP));
      line.setAttribute("y2", String(bbox.height - MARGIN_BOTTOM));
      line.setAttribute("opacity", "1");
      dot.setAttribute("cx", String(relX));
      dot.setAttribute("cy", String(py));
      dot.setAttribute("opacity", Number.isFinite(py) ? "1" : "0");
    }

    function onMove(ev: MouseEvent): void {
      const mx = ev.clientX, my = ev.clientY;
      if (rafId) cancelAnimationFrame(rafId);
      rafId = requestAnimationFrame(() => drawCursor(mx, my));
    }
    function onLeave(): void {
      if (rafId) cancelAnimationFrame(rafId);
      const line = overlay!.querySelector<SVGLineElement>("[data-role=cursor-line]")!;
      const dot  = overlay!.querySelector<SVGCircleElement>("[data-role=cursor-dot]")!;
      line.setAttribute("opacity", "0");
      dot.setAttribute("opacity", "0");
    }

    host.addEventListener("mousemove", onMove);
    host.addEventListener("mouseleave", onLeave);
    return () => {
      host.removeEventListener("mousemove", onMove);
      host.removeEventListener("mouseleave", onLeave);
      if (rafId) cancelAnimationFrame(rafId);
    };
  }, [trace]);

  return (
    <div ref={hostRef} className="w-full h-full relative anim-overlay" data-testid="trace-viewer">
      <div ref={plotContainer} className="w-full h-full" />
      <svg
        ref={overlayRef}
        className="absolute inset-0 pointer-events-none"
        aria-hidden="true"
      >
        <g data-role="tick-root" />
        <g data-role="peak-root" />
        <line
          data-role="cursor-line"
          stroke="var(--color-fg-dim)"
          strokeWidth={1}
          strokeOpacity={0.7}
          opacity={0}
        />
        <circle
          data-role="cursor-dot"
          r={3.5}
          fill="none"
          stroke="var(--color-fg-dim)"
          strokeWidth={1.5}
          opacity={0}
        />
      </svg>
    </div>
  );
}

function insideInterior(x: number, y: number, w: number, h: number): boolean {
  return (
    x >= MARGIN_LEFT && x <= w - MARGIN_RIGHT &&
    y >= MARGIN_TOP  && y <= h - MARGIN_BOTTOM
  );
}
