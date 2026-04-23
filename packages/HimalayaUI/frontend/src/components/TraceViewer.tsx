import { useEffect, useRef } from "react";
import * as Plot from "@observablehq/plot";
import type { Trace, Peak } from "../api";

function withIntensity(peaks: Peak[], trace: Trace): Array<{ q: number; I: number }> {
  return peaks.map((p) => ({ q: p.q, I: interpolateI(p.q, trace) }));
}

function interpolateI(q: number, trace: Trace): number {
  let nearest = 0;
  for (let i = 1; i < trace.q.length; i++) {
    if (Math.abs(trace.q[i]! - q) < Math.abs(trace.q[nearest]! - q)) nearest = i;
  }
  return trace.I[nearest]!;
}

export function snapToLocalMax(q: number[], I: number[], qClick: number, K = 3): number {
  let nearest = 0;
  for (let i = 1; i < q.length; i++) {
    if (Math.abs(q[i]! - qClick) < Math.abs(q[nearest]! - qClick)) nearest = i;
  }
  const lo = Math.max(0, nearest - K);
  const hi = Math.min(q.length - 1, nearest + K);
  let best = lo;
  for (let i = lo + 1; i <= hi; i++) {
    if (I[i]! > I[best]!) best = i;
  }
  return q[best]!;
}

export function findNearestPeak(
  peaks: Peak[], qClick: number, tolerance: number,
): Peak | null {
  let best: Peak | null = null;
  let bestDist = Infinity;
  for (const p of peaks) {
    const d = Math.abs(p.q - qClick);
    if (d < bestDist) { best = p; bestDist = d; }
  }
  return best && bestDist <= tolerance ? best : null;
}

export interface TraceViewerProps {
  trace: Trace;
  peaks: Peak[];
  onAddPeak: (q: number) => void;
  onRemovePeak: (peakId: number) => void;
}

export function TraceViewer({
  trace, peaks, onAddPeak, onRemovePeak,
}: TraceViewerProps): JSX.Element {
  const hostRef = useRef<HTMLDivElement>(null);

  useEffect(() => {
    const host = hostRef.current;
    if (!host) return;

    const data = trace.q.map((q, i) => ({
      q, I: trace.I[i]!,
      lo: Math.max(1e-12, trace.I[i]! - trace.sigma[i]!),
      hi: trace.I[i]! + trace.sigma[i]!,
    }));

    const el = Plot.plot({
      width: host.clientWidth || 400,
      height: host.clientHeight || 300,
      marginLeft: 50, marginBottom: 40,
      x: { type: "log", label: "q (nm⁻¹)" },
      y: { type: "log", label: "I (a.u.)" },
      marks: [
        Plot.areaY(data, { x: "q", y1: "lo", y2: "hi",
          fill: "var(--color-accent)", fillOpacity: 0.15 }),
        Plot.line(data, { x: "q", y: "I",
          stroke: "var(--color-fg)", strokeWidth: 1 }),
        Plot.dot(withIntensity(peaks.filter((p) => p.source === "auto"), trace),
          { x: "q", y: "I",
            fill: "var(--color-accent)", r: 4 }),
        Plot.dot(withIntensity(peaks.filter((p) => p.source === "manual"), trace),
          { x: "q", y: "I",
            stroke: "var(--color-warning)", strokeWidth: 2, fill: "none", r: 5 }),
      ],
    });

    host.replaceChildren(el);

    function handleClick(ev: Event): void {
      const me = ev as MouseEvent;
      const xScale = (el as unknown as { scale: (name: string) => { invert?: (v: number) => number } | undefined }).scale("x");
      if (!xScale?.invert) return;
      const rect = el.getBoundingClientRect();
      const qClick = xScale.invert(me.clientX - rect.left);

      // Tolerance: 2% of the click q (relative, log-friendly)
      const tolerance = Math.max(qClick * 0.02, 1e-6);
      const existing = findNearestPeak(peaks, qClick, tolerance);
      if (existing) onRemovePeak(existing.id);
      else          onAddPeak(snapToLocalMax(trace.q, trace.I, qClick));
    }
    (el as unknown as EventTarget).addEventListener("click", handleClick);

    return () => {
      (el as unknown as EventTarget).removeEventListener("click", handleClick);
      host.replaceChildren();
    };
  }, [trace, peaks, onAddPeak, onRemovePeak]);

  return <div ref={hostRef} className="w-full h-full" data-testid="trace-viewer" />;
}
