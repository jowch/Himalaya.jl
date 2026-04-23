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

export interface TraceViewerProps {
  trace: Trace;
  peaks: Peak[];
  onAddPeak: (q: number) => void;
  onRemovePeak: (peakId: number) => void;
}

export function TraceViewer({
  trace, peaks, onAddPeak: _onAddPeak, onRemovePeak: _onRemovePeak,
}: TraceViewerProps): JSX.Element {
  // onAddPeak/onRemovePeak are wired in Tasks 7–8; parameters are aliased to
  // underscore-prefixed locals here to silence TS6133 "declared but never used".
  void _onAddPeak; void _onRemovePeak;
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
    return () => { host.replaceChildren(); };
  }, [trace, peaks]);

  return <div ref={hostRef} className="w-full h-full" data-testid="trace-viewer" />;
}
