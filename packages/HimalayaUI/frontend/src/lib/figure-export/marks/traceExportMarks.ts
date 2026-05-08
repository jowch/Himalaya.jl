// traceExportMarks.ts — Plot marks for the TraceViewer export.
// Duplicates on-screen overlay geometry (peaks + predicted-q ticks) — the
// on-screen overlay is load-bearing for hit-testing/hover and must not be
// touched. See spec §Approach for the accepted-tradeoff rationale.
import * as Plot from "@observablehq/plot";
import type { Trace, Peak, IndexEntry } from "../../../api";
import { phaseColor } from "../../../phases";
import {
  LIGHT_PALETTE,
  TRACE_STROKE_PX,
  PEAK_TICK_STROKE_PX,
  PREDICTED_Q_STROKE_PX,
} from "../presets";

/** Returns Plot marks. Caller threads them into ExportSpec.plot.marks. */
export function buildTraceExportMarks(args: {
  trace: Trace;
  peaks: Peak[];
  activeGroupIndices: IndexEntry[];
}): Plot.Mark[] {
  const { trace, peaks, activeGroupIndices } = args;

  const marks: Plot.Mark[] = [];

  // Trace line. Map (q,I) zip → array-of-objects for Plot.
  const points = trace.q.map((q, i) => ({ q, I: trace.I[i] ?? 0 }));
  marks.push(
    Plot.line(points, {
      x: "q",
      y: "I",
      stroke: LIGHT_PALETTE.trace,
      strokeWidth: TRACE_STROKE_PX,
    }),
  );

  // Peak triangles. Color by source/excluded.
  const peakRows = peaks.map((p) => ({
    q: p.q,
    I: p.intensity ?? 0,
    fill: p.excluded
      ? LIGHT_PALETTE.peakExcluded
      : p.source === "manual"
        ? LIGHT_PALETTE.peakManual
        : LIGHT_PALETTE.peakAuto,
  }));
  if (peakRows.length > 0) {
    marks.push(
      Plot.dot(peakRows, {
        x: "q",
        y: "I",
        fill: "fill",
        symbol: "triangle2",
        r: 4,
        stroke: "fill",
        strokeWidth: PEAK_TICK_STROKE_PX,
      }),
    );
  }

  // Predicted-q ticks per active-group index, phase-coloured.
  for (const idx of activeGroupIndices) {
    const color = phaseColor(idx.phase);
    const ticks = idx.predicted_q.map((q) => ({ q, phase: idx.phase }));
    marks.push(
      Plot.ruleX(ticks, {
        x: "q",
        stroke: color,
        strokeWidth: PREDICTED_Q_STROKE_PX,
        strokeOpacity: 0.85,
      }),
    );
  }

  return marks;
}
