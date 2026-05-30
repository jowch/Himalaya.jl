// traceExportMarks.ts — Plot marks for the TraceViewer export.
// Duplicates on-screen overlay geometry (peaks + predicted-q ticks) — the
// on-screen overlay is load-bearing for hit-testing/hover and must not be
// touched. See spec §Approach for the accepted-tradeoff rationale.
import * as Plot from "@observablehq/plot";
import type { Trace, Peak, IndexEntry } from "../../../api";
import { phaseColor } from "../../../phases";
import { peakMark } from "../../../components/ui/peakMark";
import {
  LIGHT_PALETTE,
  TRACE_STROKE_PX,
  PREDICTED_Q_STROKE_PX,
} from "../presets";

/** Returns Plot marks. Caller threads them into ExportSpec.plot.marks.
 *
 *  Return type is `Plot.Markish[]` (not `Plot.Mark[]`) to match Plot's
 *  `PlotOptions.marks` field — `Markish` is the public-facing type Plot
 *  expects for caller-supplied marks. The runtime values produced by
 *  `Plot.line / Plot.dot / Plot.ruleX` satisfy both. */
export function buildTraceExportMarks(args: {
  trace: Trace;
  peaks: Peak[];
  activeGroupIndices: IndexEntry[];
}): Plot.Markish[] {
  const { trace, peaks, activeGroupIndices } = args;

  const marks: Plot.Markish[] = [];

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

  // Peak glyphs via the shared peakMark builder (Plan C plot spine). Manual
  // peaks now read as DIAMONDS, auto as triangles. Export colour stays
  // BY-SOURCE (LIGHT_PALETTE), NOT by-phase — the printable palette is
  // deliberately source-coded, so we thread the LIGHT_PALETTE literal into
  // peakMark's resolved `color` param (it carries literals fine). Excluded
  // peaks ghost via the builder's `excluded` flag.
  const peakRows = peaks.map((p) => ({
    q: p.q,
    y: p.intensity ?? 0,
    color: p.source === "manual" ? LIGHT_PALETTE.peakManual : LIGHT_PALETTE.peakAuto,
    source: p.source,
    ...(p.excluded ? { excluded: true } : {}),
  }));
  if (peakRows.length > 0) {
    marks.push(peakMark(peakRows, { y: "y" }));
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
