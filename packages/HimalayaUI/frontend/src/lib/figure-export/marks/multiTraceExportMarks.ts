// multiTraceExportMarks.ts — Plot marks for the MultiTracePlot export.
// Stacked traces in display_order, per-member labels at each band's
// y-position, peak ticks/labels gated on showPeakTicks / showPeakLabels.
import * as Plot from "@observablehq/plot";
import type { ComparisonMember } from "../../../api";
import type { Trace } from "../../../api";
import { computeYBands } from "../../comparison/yBands";
import {
  LIGHT_PALETTE,
  TRACE_STROKE_PX,
  PEAK_TICK_STROKE_PX,
} from "../presets";

export interface MultiTraceMarksArgs {
  /** Already filtered (null-exposure_id removed) and sorted by display_order. */
  members: ComparisonMember[];
  traces: Map<number, Trace>;
  /** Pre-resolved member labels (post-#73 contract). */
  displayLabelByMemberId: Map<number, string>;
  /** Pre-resolved colors keyed by member.id (computed on the UNFILTERED list
   *  per spec §MultiTracePlot export — distinct mode keys off the original
   *  index). */
  colorByMember: Map<number, string>;
  showPeakTicks: boolean;
  showPeakLabels: boolean;
  panelHeight: number;
}

export function buildMultiTraceExportMarks(args: MultiTraceMarksArgs): Plot.Markish[] {
  const {
    members, traces, displayLabelByMemberId, colorByMember,
    showPeakTicks, showPeakLabels, panelHeight,
  } = args;

  const ratios = members.map((m) => m.band_height || 1);
  const yBands = computeYBands(ratios, panelHeight);

  const marks: Plot.Markish[] = [];

  for (let i = 0; i < members.length; i++) {
    const member = members[i]!;
    if (member.exposure_id === null) continue; // defensive (caller pre-filtered)
    const trace = traces.get(member.exposure_id);
    if (!trace) continue;
    const band = yBands[i]!;
    const [bandTop, bandBottom] = band;
    const bandH = Math.max(1, bandBottom - bandTop);
    const color = colorByMember.get(member.id) ?? LIGHT_PALETTE.trace;

    // Trace line — y-mapped into this band. We compute log(I) ourselves to
    // place the line within the band (Plot's y-domain is per-figure, not
    // per-band). Each band has its own (per-member) y-range.
    const positiveI = trace.I.filter((v) => Number.isFinite(v) && v > 0);
    if (positiveI.length === 0) continue;
    // reduce() rather than Math.min(...positiveI) — V8's spread-arg cap
    // (~65k–500k depending on build) crashes on long SAXS traces seen via
    // /compare/all with high-resolution detectors. reduce() is unbounded.
    let minI = Infinity;
    let maxI = -Infinity;
    for (const v of positiveI) {
      if (v < minI) minI = v;
      if (v > maxI) maxI = v;
    }
    const logMin = Math.log10(minI);
    const logMax = Math.log10(maxI);
    const logRange = Math.max(1e-6, logMax - logMin);
    const points = trace.q.map((q, idx) => {
      const I = trace.I[idx] ?? 0;
      const li = I > 0 ? Math.log10(I) : logMin;
      // Map li ∈ [logMin, logMax] → y ∈ [bandBottom, bandTop] (inverted: low at bottom).
      const y = bandBottom - ((li - logMin) / logRange) * bandH;
      return { q, y };
    });
    marks.push(
      Plot.line(points, {
        x: "q",
        y: "y",
        stroke: color,
        strokeWidth: TRACE_STROKE_PX,
      }),
    );

    // Per-member label at the band's vertical midpoint, anchored to the left.
    const label = displayLabelByMemberId.get(member.id) ?? "";
    if (label) {
      marks.push(
        Plot.text(
          [{ q: trace.q[0] ?? 0, y: bandTop + bandH * 0.15, label }],
          {
            x: "q",
            y: "y",
            text: "label",
            textAnchor: "start",
            fill: LIGHT_PALETTE.text,
            fontSize: 11,
            dx: 4,
          },
        ),
      );
    }

    // Peaks — honour show toggles.
    const peaks = member.snapshot?.effective_peaks ?? [];
    if (showPeakTicks && peaks.length > 0) {
      const tickPoints = peaks.map((p) => ({
        q: p.q,
        y: bandTop + bandH * 0.05,
      }));
      marks.push(
        Plot.ruleX(tickPoints, {
          x: "q",
          stroke: color,
          strokeWidth: PEAK_TICK_STROKE_PX,
          // Restrict the tick height by mapping y to a small range above the band top.
        }),
      );
      if (showPeakLabels) {
        const labelRows = peaks.map((p, idx) => ({
          q: p.q,
          y: bandTop + bandH * 0.02,
          label: String(idx + 1),
        }));
        marks.push(
          Plot.text(labelRows, {
            x: "q",
            y: "y",
            text: "label",
            fill: color,
            fontSize: 9,
            dy: -2,
          }),
        );
      }
    }
  }

  return marks;
}
