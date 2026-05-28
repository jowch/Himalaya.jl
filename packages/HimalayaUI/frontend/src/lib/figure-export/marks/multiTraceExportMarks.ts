// multiTraceExportMarks.ts — Plot marks for the MultiTracePlot export.
// Stacked traces in display_order, per-member labels at each band's
// y-position, peak ticks/labels gated on showPeakTicks / showPeakLabels.
// Heatmap representation + cross-trace tracking layer mirror the on-screen
// MultiTracePlot (#251 r1 / B1).
import * as Plot from "@observablehq/plot";
import type { SeriesMember } from "../../../api";
import type { Trace } from "../../../api";
import { computeYBands } from "../../comparison/yBands";
import {
  LIGHT_PALETTE,
  TRACE_STROKE_PX,
  PEAK_TICK_STROKE_PX,
} from "../presets";
import { buildMemberHeatmapMarks } from "../../../components/MemberHeatmapLayer";
import { buildCrossTraceTrackingMarks } from "../../../components/CrossTraceTrackingLayer";
import type { Representation } from "../../../components/RepresentationToggle";

export interface MultiTraceMarksArgs {
  /** Already filtered (null-exposure_id removed) and sorted by display_order. */
  members: SeriesMember[];
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
  /** Render mode (#251 r1 / B1). Defaults to `"waterfall"` — the legacy
   *  export shape. `"heatmap"` swaps the per-row mark vocabulary to binned
   *  intensity cells while keeping the same y-band envelope. */
  representation?: Representation;
  /** Emit the cross-trace tracking polylines (#251 r1 / B1). */
  showCrossTraceTracking?: boolean;
  /** Visible q-domain, threaded to the heatmap binner. `null` → derive from
   *  the underlying traces (rare; the adapter normally supplies xDomain). */
  xDomain?: [number, number] | null;
}

/**
 * Pure white plate override for heatmap cells in the export. The on-screen
 * heatmap mixes intensity into the warm `--plate` paper hue; on a white-bg
 * export the cells should mix into pure white so the un-tinted background
 * actually reads as background, not as a faint warm cast.
 */
const EXPORT_PLATE_WHITE = "oklch(1 0 0)";

export function buildMultiTraceExportMarks(args: MultiTraceMarksArgs): Plot.Markish[] {
  const {
    members, traces, displayLabelByMemberId, colorByMember,
    showPeakTicks, showPeakLabels, panelHeight,
    representation = "waterfall",
    showCrossTraceTracking = false,
    xDomain,
  } = args;

  const ratios = members.map((m) => m.band_height || 1);
  const yBands = computeYBands(ratios, panelHeight);

  const marks: Plot.Markish[] = [];

  // Derive the heatmap q-domain once. Prefer the caller-supplied xDomain
  // (matches the on-screen brush). Fall back to the min/max across all
  // traces — only used when the adapter is called without a domain.
  const heatmapDomain = ((): [number, number] => {
    if (xDomain) return xDomain;
    let lo = Infinity;
    let hi = -Infinity;
    for (const m of members) {
      if (m.exposure_id === null) continue;
      const t = traces.get(m.exposure_id);
      if (!t || t.q.length === 0) continue;
      const q0 = t.q[0]!;
      const qN = t.q[t.q.length - 1]!;
      if (q0 < lo) lo = q0;
      if (qN > hi) hi = qN;
    }
    if (!Number.isFinite(lo) || !Number.isFinite(hi) || lo >= hi) {
      return [0.01, 1] as [number, number];
    }
    return [lo, hi];
  })();

  for (let i = 0; i < members.length; i++) {
    const member = members[i]!;
    if (member.exposure_id === null) continue; // defensive (caller pre-filtered)
    const trace = traces.get(member.exposure_id);
    if (!trace) continue;
    const band = yBands[i]!;
    const [bandTop, bandBottom] = band;
    const bandH = Math.max(1, bandBottom - bandTop);
    const color = colorByMember.get(member.id) ?? LIGHT_PALETTE.trace;

    if (representation === "heatmap") {
      // Reuse the on-screen heatmap factory so binning + contrast curve stay
      // in lockstep. Pre-resolved fill bypasses the layer's COMPARE_PALETTE
      // lookup (the export adapter has already walked COMPARE_PALETTE_LIGHT
      // into colorByMember). Pure-white plate for the export's white bg.
      const heatmapMarks = buildMemberHeatmapMarks({
        member,
        trace,
        yBand: [bandTop, bandBottom],
        qDomain: heatmapDomain,
        fillBaseOverride: color,
        plateOverride: EXPORT_PLATE_WHITE,
      });
      for (const mk of heatmapMarks) marks.push(mk as Plot.Markish);

      // Heatmap rows still get the per-member label at the band midpoint —
      // matches the on-screen left-margin label (the in-band placement is
      // the only equivalent without an axis margin).
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
      // Heatmap representation folds peaks into the intensity field; skip
      // the peak ticks / labels block.
      continue;
    }

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

  // Cross-trace peak-tracking polylines (#251 r1 / B1). Pushed AFTER the
  // per-member marks so the connectors render on top of the waterfall lines
  // / heatmap cells — same z-order as the on-screen MultiTracePlot.
  if (showCrossTraceTracking) {
    const trackingMarks = buildCrossTraceTrackingMarks({
      members,
      yBands: yBands as Array<[number, number]>,
    });
    for (const mk of trackingMarks) marks.push(mk as Plot.Markish);
  }

  return marks;
}
