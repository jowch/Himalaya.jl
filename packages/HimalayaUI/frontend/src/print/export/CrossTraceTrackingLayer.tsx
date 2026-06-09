/**
 * CrossTraceTrackingLayer — the cross-trace peak-tracking / reflection-linking
 * annotation layer in the shared `MultiTracePlot` render core (#208).
 *
 * **Shape:** A *cross*-trace annotation factory, fundamentally different from
 * `MemberTraceLayer` (which is per-member). `MemberTraceLayer` cannot host
 * it: the annotation visually links the SAME reflection across stacked
 * members, so the marks need the whole-stack view.
 *
 * For each phase that appears in two or more members' confirmed indices,
 * connect the k-th Miller-order peak across the carrying members. The
 * member's `confirmed_index.peak_ids` lists peaks in Miller-order order (the
 * authoritative ordering supplied by the Julia indexer), and the
 * `effective_peaks` snapshot gives the q value for each peak id — so the
 * k-th reflection's q is read directly from snapshot data without
 * duplicating the phase-ratio table from the Julia core. A polyline is
 * emitted per (phase, k) when at least two carrying members exist.
 *
 * Visual intent matches the mockup's `trackOn` branch
 * (`series-builder.html:730-755`): a thin coloured polyline per Miller
 * order, dimmed enough that it reads as a connector, not a primary mark.
 *
 * Y position: each connected vertex sits at the *baseline* of that member's
 * y-band (bottom edge in pixel-space — large y = bottom in the plot's
 * reversed y-axis), where peak ticks plant their root in
 * `MemberTraceLayer`'s waterfall layout. Heatmap rows use the same band
 * envelope, so the baseline-y here is the row's lower edge.
 */
import * as Plot from "@observablehq/plot";
import type { SeriesMember } from "../../api";
import { phaseColor } from "../../phases";

/**
 * Number of Miller orders to render per phase. Cap at 4 to match the
 * mockup; deeper orders typically don't add cross-trace insight and
 * crowd the layer.
 */
const MAX_MILLER_ORDERS = 4;

/**
 * Stroke width of the tracking polyline. Tuned to read clearly without
 * overpowering the underlying trace.
 */
const TRACK_STROKE_WIDTH = 1.0;

/**
 * Opacity ramp by Miller order — the strongest (lowest-q) reflection is the
 * most prominent connector. Mockup uses a 0.55..0.9 range.
 */
function opacityFor(order: number): number {
  return 0.55 + 0.4 * Math.pow(0.72, order);
}

export interface CrossTraceTrackingProps {
  /** Members in render order (top → bottom) — caller already sorts. */
  members: SeriesMember[];
  /**
   * Pixel y-bands per member (the parent plot's `computeYBands` output).
   * Indexed by member position (aligned with `members`). The tracking
   * vertex y is the *baseline* of each band — bottom edge in pixel-space.
   */
  yBands: Array<[number, number]>;
}

/**
 * Per-vertex row used by `Plot.line(..., { z: "key" })` so each polyline is
 * its own connected series. `key` groups the polyline; `q` + `y` are the
 * rendered coordinates; `phase` + `order` carry through for downstream
 * inspection / hover.
 */
interface TrackVertex {
  key: string;
  q: number;
  y: number;
  phase: string;
  order: number;
  stroke: string;
  opacity: number;
}

/**
 * Build the cross-trace tracking marks. Returns at most one `Plot.line`
 * mark grouped by `key` (phase + Miller order). Empty array when no phase
 * is shared by two-or-more indexed members — the layer renders only when
 * there's something to connect.
 *
 * For each carrying member: read `confirmed_index.peak_ids[k]`, look it up
 * in `effective_peaks` to get its q. Two-or-more carrying members per
 * (phase, k) form a polyline.
 */
export function buildCrossTraceTrackingMarks(
  props: CrossTraceTrackingProps,
): unknown[] {
  const { members, yBands } = props;

  // Phase → list of carrying-member rows. A row is the member's position +
  // its peak_ids array (peak_ids is already in Miller-order order from the
  // Julia indexer).
  const phaseRows = new Map<string, Array<{
    pos: number;
    peakIds: number[];
    peaks: ReadonlyArray<{ id: number; q: number }>;
  }>>();
  for (let i = 0; i < members.length; i++) {
    const m = members[i]!;
    const idx = m.snapshot?.confirmed_index;
    if (!idx) continue;
    const phase = idx.phase;
    const peaks = m.snapshot?.effective_peaks ?? [];
    const peakIds = idx.peak_ids ?? [];
    if (peakIds.length === 0 || peaks.length === 0) continue;
    if (!phaseRows.has(phase)) phaseRows.set(phase, []);
    phaseRows.get(phase)!.push({ pos: i, peakIds, peaks });
  }

  const vertices: TrackVertex[] = [];
  for (const [phase, rows] of phaseRows) {
    if (rows.length < 2) continue;
    // Phase colour, NOT terracotta. The migration polyline is a *phase
    // relation* (the same reflection tracked across members of one phase),
    // not an interaction mark — so it carries the phase hue, like every other
    // phase-coloured mark on the surface. DESIGN.md §2 Phase-Carries-the-
    // Surface Rule reserves terracotta (the grease pencil) for interaction
    // (reject ✕, primary action, q-link cross-highlight). Do not "fix" this to
    // the accent (cf. R3-F01: terracotta over-use is the recurring miss).
    const stroke = phaseColor(phase);
    // The longest peak_ids array bounds how many Miller orders we can connect
    // — orders k where fewer than 2 members carry a k-th peak are dropped.
    const maxK = Math.min(
      MAX_MILLER_ORDERS,
      rows.reduce((acc, r) => Math.max(acc, r.peakIds.length), 0),
    );
    for (let k = 0; k < maxK; k++) {
      const opacity = opacityFor(k);
      const key = `${phase}-${k}`;
      const rowVertices: TrackVertex[] = [];
      for (const r of rows) {
        if (k >= r.peakIds.length) continue;
        const targetId = r.peakIds[k]!;
        const peak = r.peaks.find((p) => p.id === targetId);
        if (!peak) continue;
        const band = yBands[r.pos];
        if (!band) continue;
        const yBaseline = band[1]; // bottom-of-band in pixel space
        rowVertices.push({
          key,
          q: peak.q,
          y: yBaseline,
          phase,
          order: k,
          stroke,
          opacity,
        });
      }
      // A polyline needs ≥ 2 points to read as a line.
      if (rowVertices.length >= 2) {
        for (const v of rowVertices) vertices.push(v);
      }
    }
  }

  if (vertices.length === 0) return [];

  return [
    Plot.line(vertices, {
      x: "q",
      y: "y",
      z: "key",
      stroke: (d: unknown) => (d as TrackVertex).stroke,
      strokeOpacity: (d: unknown) => (d as TrackVertex).opacity,
      strokeWidth: TRACK_STROKE_WIDTH,
    }),
  ];
}
