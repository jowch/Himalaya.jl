import type { SeriesMember, Trace, AssignmentState } from "../../api";
import { traceIntensityAt } from "../../lib/plot/traceIntensity";
import {
  dominantPhase,
  indexedAnchorPeaks,
  memberRowLabel,
  resolveState,
} from "../../lib/series/memberRead";

// Single-shared-predicate (BU-EXPORTDIVERGE): the plate's identity reads live
// in lib/series/memberRead so the export pipeline resolves phase / label /
// state / peak set through the exact same chains. Re-exported here for the
// existing waterfallModel importers (builderAdapters, folioAdapters, …).
export { dominantPhase, resolveState } from "../../lib/series/memberRead";

const EMPTY_TRACE: Trace = { q: [], I: [], sigma: [] };

/** One indexed-peak bead position on a member's trace. */
export interface WaterfallAnchor {
  /** Effective-peak id (from peak_ids) — maps 1:1 to PlotPeak.id. */
  id: number;
  q: number;
  /** Measured peak intensity (from effective_peaks) — places the bead on the curve. */
  intensity: number | null;
  /** Member dominant phase (drives bead colour); null = unindexed. */
  phase: string | null;
}

/** One member row, fully resolved for rendering (low→high variable order). */
export interface WaterfallRow {
  key: string;
  label: string;
  trace: Trace;
  phase: string | null;
  state: AssignmentState;
  anchors: WaterfallAnchor[];
  bandHeight: number;
  yOffset: number;
}

/** Map API members + a trace lookup → rows in input (display) order. */
export function toWaterfallRows(
  members: SeriesMember[],
  tracesById: Record<number, Trace>,
): WaterfallRow[] {
  return members.map((member) => {
    const state = resolveState(member);
    const phase = dominantPhase(member);
    const trace = (member.exposure_id != null && tracesById[member.exposure_id]) || EMPTY_TRACE;

    // indexedAnchorPeaks is the shared annotated-peak predicate (also drives
    // the export's peak marks) — indexed anchors only, ascending q.
    // Anchor every bead to the curve. A manually-added peak has no measured
    // intensity (`null`) — without this it would drop to the row baseline
    // (PlotPeaks `intensity ?? baselineI`) and render flat on the floor instead
    // of on its trace. Mirrors the Focus `toTraceModel` adapter via the shared
    // `traceIntensityAt`.
    const anchors: WaterfallAnchor[] = indexedAnchorPeaks(member).map((p) => ({
      id: p.id,
      q: p.q,
      intensity: p.intensity != null ? p.intensity : traceIntensityAt(p.q, trace),
      phase,
    }));

    return {
      key: String(member.id),
      label: memberRowLabel(member),
      trace,
      phase,
      state,
      anchors,
      bandHeight: member.band_height,
      yOffset: member.y_offset,
    };
  });
}

/** Positive q-extent across all rows, padded ×0.97 / ×1.03. */
export function waterfallQDomain(rows: WaterfallRow[]): [number, number] {
  let lo = Infinity;
  let hi = 0;
  for (const r of rows) {
    for (const q of r.trace.q) {
      if (q > 0 && Number.isFinite(q)) {
        if (q < lo) lo = q;
        if (q > hi) hi = q;
      }
    }
  }
  if (!Number.isFinite(lo) || hi <= 0) return [0.01, 1]; // no positive data fallback
  return [lo * 0.97, hi * 1.03];
}
