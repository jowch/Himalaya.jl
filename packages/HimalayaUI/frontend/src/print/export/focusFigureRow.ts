// focusFigureRow.ts — turn the Focus page's carried state (trace + observed
// peaks + active assignment) into the SINGLE WaterfallRow the clean export
// figure renders. This is the Focus analog of toWaterfallRows: the series figure
// anchors a member's `indexedAnchorPeaks`; the Focus figure anchors the peaks
// the ACTIVE assignment claims (activeIndices[].peaks[].peak_id), so the
// exported figure shows the confirmed assignment, not every detected peak the
// working surface paints. Both feed the same buildCleanFigureSvg builder.
import type { Trace, Peak, IndexEntry } from "../../api";
import type { WaterfallRow } from "../waterfall/waterfallModel";
import { traceIntensityAt } from "../../lib/plot/traceIntensity";

export interface FocusFigureInput {
  trace: Trace;
  /** All observed peaks on the exposure (the byId lookup source). */
  peaks: Peak[];
  /** The active assignment's index entries (claim the figure's anchors). */
  activeIndices: IndexEntry[];
  /** Primary assigned phase (drives trace + anchor colour); null = unindexed. */
  phase: string | null;
  /** Right-gutter row label. */
  label: string;
}

/** Build the one waterfall row the clean Focus figure renders. */
export function buildFocusFigureRow(input: FocusFigureInput): WaterfallRow {
  const { trace, peaks, activeIndices, phase, label } = input;

  const claimed = new Set<number>();
  for (const ix of activeIndices) for (const p of ix.peaks) claimed.add(p.peak_id);
  const byId = new Map(peaks.map((p) => [p.id, p]));

  // Anchor each claimed peak to the curve. A manual/curation peak carries no
  // measured intensity — read it off the trace at its q (mirrors toWaterfallRows
  // / toTraceModel) so the glyph rides the data, not the baseline.
  const anchors = [...claimed]
    .map((id) => byId.get(id))
    .filter((p): p is Peak => p != null)
    .map((p) => ({
      id: p.id,
      q: p.q,
      intensity: p.intensity != null ? p.intensity : traceIntensityAt(p.q, trace),
      phase,
    }))
    .sort((a, b) => a.q - b.q);

  return {
    key: "focus",
    label,
    trace,
    phase,
    state: phase != null ? "indexed" : "null",
    anchors,
    bandHeight: 1,
    yOffset: 0,
  };
}
