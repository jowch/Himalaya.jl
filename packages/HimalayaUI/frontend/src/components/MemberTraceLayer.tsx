/**
 * MemberTraceLayer — per-member mark factory for the Compare page multi-trace
 * plot (Plan §Phase 6, Task 6.1; spec §Plot rendering).
 *
 * **Shape:** This is NOT a React component that renders its own DOM. It's a
 * pure mark factory: given a comparison member, its live `(q, I)` trace, and
 * the y-band allocated to it, it produces an array of Observable Plot marks
 * (line + peak dots + optional labels) that the parent `MultiTracePlot`
 * composes into a single `<Plot>`. One shared plot lets brush + zoom apply
 * uniformly across all traces — see spec §Plot rendering.
 *
 * Peak rendering rules (spec §Peak rendering and hover-driven coloring):
 * - Triangles render in BLACK by default — calm waterfall aesthetic.
 * - When `highlightedIndexId` matches `snapshot.confirmed_index.id`, that
 *   member's index peaks recolor to the phase color from `phases.ts`;
 *   non-index peaks stay black.
 * - `peak_display.hidden` peaks are omitted entirely.
 * - `peak_display.labeled` peaks get a text label above the triangle.
 *   v6.1 places the text directly above; Phase 8 adds leader-line dodge.
 *
 * Peaks come from `member.snapshot.effective_peaks` (the curation set
 * snapshotted at submit time) — NOT a live peaks query. See spec §Derived
 * analysis state and staleness for why we render against the snapshot.
 *
 * Y-positioning uses `applyNormalization` for the line trace; peak dots are
 * positioned at the same y as the underlying line at the peak's q (via
 * intensity lookup against the snapshot peak — manual peaks have null
 * intensity and fall back to the normalized line y at that q).
 *
 * **React wrapper not used yet.** v6.1 exposes only `buildMemberMarks` as the
 * tested API. `MultiTracePlot` calls it directly. A `<MemberTraceLayer>` JSX
 * component is exported as a thin convenience for callers that want a
 * declarative wrapper, but it just wraps `buildMemberMarks`.
 */
import * as Plot from "@observablehq/plot";
import type { Trace, ComparisonMember, MemberSnapshotPeak } from "../api";
import { phaseColor } from "../phases";
import {
  applyNormalization,
  computeReference,
  type Normalization,
  type QWindow,
} from "../lib/comparison/normalization";

const PEAK_COLOR_DEFAULT = "black";

/** Pixel offset above the trace at which peak triangles render. */
const PEAK_OFFSET_PX = 5;

/** Pixel offset above the peak triangle at which labels render (v6.1). */
const LABEL_OFFSET_PX = 12;

export interface MemberMarksProps {
  member: ComparisonMember;
  /**
   * Live `(q, I)` trace from `useTrace(exposureId)`. May be undefined while
   * the trace is loading; in that case no marks are emitted.
   */
  trace: Trace | undefined;
  /** Numeric y-band envelope `[topPx, bottomPx]` for this member. */
  yBand: [number, number];
  /** Per-member display state (from `member.peak_display` or draft). */
  peakDisplay?: { hidden: number[]; labeled: number[] } | undefined;
  /**
   * When set and matches `member.snapshot.confirmed_index.id`, peaks
   * belonging to that index recolor to the phase color.
   */
  highlightedIndexId?: number | undefined;
}

interface PeakRow {
  q: number;
  y: number;
  peakId: number;
  color: string;
}

/**
 * Build the array of Observable Plot marks for one member. Pure function —
 * safe to call inside a render. `Plot.line / Plot.dot / Plot.text` are
 * imported lazily through dynamic dispatch to keep the unit-test mock
 * surface minimal.
 */
export function buildMemberMarks(props: MemberMarksProps): unknown[] {
  const marks: unknown[] = [];
  const { member, trace, yBand, peakDisplay, highlightedIndexId } = props;

  if (!trace || trace.q.length === 0) return marks;

  const snapshot = member.snapshot;
  const peaks: MemberSnapshotPeak[] = snapshot ? snapshot.effective_peaks : [];
  const qWindow: QWindow =
    member.q_window_min !== null && member.q_window_max !== null
      ? [member.q_window_min, member.q_window_max]
      : null;
  const normalization = (member.normalization as Normalization) ?? "qwindow";

  const reference = computeReference(
    { q: trace.q, I: trace.I },
    peaks.map((p) => ({ q: p.q, intensity: p.intensity })),
    qWindow,
    normalization,
  );
  const linePoints = applyNormalization(
    { q: trace.q, I: trace.I },
    reference,
    yBand,
  );

  marks.push(
    Plot.line(linePoints, {
      x: "q",
      y: "y",
      stroke: member.color_override ?? "var(--color-fg)",
      strokeWidth: 1,
    }),
  );

  // Peaks. `effective_peaks` is rendered against the snapshot intensities;
  // the y-coordinate is derived from the same normalization mapping as the
  // line, so peak dots sit on (or just above) the line at their q.
  if (snapshot && peaks.length > 0) {
    const hidden = new Set<number>(peakDisplay?.hidden ?? []);
    const labeled = new Set<number>(peakDisplay?.labeled ?? []);

    // Highlighted = those belonging to the confirmed_index when its id matches.
    const highlightedPeakIds = new Set<number>(
      highlightedIndexId !== undefined
      && snapshot.confirmed_index !== null
      && snapshot.confirmed_index.id === highlightedIndexId
        ? snapshot.confirmed_index.peak_ids
        : [],
    );
    const highlightColor = snapshot.confirmed_index
      ? phaseColor(snapshot.confirmed_index.phase)
      : PEAK_COLOR_DEFAULT;

    const visiblePeaks: PeakRow[] = [];
    for (const p of peaks) {
      if (hidden.has(p.id)) continue;
      const lineY = interpolateLineY(linePoints, p.q);
      const y = Math.max(yBand[0], lineY - PEAK_OFFSET_PX);
      const color = highlightedPeakIds.has(p.id)
        ? highlightColor
        : PEAK_COLOR_DEFAULT;
      visiblePeaks.push({ q: p.q, y, peakId: p.id, color });
    }

    if (visiblePeaks.length > 0) {
      marks.push(
        Plot.dot(visiblePeaks, {
          x: "q",
          y: "y",
          symbol: "triangle",
          // Per-row color via channel — Observable Plot reads the `color`
          // field as the fill via the `fill` accessor below.
          fill: (d: unknown) => (d as PeakRow).color,
          stroke: "var(--color-bg)",
          strokeWidth: 0.75,
          r: 4,
        }),
      );
    }

    // Labels: v6.1 renders directly above the triangle. Phase 8 adds the
    // leader-line dodge layout.
    const labels = visiblePeaks.filter((p) => labeled.has(p.peakId));
    if (labels.length > 0) {
      marks.push(
        Plot.text(
          labels.map((p) => ({
            ...p,
            // Display label = q rounded to 3 sig digits — matches the q
            // label format used elsewhere; Phase 8 may swap for lattice-d
            // or Miller index later.
            label: p.q.toPrecision(3),
            y: Math.max(yBand[0], p.y - LABEL_OFFSET_PX),
          })),
          {
            x: "q",
            y: "y",
            text: "label",
            fill: "var(--color-fg)",
            fontSize: 10,
            textAnchor: "middle",
          },
        ),
      );
    }
  }

  return marks;
}

/**
 * Linear interpolation of the normalized line's y value at the given q.
 * Falls back to the nearest endpoint when the query is outside the trace.
 */
function interpolateLineY(line: Array<{ q: number; y: number }>, qTarget: number): number {
  if (line.length === 0) return 0;
  if (qTarget <= line[0]!.q) return line[0]!.y;
  if (qTarget >= line[line.length - 1]!.q) return line[line.length - 1]!.y;
  // Linear scan — N is the trace sample count; for a typical SAXS trace this
  // is ~1000 and the call frequency is per-render-of-one-member; binary
  // search would be marginal optimization.
  for (let i = 1; i < line.length; i++) {
    const a = line[i - 1]!;
    const b = line[i]!;
    if (qTarget >= a.q && qTarget <= b.q) {
      const t = (qTarget - a.q) / (b.q - a.q);
      return a.y + t * (b.y - a.y);
    }
  }
  return line[line.length - 1]!.y;
}

/**
 * Thin React component wrapper. `MultiTracePlot` doesn't actually mount this
 * — it calls `buildMemberMarks` directly because Observable Plot composes
 * marks at the plot factory level, not via JSX children. Exported so test
 * harnesses or future declarative renderers have a handle.
 */
export function MemberTraceLayer(props: MemberMarksProps): null {
  // A mark factory is not a renderable component on its own. Returning null
  // keeps it usable as a JSX placeholder if a parent wants to enumerate
  // children for some reason — but the actual marks must be produced via
  // `buildMemberMarks` and passed to `Plot.plot({ marks: ... })`.
  void props;
  return null;
}
