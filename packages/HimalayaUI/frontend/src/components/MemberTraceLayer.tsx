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
import type { Trace, SeriesMember, MemberSnapshotPeak } from "../api";
import { phaseColor } from "../phases";
import { peakMark } from "./ui/peakMark";
import {
  applyNormalization,
  computeReference,
  type Normalization,
  type QWindow,
} from "../lib/comparison/normalization";
import { layoutPeakLabels } from "../lib/plot/labelDodge";
import {
  colorFor,
  COMPARE_PALETTE,
  type GroupingMode,
} from "../lib/comparison/coloring";

const PEAK_COLOR_DEFAULT = "black";

/** Pixel offset above the trace at which peak triangles render. */
const PEAK_OFFSET_PX = 5;

/** Pixel offset above the peak triangle at which labels render (v6.1). */
const LABEL_OFFSET_PX = 12;

/**
 * Default label width estimate (px). Tuned for the q→3-sig-fig label
 * format ("0.123") at fontSize 10. Used as the minimum gap between two
 * adjacent labels by `layoutPeakLabels`.
 */
const DEFAULT_LABEL_WIDTH_PX = 32;

export interface MemberMarksProps {
  member: SeriesMember;
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
  /**
   * x-scale apply/invert pair for label dodge layout (Phase 8.2). When
   * provided, labels for crowded peaks are spread horizontally and a leader
   * line is emitted from the dodged label back to the triangle. Omitted on
   * the first render pass (before the plot exists); the dodge re-runs on
   * the second pass.
   */
  xScale?: { toPx: (q: number) => number; fromPx: (px: number) => number } | undefined;
  /**
   * Global annotation toggles (Phase 9.3). When `showPeakTicks` is false,
   * NO peak triangles render for this member regardless of `peak_display`
   * (`peak_display.hidden` per-peak hides still take effect when the
   * global flag is on). When `showPeakLabels` is false, no labels render
   * even for peaks marked `peak_display.labeled`. Both default to `true`
   * when omitted (review-mode default; edit mode passes `true` for both).
   */
  showPeakTicks?: boolean;
  showPeakLabels?: boolean;
  /**
   * Grouping mode for the per-member line stroke (Phase 9 gap-fix; spec
   * §Trace coloring). Wires the line color through the shared `colorFor`
   * library so toggling between bySample / byPhase / distinct visibly
   * recolors traces. `member.color_override` always wins (resolution order
   * is enforced inside `colorFor`).
   *
   * When `groupingMode` is omitted, the layer falls back to the legacy
   * "color_override → var(--color-ink)" behaviour. This keeps the few
   * non-Compare callers (and minimal test fixtures) working without
   * forcing them to wire the full `allMembers` + `sampleIdFor` context.
   */
  groupingMode?: GroupingMode;
  /**
   * All members in the same comparison (in display order). Required for
   * `bySample` and `distinct` palette indexing in `colorFor`. Ignored
   * when `groupingMode` is undefined.
   */
  allMembers?: ReadonlyArray<SeriesMember>;
  /**
   * Sample-id resolver supplied by the caller (typically wired against the
   * TanStack exposure cache at the page level). Returns `null` when the
   * exposure is unknown (cache miss / orphan). Ignored when `groupingMode`
   * is undefined.
   */
  sampleIdFor?: (m: SeriesMember) => number | null;
  /**
   * Fraction of each total band occupied by the working band (R8 offset
   * slider, #231). Forwarded to `applyNormalization`; when omitted the
   * library default (DEFAULT_WORKING_BAND_FRACTION = 0.7) applies. A larger
   * value makes traces taller within their band — the "tighter waterfall"
   * the builder's offset slider composes.
   */
  workingBandFraction?: number;
}

/**
 * One visible peak after hidden-filtering and y-positioning. Exported so the
 * parent `MultiTracePlot` can re-derive identical rows for click hit-testing
 * (Phase 8.1) and tooltip lookup (Phase 8.3) without re-implementing the
 * normalization math.
 */
export interface PeakRow {
  q: number;
  y: number;
  peakId: number;
  color: string;
  /** Provenance silhouette: manual → diamond, auto → triangle (peakMark). */
  source: "auto" | "manual";
}

/**
 * Compute the per-member visible-peak rows + the normalized line points used
 * by `buildMemberMarks`. Extracted so the parent `MultiTracePlot` can re-use
 * the same data for click hit-testing (Phase 8.1) and tooltip lookup
 * (Phase 8.3) without re-implementing the normalization math.
 *
 * Returns `null` for `peaks` (and an empty `linePoints`) when the trace is
 * missing or empty — matches `buildMemberMarks`'s "emit no marks" path.
 */
export function buildMemberPeakRows(props: MemberMarksProps): {
  peaks: PeakRow[];
  linePoints: Array<{ q: number; y: number }>;
} {
  const { member, trace, yBand, peakDisplay, highlightedIndexId, workingBandFraction } = props;

  if (!trace || trace.q.length === 0) return { peaks: [], linePoints: [] };

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
    workingBandFraction,
  );

  if (!snapshot || peaks.length === 0) return { peaks: [], linePoints };

  const hidden = new Set<number>(peakDisplay?.hidden ?? []);

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
    visiblePeaks.push({ q: p.q, y, peakId: p.id, color, source: p.source });
  }

  return { peaks: visiblePeaks, linePoints };
}

/**
 * Build the array of Observable Plot marks for one member. Pure function —
 * safe to call inside a render. `Plot.line / Plot.dot / Plot.text` are
 * imported lazily through dynamic dispatch to keep the unit-test mock
 * surface minimal.
 */
export function buildMemberMarks(props: MemberMarksProps): unknown[] {
  const marks: unknown[] = [];
  const { member, trace, yBand, peakDisplay } = props;
  const showPeakTicks  = props.showPeakTicks  ?? true;
  const showPeakLabels = props.showPeakLabels ?? true;

  if (!trace || trace.q.length === 0) return marks;

  const { peaks: visiblePeaks, linePoints } = buildMemberPeakRows(props);

  // Phase 9 gap-fix: route the line stroke through the shared `colorFor`
  // library so toggling the grouping mode visibly recolors traces. The
  // resolution order (`color_override` → grouping default → orphan
  // fallback) lives inside `colorFor` — keep the call site dumb.
  //
  // Fallback: when no grouping context is supplied, retain the legacy
  // behaviour so out-of-tree callers (and minimal test fixtures) don't
  // break. The Compare page always supplies the full triplet.
  const stroke = props.groupingMode !== undefined
    && props.allMembers !== undefined
    && props.sampleIdFor !== undefined
    ? colorFor(member, props.groupingMode, COMPARE_PALETTE, {
        allMembers: props.allMembers,
        sampleIdFor: props.sampleIdFor,
      })
    : (member.color_override ?? "var(--color-ink)");

  marks.push(
    Plot.line(linePoints, {
      x: "q",
      y: "y",
      stroke,
      strokeWidth: 1,
    }),
  );

  if (showPeakTicks && visiblePeaks.length > 0) {
    // No per-peak `data-testid="peak-tick"`: Plot.dot doesn't pass arbitrary
    // data-* attributes through to the rendered SVG. Per-peak interaction is
    // unit-tested via peakCycle.test.ts + this file's mark-arg capture tests;
    // the band-level `member-trace` overlay in MultiTracePlot covers E2E.
    // Converged onto the shared peakMark builder (Plan C plot spine): manual
    // peaks read as diamonds, auto as triangles, all in the per-row resolved
    // colour. `peakMark` reads `q`/`y`/`color`/`source` off each row.
    marks.push(peakMark(visiblePeaks));

    // Labels: when `xScale` is provided, run `layoutPeakLabels` to spread
    // crowded labels horizontally + emit leader lines from the dodged
    // text back to the triangle anchor. Without xScale (first render pass
    // before Plot creates the scale), fall back to direct-above placement;
    // the parent re-renders with the scale on the next pass.
    //
    // Phase 9.3: also gated on the global `showPeakLabels` flag — toggling
    // it off hides labels everywhere even when individual peaks are
    // `peak_display.labeled`.
    const labeled = new Set<number>(peakDisplay?.labeled ?? []);
    const labels = showPeakLabels
      ? visiblePeaks.filter((p) => labeled.has(p.peakId))
      : [];
    if (labels.length > 0) {
      if (props.xScale) {
        const dodged = layoutPeakLabels(
          labels.map((p) => ({ q: p.q, y: p.y, peakId: p.peakId })),
          {
            toPx: props.xScale.toPx,
            fromPx: props.xScale.fromPx,
            labelWidthPx: DEFAULT_LABEL_WIDTH_PX,
          },
        );
        marks.push(
          Plot.text(
            dodged.map((d) => ({
              q: d.qLabel,
              y: Math.max(yBand[0], d.yLabel),
              label: d.label,
              peakId: d.peakId,
            })),
            {
              x: "q",
              y: "y",
              text: "label",
              fill: "var(--color-ink)",
              fontSize: 10,
              textAnchor: "middle",
            },
          ),
        );
        // Leader lines: only emit for labels the dodge had to push sideways
        // (qLabel ≠ qPeak). Sparse labels sit directly above and need no
        // string. The link mark connects (qPeak, yPeak) → (qLabel, yLabel).
        const links = dodged.filter((d) => d.qLabel !== d.qPeak);
        if (links.length > 0) {
          marks.push(
            Plot.link(
              links.map((d) => ({
                qPeak: d.qPeak,
                yPeak: d.yPeak,
                qLabel: d.qLabel,
                yLabel: Math.max(yBand[0], d.yLabel),
              })),
              {
                x1: "qPeak", y1: "yPeak",
                x2: "qLabel", y2: "yLabel",
                stroke: "var(--color-ink-soft)",
                strokeWidth: 0.5,
                strokeOpacity: 0.7,
              },
            ),
          );
        }
      } else {
        marks.push(
          Plot.text(
            labels.map((p) => ({
              ...p,
              label: p.q.toPrecision(3),
              y: Math.max(yBand[0], p.y - LABEL_OFFSET_PX),
            })),
            {
              x: "q",
              y: "y",
              text: "label",
              fill: "var(--color-ink)",
              fontSize: 10,
              textAnchor: "middle",
            },
          ),
        );
      }
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
