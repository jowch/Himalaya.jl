/**
 * MultiTracePlot — shared Observable Plot host for the Compare page
 * (Plan §Phase 6, Task 6.2; spec §Plot rendering).
 *
 * Single `<Plot>` instance with one shared q-scale and one shared zoom
 * domain so brush + double-click reset uniformly across all traces. Marks
 * are produced per member by `MemberTraceLayer.buildMemberMarks`, then
 * concatenated into the plot's `marks` array.
 *
 * Y-bands are derived from each member's `band_height` ratio per spec
 * formula:
 *
 *   panel_height = container.clientHeight (from a ResizeObserver)
 *   band_i.height = (band_height_i / Σ band_heights) × panel_height
 *
 * Members are stacked top-down in the order they're passed to this
 * component (the parent already sorts by `display_order`).
 *
 * Aspect ratio is hardcoded at `COMPARE_PLOT_ASPECT = 0.3` (3:10 W:H) per
 * spec; exported so a future surface can flip it. The aspect is enforced
 * by the parent layout via `aspect-ratio` CSS — this component does not
 * impose a width.
 *
 * Brush/zoom semantics mirror the existing `TraceViewer` (mouse-wheel
 * zoom in q, double-click to reset). Y is fixed by bands; this plot does
 * not zoom y. Peak click semantics (edit-mode cycle through shown / labeled
 * / hidden) land in Phase 8 — for v6.2 we host the plot only.
 */
import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import type { Trace, SeriesMember } from "../api";
import { buildMemberMarks, buildMemberPeakRows } from "./MemberTraceLayer";
import type { PeakRow } from "./MemberTraceLayer";
import { buildMemberHeatmapMarks } from "./MemberHeatmapLayer";
import { buildCrossTraceTrackingMarks } from "./CrossTraceTrackingLayer";
import { PlotSurface } from "./PlotSurface";
import type { PlotOverlayContext } from "./PlotSurface";
import { SeriesTrackingOverlay } from "./SeriesTrackingOverlay";
import { buildAnchorMap } from "../lib/series/anchors";
import type { GroupingMode } from "../lib/comparison/coloring";
import { computeYBands } from "../lib/comparison/yBands";
import { useActiveBand } from "./ActiveBandContext";

/**
 * Plot layout vocabulary (#208 — render-core finish). Today the render core
 * has only the waterfall layout; the heatmap is a parallel representation
 * built from `MemberHeatmapLayer`. Unlike `GroupingMode` (a coloring axis),
 * `Representation` is a layout mode — picking one swaps the per-member mark
 * vocabulary (line+peaks vs binned rects).
 */
export type Representation = "waterfall" | "heatmap";

/**
 * Pixel hit radius for peak click hit-testing in edit mode (Phase 8.1).
 * Triangle radius is `r = 4` so a 10 px radius covers the marker plus a
 * comfortable interaction buffer.
 */
const PEAK_HIT_PX = 10;

/**
 * Synthesize an x-scale `(toPx, fromPx)` matching the active plot's domain
 * and width. Used by the label-dodge layout (Phase 8.2) so we can dodge in
 * the same render pass as mark building, without needing two `Plot.plot()`
 * calls. Match's Plot's clamping behaviour at the domain edges.
 */
function makeXScale(
  type: "log" | "linear",
  xMin: number,
  xMax: number,
  innerLeftPx: number,
  innerWidthPx: number,
): { toPx: (q: number) => number; fromPx: (px: number) => number } {
  if (type === "log") {
    const lo = Math.log(Math.max(xMin, 1e-12));
    const hi = Math.log(Math.max(xMax, lo + 1e-12));
    const span = hi - lo;
    return {
      toPx: (q) => innerLeftPx + ((Math.log(Math.max(q, 1e-12)) - lo) / span) * innerWidthPx,
      fromPx: (px) => Math.exp(lo + ((px - innerLeftPx) / innerWidthPx) * span),
    };
  }
  const span = xMax - xMin;
  return {
    toPx: (q) => innerLeftPx + ((q - xMin) / span) * innerWidthPx,
    fromPx: (px) => xMin + ((px - innerLeftPx) / innerWidthPx) * span,
  };
}

/** Hardcoded plot aspect ratio (W / H) per spec §Plot rendering. */
export const COMPARE_PLOT_ASPECT = 0.3;

/**
 * Per-band aspect target (W / H). 1.0 = each member's band is square. Issue
 * #81: SAXS stack-plot convention is taller-than-wide bands so peak features
 * aren't visually compressed. The acceptance bar is ≤ 1.5:1; 1.0 leaves
 * headroom while remaining sensible at small member counts.
 */
const BAND_ASPECT_TARGET = 1.0;

/**
 * Lower bound on the actual plot width (px). Without this floor, a single-
 * member case would shrink the plot to band-height-sized which is too narrow
 * to read q-axis ticks. ~280 px keeps default tick spacing legible and lets
 * the per-band aspect target stay reachable at N≥4 in typical panel heights
 * (a higher floor consistently overrode the target for large member counts).
 */
export const MIN_PLOT_WIDTH = 280;

/**
 * Compute the plot's `maxWidth` cap in CSS pixels for the inner wrapper that
 * centers the plot inside the slot (issue #81). Pure helper so the math is
 * unit-testable without rendering the full component.
 *
 * Each band is `panelHeight / memberCount` tall; we target a per-band aspect
 * of `BAND_ASPECT_TARGET` (1.0 = square), so `plotW = bandH × target`. The
 * `MIN_PLOT_WIDTH` floor keeps small panels / single-member cases legible.
 *
 * Degenerate inputs (`memberCount === 0`, `panelHeight ≤ 0`) fall back to
 * the floor so initial mount + JSDOM tests don't produce NaN/negative.
 */
export function computeMaxPlotWidth(
  panelHeight: number,
  memberCount: number,
): number {
  if (memberCount <= 0 || panelHeight <= 0) return MIN_PLOT_WIDTH;
  const bandH = panelHeight / memberCount;
  const target = bandH * BAND_ASPECT_TARGET;
  return Math.max(MIN_PLOT_WIDTH, target);
}

const MARGIN_LEFT   = 50;
const MARGIN_RIGHT  = 14;
const MARGIN_TOP    = 8;
const MARGIN_BOTTOM = 32;

type Scale = { invert?: (v: number) => number; apply?: (v: number) => number } | undefined;

/**
 * Map the builder's trace-offset slider (0.4..1.4, mockup default 1.2) to the
 * working-band fraction consumed by `applyNormalization` (R8, #231). A larger
 * offset = taller traces filling more of each band = a tighter waterfall.
 * Clamped so the slider can never collapse a band (min 0.45) or overflow it
 * (max 0.95).
 */
export function offsetToBandFraction(offset: number): number {
  const t = Math.min(1, Math.max(0, (offset - 0.4) / (1.4 - 0.4)));
  return 0.45 + t * (0.95 - 0.45);
}

export interface MultiTracePlotProps {
  /** Members in render order (top → bottom). Caller sorts by `display_order`. */
  members: SeriesMember[];
  /** Live traces keyed by exposure_id. Members without a trace render no line. */
  traces: Map<number, Trace>;
  /** Visible q-range. null = full data range. */
  xDomain: [number, number] | null;
  onXDomain: (d: [number, number] | null) => void;
  /** Per-member display state overrides (edit-mode draft, optional). */
  peakDisplayByMemberId?: Map<number, { hidden: number[]; labeled: number[] }>;
  /**
   * When set, peaks belonging to that member's confirmed_index recolor to
   * the phase color. Driven by hover/click on the metadata gutter row
   * (Phase 7 wires this; v6.2 just plumbs the prop).
   */
  highlightedMemberId?: number | undefined;
  /** Q-axis units label (defaults to "Å⁻¹"). */
  qUnits?: string;
  /** X-axis scale. Defaults to "log". */
  xType?: "log" | "linear";
  /**
   * Global annotation toggles (Phase 9.3). When `showPeakTicks` is false,
   * NO peak triangles render across all members; when `showPeakLabels` is
   * false, no labels render either. Both default to `true`. The parent
   * page (review = ComparePage) reads them from Zustand and passes them
   * down; edit mode (ComparePageEdit) leaves them at the defaults so the
   * peak-click cycle stays usable.
   */
  showPeakTicks?: boolean;
  showPeakLabels?: boolean;
  /**
   * Edit-mode peak click handler (Phase 8.1). When set, the plot installs a
   * click listener that does pixel-distance hit-testing against the rendered
   * peak rows (per member band) and dispatches `onPeakClick(memberId,
   * peakId, altKey)` for the closest peak within `PEAK_HIT_PX`. Omitting
   * this prop disables click handling entirely (review-mode default).
   */
  onPeakClick?: (memberId: number, peakId: number, altKey: boolean) => void;
  /**
   * Trace coloring grouping mode (Phase 9 gap-fix; spec §Trace coloring).
   * Forwarded to `MemberTraceLayer.buildMemberMarks` so the per-member line
   * stroke routes through `colorFor`. Omit to fall back to the legacy
   * "color_override → var(--color-ink)" behaviour (e.g. tests that don't
   * exercise grouping).
   */
  groupingMode?: GroupingMode;
  /**
   * Sample-id resolver supplied by the parent page (typically wired against
   * the TanStack exposure cache). `colorFor` consults this in `bySample`
   * mode. Required when `groupingMode` is set; ignored otherwise.
   */
  sampleIdFor?: (m: SeriesMember) => number | null;
  /**
   * Working-band fraction for the waterfall offset slider (R8, #231). When
   * omitted the layer falls back to `DEFAULT_WORKING_BAND_FRACTION`. The page
   * derives this from its offset slider via `offsetToBandFraction`.
   */
  workingBandFraction?: number;
  /**
   * Plot layout vocabulary (#208 — render-core finish). Default
   * `"waterfall"` preserves the legacy line+peaks behaviour. `"heatmap"`
   * swaps the per-member mark factory to `MemberHeatmapLayer`, rendering each
   * member as a horizontal row of intensity-binned rectangles instead of a
   * stacked trace. Peak ticks, peak labels, click hit-testing, and the
   * hover tooltip are all waterfall-only: peaks are folded into the
   * intensity field in the heatmap and have no per-peak target.
   */
  representation?: Representation;
  /**
   * Cross-trace peak-tracking layer (#208). When `true`, draw a thin
   * coloured polyline per (phase, Miller-order) connecting the same
   * reflection across every member whose confirmed_index carries that
   * phase. The connector reads as a migration line as q drifts with
   * the ordering variable. Default `false`. Renders in both waterfall
   * and heatmap; peaks come from the member snapshot, the indexer
   * supplies the Miller-order ordering via `confirmed_index.peak_ids`.
   */
  showCrossTraceTracking?: boolean;
}

export function MultiTracePlot(props: MultiTracePlotProps): JSX.Element {
  const {
    members, traces, xDomain, onXDomain,
    peakDisplayByMemberId, highlightedMemberId,
    qUnits, xType = "log",
    onPeakClick,
    showPeakTicks = true, showPeakLabels = true,
    groupingMode, sampleIdFor,
    workingBandFraction,
    representation = "waterfall",
    showCrossTraceTracking = false,
  } = props;

  const hostRef       = useRef<HTMLDivElement>(null);
  const plotContainer = useRef<HTMLDivElement>(null);

  const [_resizeKey, setResizeKey] = useState(0);
  // Tracked panel height drives the per-band overlay positions + the band math.
  const [panelHeight, setPanelHeight] = useState(0);
  const [panelWidth, setPanelWidth] = useState(0);

  // Phase 8.3 — peak hover tooltip state (unchanged contract).
  const [tooltip, setTooltip] = useState<{
    q: number; peakId: number; xPx: number; yPx: number;
  } | null>(null);
  const showPeakIds = typeof window !== "undefined"
    && new URLSearchParams(window.location.search).has("showPeakIds");

  // E-3: the ephemeral hover/focus-tracked (phase,order) anchor key. Hovering
  // an anchor bead threads the terracotta migration connector + ghost rings
  // across the stack; leaving clears it. Local UI state (per the plan).
  const [trackedKey, setTrackedKey] = useState<string | null>(null);

  useEffect(() => {
    const el = plotContainer.current;
    if (!el) return;
    setPanelHeight(el.clientHeight);
    setPanelWidth(el.clientWidth);
    const obs = new ResizeObserver(() => {
      setResizeKey((k) => k + 1);
      if (plotContainer.current) {
        setPanelHeight(plotContainer.current.clientHeight);
        setPanelWidth(plotContainer.current.clientWidth);
      }
    });
    obs.observe(el);
    return () => obs.disconnect();
  }, []);

  // Derive q-domain from union of all member traces (full-range default).
  const qExtent = useCallback((): [number, number] | null => {
    let lo = Infinity;
    let hi = -Infinity;
    for (const m of members) {
      if (m.exposure_id === null) continue;
      const t = traces.get(m.exposure_id);
      if (!t || t.q.length === 0) continue;
      const a = t.q[0]!;
      const b = t.q[t.q.length - 1]!;
      if (a < lo) lo = a;
      if (b > hi) hi = b;
    }
    if (!Number.isFinite(lo) || !Number.isFinite(hi) || hi <= lo) return null;
    return [lo, hi];
  }, [members, traces]);

  // Build the marks + per-member hit-test index in one memoized pass so they
  // never drift from each other (same `buildMemberPeakRows` source). `_resizeKey`
  // is a dep so the synthesized x-scale (used for label dodge) tracks resize.
  const { allMarks, hoverPeakIndex } = useMemo(() => {
    void _resizeKey;
    const panelW = panelWidth || 400;
    const panelH = panelHeight || 300;
    const ratios = members.map((m) => m.band_height || 1);
    const yBands = computeYBands(ratios, panelH);

    const ext = qExtent();
    const xMin = xDomain ? xDomain[0] : ext?.[0] ?? 1;
    const xMax = xDomain ? xDomain[1] : ext?.[1] ?? 10;
    const innerW = Math.max(1, panelW - MARGIN_LEFT - MARGIN_RIGHT);
    const xScale = makeXScale(xType, xMin, xMax, MARGIN_LEFT, innerW);

    const marks: unknown[] = [];
    const index: Array<{ memberId: number; yBand: [number, number]; peaks: PeakRow[] }> = [];
    for (let i = 0; i < members.length; i++) {
      const m = members[i]!;
      const yBand = yBands[i] ?? [0, panelH];
      const trace = m.exposure_id !== null ? traces.get(m.exposure_id) : undefined;
      const peakDisplay = peakDisplayByMemberId?.get(m.id);
      const layerProps = {
        member: m,
        trace,
        yBand: yBand as [number, number],
        peakDisplay,
        highlightedIndexId:
          highlightedMemberId === m.id && m.snapshot?.confirmed_index
            ? m.snapshot.confirmed_index.id
            : undefined,
        xScale,
        showPeakTicks,
        showPeakLabels,
        ...(groupingMode !== undefined && sampleIdFor !== undefined
          ? { groupingMode, allMembers: members, sampleIdFor }
          : {}),
        ...(workingBandFraction !== undefined ? { workingBandFraction } : {}),
      };
      if (representation === "heatmap") {
        const heatmapMarks = buildMemberHeatmapMarks({
          member: m,
          trace,
          yBand: yBand as [number, number],
          qDomain: [xMin, xMax],
          ...(groupingMode !== undefined && sampleIdFor !== undefined
            ? { groupingMode, allMembers: members, sampleIdFor }
            : {}),
        });
        for (const mk of heatmapMarks) marks.push(mk);
        continue;
      }
      const memberMarks = buildMemberMarks(layerProps);
      for (const mk of memberMarks) marks.push(mk);
      const { peaks } = buildMemberPeakRows(layerProps);
      index.push({ memberId: m.id, yBand: yBand as [number, number], peaks });
    }
    if (showCrossTraceTracking) {
      const trackingMarks = buildCrossTraceTrackingMarks({
        members,
        yBands: yBands as Array<[number, number]>,
      });
      for (const mk of trackingMarks) marks.push(mk);
    }
    return { allMarks: marks as Plot.Markish[], hoverPeakIndex: index };
  }, [
    members, traces, xDomain, xType,
    peakDisplayByMemberId, highlightedMemberId,
    showPeakTicks, showPeakLabels, groupingMode, sampleIdFor,
    workingBandFraction, representation, showCrossTraceTracking,
    qExtent, panelWidth, panelHeight, _resizeKey,
  ]);

  // Per-band peak click hit-test (edit mode). Routes through PlotSurface's
  // onPointerClick escape hatch: picks the member whose y-band contains the
  // click, then the closest peak within PEAK_HIT_PX. Returns `true` when a
  // peak is hit so PlotSurface skips its own (single-id) peak/add handlers.
  const handlePointerClick = useCallback(
    (ctx: PlotOverlayContext, clickX: number, clickY: number, altKey: boolean): boolean => {
      if (!onPeakClick) return false;
      let best: { memberId: number; peakId: number; dist: number } | null = null;
      for (const band of hoverPeakIndex) {
        const [top, bottom] = band.yBand;
        if (clickY < top - PEAK_HIT_PX || clickY > bottom + PEAK_HIT_PX) continue;
        for (const p of band.peaks) {
          const px = ctx.applyQ(p.q);
          if (px === null) continue;
          const d = Math.abs(px - clickX);
          if (d <= PEAK_HIT_PX && (best === null || d < best.dist)) {
            best = { memberId: band.memberId, peakId: p.peakId, dist: d };
          }
        }
      }
      if (best !== null) {
        onPeakClick(best.memberId, best.peakId, altKey);
        return true;
      }
      return false;
    },
    [onPeakClick, hoverPeakIndex],
  );

  // Peak hover tooltip (unconditional — informational). Mirrors the legacy
  // handleHoverMove: hit-test against ALL members, translate to host-relative
  // coords for the absolutely-positioned tooltip.
  const handlePointerMove = useCallback(
    (ctx: PlotOverlayContext, cursorX: number, cursorY: number): void => {
      let best: { memberId: number; peakId: number; q: number; dist: number } | null = null;
      for (const band of hoverPeakIndex) {
        const [top, bottom] = band.yBand;
        if (cursorY < top - PEAK_HIT_PX || cursorY > bottom + PEAK_HIT_PX) continue;
        for (const p of band.peaks) {
          const px = ctx.applyQ(p.q);
          if (px === null) continue;
          const d = Math.abs(px - cursorX);
          if (d <= PEAK_HIT_PX && (best === null || d < best.dist)) {
            best = { memberId: band.memberId, peakId: p.peakId, q: p.q, dist: d };
          }
        }
      }
      if (best === null) {
        setTooltip(null);
        return;
      }
      // cursorX/Y are plotContainer-relative; the tooltip is positioned against
      // hostRef. After the #81 width cap the two origins differ, so translate.
      const host = hostRef.current;
      const cont = plotContainer.current;
      const hostRect = host?.getBoundingClientRect();
      const contRect = cont?.getBoundingClientRect();
      const dx = (contRect?.left ?? 0) - (hostRect?.left ?? 0);
      const dy = (contRect?.top ?? 0) - (hostRect?.top ?? 0);
      setTooltip({
        q: best.q,
        peakId: best.peakId,
        xPx: cursorX + dx,
        yPx: cursorY + dy,
      });
    },
    [hoverPeakIndex],
  );

  const handlePointerLeave = useCallback(() => setTooltip(null), []);

  const margins = useMemo(
    () => ({ left: MARGIN_LEFT, right: MARGIN_RIGHT, top: MARGIN_TOP, bottom: MARGIN_BOTTOM }),
    [],
  );
  const ext = qExtent();


  // Per-member invisible overlays carrying `data-testid="member-trace"` and
  // `data-member-id={id}` for E2E selectors (Plan §"E2E selector and
  // accessibility strategy"). MemberTraceLayer is a mark factory that
  // returns null, so the per-band selector cannot live on a JSX element it
  // owns. The overlays mount at each band's y-position with
  // `pointer-events: none` so they don't interfere with the plot's own
  // wheel/brush listeners. Phase 9's hover-driven phase coloring will hang
  // hover affordances on these same nodes.
  // Always emit the per-member overlays (one per member) so E2E selectors
  // resolve regardless of layout state. Positions degenerate to [0, 0] when
  // panelHeight is 0 (initial mount before ResizeObserver fires, JSDOM
  // tests) — that's fine; the boxes still exist and carry the right
  // `data-member-id`. Once layout settles, the next render places them at
  // the correct y-band envelope.
  const overlayBands = computeYBands(members.map((m) => m.band_height || 1), panelHeight);

  // E-2/E-3: the (phase,order) anchor map for the migration-tracking overlay.
  // Derived from member snapshots only — bounded by member count, recomputed
  // when the snapshot set changes. Gated on `showCrossTraceTracking` so the
  // overlay (and its hover handles) mount only when tracking is enabled.
  const anchorMap = useMemo(() => buildAnchorMap(members), [members]);

  // Issue #81 — self-constrain the plot width so each member's band lands at
  // a square aspect ratio (default 1.0). Without this, 4 members in an
  // ~810 px-wide column produces ~5.4:1 W:H bands which visually crushes peak
  // features. The math lives in `computeMaxPlotWidth` for unit-testability —
  // see `MultiTracePlot.test.tsx` for the test cases that pin it.
  const maxPlotWidth = useMemo(
    () => computeMaxPlotWidth(panelHeight, members.length),
    [members.length, panelHeight],
  );

  return (
    <div
      ref={hostRef}
      className="w-full h-full relative"
      data-testid="multi-trace-plot"
      data-representation={representation}
    >
      {/*
        Inner wrapper carries the width cap + horizontal centering. Both the
        plot host and the per-band overlays are children of this wrapper so
        their bounding boxes coincide pixel-for-pixel. The plot host gets
        replaceChildren'd by Plot.plot() — overlays must be a SIBLING, not a
        child, or they'd get wiped on every render.
      */}
      <div
        className="h-full mx-auto relative"
        style={{ maxWidth: `${maxPlotWidth}px`, width: "100%" }}
      >
        <div ref={plotContainer} className="w-full h-full absolute inset-0">
          <PlotSurface
            marks={allMarks}
            xType={xType}
            xDomain={xDomain}
            yType="linear"
            yDomain={[0, panelHeight || 300]}
            yReversed
            hideYAxis
            margins={margins}
            {...(qUnits !== undefined ? { qUnits } : {})}
            {...(ext !== null ? { xExtent: ext } : {})}
            onXDomain={onXDomain}
            onReset={() => onXDomain(null)}
            onBrush={(d) => onXDomain(d)}
            {...(onPeakClick ? { onPointerClick: handlePointerClick } : {})}
            onPointerMove={handlePointerMove}
            onPointerLeave={handlePointerLeave}
            data-testid="multi-trace-plot-surface"
            {...(showCrossTraceTracking
              ? {
                  overlay: (ctx: PlotOverlayContext) => (
                    <SeriesTrackingOverlay
                      ctx={ctx}
                      members={members}
                      anchorMap={anchorMap}
                      yBands={overlayBands as Array<[number, number]>}
                      trackedKey={trackedKey}
                      onTrack={setTrackedKey}
                    />
                  ),
                }
              : {})}
          />
        </div>
        <div
          aria-hidden="true"
          className="absolute inset-0 pointer-events-none"
          data-testid="member-trace-overlays"
        >
          {members.map((m, i) => {
            const band = overlayBands[i];
            const top = band ? band[0] : 0;
            const height = band ? band[1] - band[0] : 0;
            return (
              <MemberBandOverlay
                key={m.id}
                memberId={m.id}
                top={top}
                height={height}
              />
            );
          })}
        </div>
      </div>
      {tooltip !== null && (
        <div
          data-testid="peak-tooltip"
          role="tooltip"
          className="absolute pointer-events-none rounded border border-hair
                     bg-plate text-ink text-xs px-2 py-1 shadow"
          style={{
            // Offset above-and-right of the cursor so the tooltip doesn't
            // sit under the pointer (which would block subsequent clicks).
            left: `${tooltip.xPx + 8}px`,
            top: `${tooltip.yPx - 22}px`,
            zIndex: 10,
          }}
        >
          q = {tooltip.q.toPrecision(3)}
          {showPeakIds ? <span className="text-ink-soft ml-2">id={tooltip.peakId}</span> : null}
        </div>
      )}
    </div>
  );
}

/**
 * Per-member invisible plot overlay (Plan §"E2E selector and accessibility
 * strategy"). Carries `data-testid="member-trace"` + `data-member-id` for
 * E2E selectors, and subscribes to `ActiveBandContext` for two distinct
 * gutter-driven gestures:
 *
 *   - Compare UX E-4 — a hovered/dragged inter-row resize gap tints the
 *     band ABOVE it: the overlay gains `data-active-band` when its member
 *     id is the published active band.
 *   - Compare UX E-5 — a pointer-drag reorder highlights the band the
 *     dragged row would land in: the overlay gains `data-drop-target` when
 *     its member id is the published drop target.
 *
 * Outside an `ActiveBandProvider` the context returns `null` for both, so
 * neither attribute is ever set.
 */
function MemberBandOverlay(props: {
  memberId: number;
  top: number;
  height: number;
}): JSX.Element {
  const { activeBandMemberId, dropTargetMemberId } = useActiveBand();
  const active = activeBandMemberId === props.memberId;
  const dropTarget = dropTargetMemberId === props.memberId;
  return (
    <div
      data-testid="member-trace"
      data-member-id={String(props.memberId)}
      {...(active ? { "data-active-band": String(props.memberId) } : {})}
      {...(dropTarget ? { "data-drop-target": "true" } : {})}
      className={
        dropTarget
          ? "ring-1 ring-inset ring-print-accent bg-print-accent/15"
          : active
            ? "bg-print-accent/10"
            : undefined
      }
      style={{
        position: "absolute",
        left: 0,
        right: 0,
        top: `${props.top}px`,
        height: `${props.height}px`,
      }}
    />
  );
}

// Silence unused-Scale warning while keeping the type next to the helpers
// for future overlay work (cursor crosshair, peak hover effects).
export type { Scale };
