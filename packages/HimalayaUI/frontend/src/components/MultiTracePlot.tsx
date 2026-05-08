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
import type { Trace, ComparisonMember } from "../api";
import { buildMemberMarks, buildMemberPeakRows } from "./MemberTraceLayer";
import type { PeakRow } from "./MemberTraceLayer";
import { invertQ, applyQ } from "../lib/plot/invertQ";
import { formatAxis } from "../lib/plot/formatAxis";
import { prettifyUnits } from "../lib/units";
import type { GroupingMode } from "../lib/comparison/coloring";

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
 * Compute the per-member y-band envelopes (top, bottom in pixels) given
 * the member band ratios and the panel pixel height. Pure function; the
 * order of `bandRatios` matches the caller's render order (display_order).
 */
export function computeYBands(
  bandRatios: number[],
  panelHeight: number,
): Array<[number, number]> {
  if (bandRatios.length === 0) return [];
  const total = bandRatios.reduce((s, r) => s + Math.max(0, r), 0);
  if (total <= 0) {
    // Degenerate: every ratio is zero. Fall back to equal slicing so the
    // plot still has well-defined bands.
    const each = panelHeight / bandRatios.length;
    return bandRatios.map((_, i) => [i * each, (i + 1) * each]);
  }
  const out: Array<[number, number]> = [];
  let cumulative = 0;
  for (const r of bandRatios) {
    const top    = (cumulative / total) * panelHeight;
    const next   = cumulative + Math.max(0, r);
    const bottom = (next / total) * panelHeight;
    out.push([top, bottom]);
    cumulative = next;
  }
  return out;
}

export interface MultiTracePlotProps {
  /** Members in render order (top → bottom). Caller sorts by `display_order`. */
  members: ComparisonMember[];
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
   * "color_override → var(--color-fg)" behaviour (e.g. tests that don't
   * exercise grouping).
   */
  groupingMode?: GroupingMode;
  /**
   * Sample-id resolver supplied by the parent page (typically wired against
   * the TanStack exposure cache). `colorFor` consults this in `bySample`
   * mode. Required when `groupingMode` is set; ignored otherwise.
   */
  sampleIdFor?: (m: ComparisonMember) => number | null;
}

export function MultiTracePlot(props: MultiTracePlotProps): JSX.Element {
  const {
    members, traces, xDomain, onXDomain,
    peakDisplayByMemberId, highlightedMemberId,
    qUnits, xType = "log",
    onPeakClick,
    showPeakTicks = true, showPeakLabels = true,
    groupingMode, sampleIdFor,
  } = props;

  const hostRef       = useRef<HTMLDivElement>(null);
  const plotContainer = useRef<HTMLDivElement>(null);
  const plotElRef     = useRef<HTMLElement | SVGElement | null>(null);

  const [_resizeKey, setResizeKey] = useState(0);
  // Tracked panel height drives the per-band overlay positions in the JSX
  // below. Read from the same `clientHeight` source the plot uses, so the
  // overlay y-bands stay in sync with the rendered plot.
  const [panelHeight, setPanelHeight] = useState(0);

  // Phase 8.3 — peak hover tooltip state. `null` = nothing hovered. The
  // `xPx` / `yPx` values are container-local pixel coordinates used to
  // position the tooltip overlay; `q` and `peakId` are the rendered fields.
  // Peak id is only displayed when the developer-only `?showPeakIds` URL
  // flag is set (read once on mount; doesn't react to history changes).
  const [tooltip, setTooltip] = useState<{
    q: number; peakId: number; xPx: number; yPx: number;
  } | null>(null);
  const showPeakIds = typeof window !== "undefined"
    && new URLSearchParams(window.location.search).has("showPeakIds");

  useEffect(() => {
    const el = plotContainer.current;
    if (!el) return;
    setPanelHeight(el.clientHeight);
    const obs = new ResizeObserver(() => {
      setResizeKey((k) => k + 1);
      if (plotContainer.current) setPanelHeight(plotContainer.current.clientHeight);
    });
    obs.observe(el);
    return () => obs.disconnect();
  }, []);

  // Derive q-domain from union of all member traces (for full-range default).
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

  // Imperative render-and-bind: builds the Plot element from current props,
  // installs the wheel/dblclick/brush listeners, and returns a cleanup that
  // detaches them. Wrapped in `useCallback` with its true deps so the effect
  // below depends on `[renderPlot, _resizeKey]` alone — no eslint-disable,
  // no hand-curated dep list (per CLAUDE.md "Imperative render functions
  // in effects: use `useCallback`").
  const renderPlot = useCallback((): (() => void) | undefined => {
    const container = plotContainer.current;
    if (!container) return;

    const panelW = container.clientWidth  || 400;
    const panelH = container.clientHeight || 300;

    const ratios = members.map((m) => m.band_height || 1);
    const yBands = computeYBands(ratios, panelH);

    // Synthesize an x-scale directly from the domain + plot width so the
    // label-dodge layout (Phase 8.2) can run in the SAME pass as mark
    // building — no two-pass `Plot.plot()` call needed. The plot interior
    // is `[MARGIN_LEFT, panelW - MARGIN_RIGHT]`; pixels outside that range
    // map to the corresponding domain edges (matches Plot's clamping).
    const ext = qExtent();
    const xMin = xDomain ? xDomain[0] : ext?.[0] ?? 1;
    const xMax = xDomain ? xDomain[1] : ext?.[1] ?? 10;
    const innerW = Math.max(1, panelW - MARGIN_LEFT - MARGIN_RIGHT);
    const xScale = makeXScale(xType, xMin, xMax, MARGIN_LEFT, innerW);

    const allMarks: unknown[] = [];
    // Per-member visible peak rows + y-band, captured for click + hover
    // hit-testing. Built in the same loop as the marks so we never drift
    // from what was rendered (the same `buildMemberPeakRows` call sources
    // both). The HOVER index is unconditional (tooltip works in review
    // mode too); the CLICK index is the same data, gated by `onPeakClick`.
    const hoverPeakIndex: Array<{
      memberId: number;
      yBand: [number, number];
      peaks: PeakRow[];
    }> = [];
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
        // Phase 9 gap-fix: forward grouping context so the per-member line
        // stroke recolors with the toggle. `MemberTraceLayer` falls back
        // to legacy single-color rendering when these are undefined.
        ...(groupingMode !== undefined && sampleIdFor !== undefined
          ? { groupingMode, allMembers: members, sampleIdFor }
          : {}),
      };
      const memberMarks = buildMemberMarks(layerProps);
      for (const mk of memberMarks) allMarks.push(mk);
      const { peaks } = buildMemberPeakRows(layerProps);
      hoverPeakIndex.push({ memberId: m.id, yBand: yBand as [number, number], peaks });
    }
    // Click-hit-test reuses the same index when `onPeakClick` is wired.
    const peakIndex = onPeakClick ? hoverPeakIndex : [];

    const el = Plot.plot({
      width:  panelW,
      height: panelH,
      marginLeft: MARGIN_LEFT, marginRight: MARGIN_RIGHT,
      marginTop: MARGIN_TOP,  marginBottom: MARGIN_BOTTOM,
      style: {
        fontFamily: "var(--font-sans)",
        color: "var(--color-fg-muted)",
        background: "transparent",
        overflow: "visible",
      },
      x: {
        type: xType,
        label: `q (${prettifyUnits(qUnits ?? "A-1")})`,
        // Plain decimal tick labels — Plot's default SI-suffix formatter
        // renders 0.040 as "40m" which is unhelpful for SAXS q values.
        // Shared with `TraceViewer` for cross-page parity (issue #80).
        tickFormat: (d: number) => formatAxis(d),
        ...(xDomain ? { domain: xDomain } : {}),
      },
      // Y-axis is in pixel-space envelope coordinates produced by
      // `applyNormalization` (small y = top of band, large y = bottom).
      // Plot's default orientation maps domain[0] → bottom, domain[1] →
      // top, so we reverse the domain to keep small-y at top — otherwise
      // every trace renders upside-down (issue #63).
      y: {
        type: "linear",
        domain: [panelH, 0],
        // Hide y axis — the y-band layout is the visualization, not the
        // numbers themselves.
        axis: null,
      },
      // Observable Plot's `Markish` type is closed over the public mark
      // factories; `buildMemberMarks` returns `unknown[]` because callers
      // shouldn't depend on the internal mark constructor shapes. Cast at
      // the boundary; the runtime contract holds (these were produced by
      // `Plot.line / Plot.dot / Plot.text`).
      marks: allMarks as Plot.Markish[],
    });

    container.replaceChildren(el);
    plotElRef.current = el as unknown as HTMLElement;

    function handleWheel(evRaw: Event): void {
      const ev = evRaw as WheelEvent;
      ev.preventDefault();
      const rect = container!.getBoundingClientRect();
      const cursorQ = invertQ(plotElRef.current, ev.clientX - rect.left);
      if (cursorQ === null) return;
      const ext = qExtent();
      const curMin = xDomain ? xDomain[0] : ext?.[0] ?? cursorQ * 0.5;
      const curMax = xDomain ? xDomain[1] : ext?.[1] ?? cursorQ * 2;
      const factor = Math.exp(ev.deltaY * 0.001);
      const q0 = ext?.[0] ?? curMin;
      const qN = ext?.[1] ?? curMax;
      let newMin: number;
      let newMax: number;
      if (xType === "log") {
        const logMin = Math.log(curMin);
        const logMax = Math.log(curMax);
        const logCur = Math.log(Math.max(cursorQ, 1e-6));
        newMin = Math.max(q0, Math.exp(logCur - (logCur - logMin) * factor));
        newMax = Math.min(qN, Math.exp(logCur + (logMax - logCur) * factor));
      } else {
        newMin = Math.max(q0, cursorQ - (cursorQ - curMin) * factor);
        newMax = Math.min(qN, cursorQ + (curMax - cursorQ) * factor);
      }
      if (newMax - newMin < (qN - q0) * 1e-4) return;
      onXDomain([newMin, newMax]);
    }
    (el as unknown as EventTarget).addEventListener("wheel", handleWheel, { passive: false } as AddEventListenerOptions);

    function handleDblClick(): void {
      onXDomain(null);
    }
    (el as unknown as EventTarget).addEventListener("dblclick", handleDblClick);

    // ── peak click hit-testing (Phase 8.1, edit mode only) ─────────────
    // Only installs when `onPeakClick` is provided. Hit-test logic:
    //   1. Convert click pixel → q via the x-scale.
    //   2. Identify the member whose y-band contains the click Y.
    //   3. Find the closest peak in that member's `peaks` within
    //      `PEAK_HIT_PX` pixels (using `applyQ` to convert peak q → px).
    //   4. Dispatch `onPeakClick(memberId, peakId, altKey)`.
    // Bails silently when nothing matches — leaves brush/zoom untouched.
    function handlePeakClick(evRaw: Event): void {
      if (!onPeakClick) return;
      const ev = evRaw as MouseEvent;
      const rect = container!.getBoundingClientRect();
      const clickX = ev.clientX - rect.left;
      const clickY = ev.clientY - rect.top;
      let best: { memberId: number; peakId: number; dist: number } | null = null;
      for (const band of peakIndex) {
        // Y-band selection — pick the member whose band contains the click Y.
        // We still allow a small tolerance so clicks on the very edge of the
        // band (e.g. on a triangle that pokes above the line) hit the right
        // member.
        const [top, bottom] = band.yBand;
        if (clickY < top - PEAK_HIT_PX || clickY > bottom + PEAK_HIT_PX) continue;
        for (const p of band.peaks) {
          const px = applyQ(plotElRef.current, p.q);
          if (px === null) continue;
          const d = Math.abs(px - clickX);
          if (d <= PEAK_HIT_PX && (best === null || d < best.dist)) {
            best = { memberId: band.memberId, peakId: p.peakId, dist: d };
          }
        }
      }
      if (best !== null) onPeakClick(best.memberId, best.peakId, ev.altKey);
    }
    (el as unknown as EventTarget).addEventListener("click", handlePeakClick);

    // ── peak hover tooltip (Phase 8.3) ─────────────────────────────────
    // Mousemove hit-tests against ALL members' visible peaks (regardless
    // of edit mode — tooltip is informational and unconditional). Same
    // hit radius as click. Mouseleave hides the tooltip.
    function handleHoverMove(evRaw: Event): void {
      const ev = evRaw as MouseEvent;
      const rect = container!.getBoundingClientRect();
      const cursorX = ev.clientX - rect.left;
      const cursorY = ev.clientY - rect.top;
      let best: { memberId: number; peakId: number; q: number; dist: number } | null = null;
      for (const band of hoverPeakIndex) {
        const [top, bottom] = band.yBand;
        if (cursorY < top - PEAK_HIT_PX || cursorY > bottom + PEAK_HIT_PX) continue;
        for (const p of band.peaks) {
          const px = applyQ(plotElRef.current, p.q);
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
      // Tooltip is positioned absolutely against `hostRef`, but `cursorX/Y`
      // above are `plotContainer`-relative (used for plot-pixel hit-testing
      // via `applyQ`). After issue #81's width cap, `plotContainer` is
      // `mx-auto`-centered inside `hostRef`, so the two origins differ by
      // `(hostWidth - maxPlotWidth) / 2`. Translate to host-relative before
      // storing — otherwise the tooltip lands offset to the left of the
      // hovered peak. `hostRef.current` is non-null whenever `plotContainer`
      // is (plotContainer is a descendant of hostRef in the JSX), so the
      // bang matches the `container!` usage above.
      const hostRect = hostRef.current!.getBoundingClientRect();
      setTooltip({
        q: best.q,
        peakId: best.peakId,
        xPx: ev.clientX - hostRect.left,
        yPx: ev.clientY - hostRect.top,
      });
    }
    function handleHoverLeave(): void {
      setTooltip(null);
    }
    (el as unknown as EventTarget).addEventListener("mousemove", handleHoverMove);
    (el as unknown as EventTarget).addEventListener("mouseleave", handleHoverLeave);

    // Brush-to-zoom: drag horizontally to set a q sub-range. Implemented as
    // mousedown→mousemove→mouseup; we track pixel coords and invert at end.
    let brushStartPx: number | null = null;
    function handleMouseDown(evRaw: Event): void {
      const ev = evRaw as MouseEvent;
      const rect = container!.getBoundingClientRect();
      brushStartPx = ev.clientX - rect.left;
    }
    function handleMouseUp(evRaw: Event): void {
      if (brushStartPx === null) return;
      const ev = evRaw as MouseEvent;
      const rect = container!.getBoundingClientRect();
      const endPx = ev.clientX - rect.left;
      const start = brushStartPx;
      brushStartPx = null;
      // Ignore tiny drags (single-click).
      if (Math.abs(endPx - start) < 4) return;
      const a = invertQ(plotElRef.current, Math.min(start, endPx));
      const b = invertQ(plotElRef.current, Math.max(start, endPx));
      if (a === null || b === null) return;
      if (b <= a) return;
      onXDomain([a, b]);
    }
    (el as unknown as EventTarget).addEventListener("mousedown", handleMouseDown);
    (el as unknown as EventTarget).addEventListener("mouseup", handleMouseUp);

    return () => {
      (el as unknown as EventTarget).removeEventListener("wheel", handleWheel);
      (el as unknown as EventTarget).removeEventListener("dblclick", handleDblClick);
      (el as unknown as EventTarget).removeEventListener("mousedown", handleMouseDown);
      (el as unknown as EventTarget).removeEventListener("mouseup", handleMouseUp);
      (el as unknown as EventTarget).removeEventListener("click", handlePeakClick);
      (el as unknown as EventTarget).removeEventListener("mousemove", handleHoverMove);
      (el as unknown as EventTarget).removeEventListener("mouseleave", handleHoverLeave);
      container.replaceChildren();
      plotElRef.current = null;
      setTooltip(null);
    };
  }, [
    members, traces, xDomain, xType, qUnits,
    peakDisplayByMemberId, highlightedMemberId,
    onXDomain, qExtent, onPeakClick,
    showPeakTicks, showPeakLabels,
    groupingMode, sampleIdFor,
  ]);

  useEffect(() => {
    return renderPlot();
    // `_resizeKey` rerenders the plot when the container resizes — the
    // ResizeObserver bumps the key so we recompute panelW/panelH. It's not
    // captured by `renderPlot` (we read `container.clientWidth` directly),
    // so include it as a primitive dep.
  }, [renderPlot, _resizeKey]);

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
        <div ref={plotContainer} className="w-full h-full" />
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
              <div
                key={m.id}
                data-testid="member-trace"
                data-member-id={String(m.id)}
                style={{
                  position: "absolute",
                  left: 0,
                  right: 0,
                  top: `${top}px`,
                  height: `${height}px`,
                }}
              />
            );
          })}
        </div>
      </div>
      {tooltip !== null && (
        <div
          data-testid="peak-tooltip"
          role="tooltip"
          className="absolute pointer-events-none rounded border border-border
                     bg-bg-elevated text-fg text-xs px-2 py-1 shadow"
          style={{
            // Offset above-and-right of the cursor so the tooltip doesn't
            // sit under the pointer (which would block subsequent clicks).
            left: `${tooltip.xPx + 8}px`,
            top: `${tooltip.yPx - 22}px`,
            zIndex: 10,
          }}
        >
          q = {tooltip.q.toPrecision(3)}
          {showPeakIds ? <span className="text-fg-muted ml-2">id={tooltip.peakId}</span> : null}
        </div>
      )}
    </div>
  );
}

// Silence unused-Scale warning while keeping the type next to the helpers
// for future overlay work (cursor crosshair, peak hover effects).
export type { Scale };
