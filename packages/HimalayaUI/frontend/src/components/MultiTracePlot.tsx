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
import { useCallback, useEffect, useRef, useState } from "react";
import * as Plot from "@observablehq/plot";
import type { Trace, ComparisonMember } from "../api";
import { buildMemberMarks } from "./MemberTraceLayer";
import { invertQ } from "../lib/plot/invertQ";
import { prettifyUnits } from "../lib/units";

/** Hardcoded plot aspect ratio (W / H) per spec §Plot rendering. */
export const COMPARE_PLOT_ASPECT = 0.3;

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
}

export function MultiTracePlot(props: MultiTracePlotProps): JSX.Element {
  const {
    members, traces, xDomain, onXDomain,
    peakDisplayByMemberId, highlightedMemberId,
    qUnits, xType = "log",
  } = props;

  const hostRef       = useRef<HTMLDivElement>(null);
  const plotContainer = useRef<HTMLDivElement>(null);
  const plotElRef     = useRef<HTMLElement | SVGElement | null>(null);

  const [_resizeKey, setResizeKey] = useState(0);
  // Tracked panel height drives the per-band overlay positions in the JSX
  // below. Read from the same `clientHeight` source the plot uses, so the
  // overlay y-bands stay in sync with the rendered plot.
  const [panelHeight, setPanelHeight] = useState(0);

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

    const allMarks: unknown[] = [];
    for (let i = 0; i < members.length; i++) {
      const m = members[i]!;
      const yBand = yBands[i] ?? [0, panelH];
      const trace = m.exposure_id !== null ? traces.get(m.exposure_id) : undefined;
      const peakDisplay = peakDisplayByMemberId?.get(m.id);
      const memberMarks = buildMemberMarks({
        member: m,
        trace,
        yBand,
        peakDisplay,
        highlightedIndexId:
          highlightedMemberId === m.id && m.snapshot?.confirmed_index
            ? m.snapshot.confirmed_index.id
            : undefined,
      });
      for (const mk of memberMarks) allMarks.push(mk);
    }

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
        ...(xDomain ? { domain: xDomain } : {}),
      },
      // Y-axis is in pixel-space envelope coordinates produced by
      // `applyNormalization`. Using a linear identity scale on
      // [0, panelH] keeps the rendered y values aligned with the y-bands.
      y: {
        type: "linear",
        domain: [0, panelH],
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
      container.replaceChildren();
      plotElRef.current = null;
    };
  }, [
    members, traces, xDomain, xType, qUnits,
    peakDisplayByMemberId, highlightedMemberId,
    onXDomain, qExtent,
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

  return (
    <div
      ref={hostRef}
      className="w-full h-full relative"
      data-testid="multi-trace-plot"
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
  );
}

// Silence unused-Scale warning while keeping the type next to the helpers
// for future overlay work (cursor crosshair, peak hover effects).
export type { Scale };
