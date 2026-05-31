/**
 * SeriesTrackingOverlay — the SVG migration-tracking layer for the Series
 * waterfall (Plan E, Tasks E-2 + E-3).
 *
 * Rendered inside `PlotSurface`'s `overlay` <svg> so it shares the live q/y
 * scales (projected via `ctx.applyQ` / `ctx.applyY`) AND so it can paint the
 * hollow predicted-absent CARET — which `Plot.dot` cannot draw, but the SVG
 * `<PeakGlyph>` path (peakMark.ts §5.1 atoms) can.
 *
 * Two jobs:
 *   E-2  — for each `(phase, order)` carried by ≥ 2 members, draw a terracotta
 *          migration connector threading the SAME reflection across the stack
 *          (routing through the predicted q where an order is absent) and a
 *          hollow ghost-ring caret at each absent order.
 *   E-3  — light plate-ringed anchor beads at every OBSERVED reflection are
 *          hover/focus targets: pointing at one tracks that `(phase, order)`
 *          across the whole stack (the connector + ghost rings light up). The
 *          tracked state is ephemeral, owned by the parent (`trackedKey` +
 *          `onTrack`). Keyboard focus drives the same highlight (a11y P1).
 *
 * The connector is terracotta (the migration is an *interaction-driven*
 * cross-highlight the user is dragging across the stack — the grease-pencil
 * accent, matching the mockup's `.wf-track-line`), distinct from the phase-hued
 * always-on `CrossTraceTrackingLayer` polylines. Anchor beads are phase-hued.
 */
import { useMemo } from "react";
import type { SeriesMember } from "../api";
import type { PlotOverlayContext } from "./PlotSurface";
import { phaseColor } from "../phases";
import { peakGlyph } from "./ui/peakMark";
import { PeakGlyph } from "./ui/PeakGlyph";
import type { AnchorMap } from "../lib/series/anchors";
import { buildSeriesTrackGeometry } from "../lib/series/trackGeometry";

export interface SeriesTrackingOverlayProps {
  ctx: PlotOverlayContext;
  /** Members in render order (aligned with `yBands`). */
  members: SeriesMember[];
  /** `(phase,order)` anchor map from `buildAnchorMap(members)`. */
  anchorMap: AnchorMap;
  /** Per-member pixel y-bands `[topPx, bottomPx]` (parent's `computeYBands`). */
  yBands: Array<[number, number]>;
  /** Currently-tracked `(phase,order)` key (hover/focus). null = none. */
  trackedKey?: string | null;
  /** Set/clear the tracked key from an anchor hover/focus. */
  onTrack?: (key: string | null) => void;
}

/** Accent (terracotta) for the active migration connector — matches the
 *  mockup `.wf-track-line`. Resolved at module scope via a CSS var token so the
 *  on-screen overlay theming stays single-sourced. */
const ACCENT = "var(--color-accent)";
const PLATE = "var(--color-plate)";

export function SeriesTrackingOverlay({
  ctx, members, anchorMap, yBands, trackedKey = null, onTrack,
}: SeriesTrackingOverlayProps): JSX.Element {
  const geom = useMemo(
    () => buildSeriesTrackGeometry(anchorMap, members, yBands),
    [anchorMap, members, yBands],
  );

  // Observed anchor beads: one per (phase,order) vertex that is NOT absent.
  // These are the hover/focus handles (E-3).
  const beads = useMemo(() => {
    const out: Array<{ key: string; q: number; y: number; color: string }> = [];
    for (const [key, vertices] of anchorMap) {
      const sep = key.lastIndexOf(":");
      const phase = key.slice(0, sep);
      const color = phaseColor(phase);
      for (const v of vertices) {
        if (v.absent || v.q === null) continue;
        const band = yBands[v.memberPos];
        if (!band) continue;
        out.push({ key, q: v.q, y: band[1], color });
      }
    }
    return out;
  }, [anchorMap, yBands]);

  return (
    <g data-role="series-tracking-overlay">
      {/* ── E-2: migration connectors + ghost rings ───────────────────── */}
      {geom.tracks.map((track) => {
        const isTracked = trackedKey === track.key;
        // Only paint the connector when this track is the hovered one (the
        // mockup shows the terracotta migration line on hover, not always-on —
        // the always-on phase polylines are the separate CrossTraceTrackingLayer).
        const pts = track.points
          .map((p) => {
            const px = ctx.applyQ(p.q);
            const py = ctx.applyY(p.y);
            return px === null || py === null ? null : `${px.toFixed(1)},${py.toFixed(1)}`;
          })
          .filter((s): s is string => s !== null);
        return (
          <g key={track.key} data-role="series-track" data-key={track.key}>
            {isTracked && pts.length >= 2 && (
              <polyline
                data-role="series-track-line"
                data-tracked="true"
                points={pts.join(" ")}
                fill="none"
                stroke={ACCENT}
                strokeWidth={1.5}
                strokeDasharray="4 3"
                strokeOpacity={0.9}
              />
            )}
            {!isTracked && pts.length >= 2 && (
              // Off-hover, the track is latent — render an invisible carrier so
              // the test (and DOM order) sees one polyline per carried track,
              // and the connector flips on instantly when hovered.
              <polyline
                data-role="series-track-line"
                data-tracked="false"
                points={pts.join(" ")}
                fill="none"
                stroke={track.color}
                strokeWidth={1}
                strokeOpacity={0}
              />
            )}
            {/* Ghost rings: a hollow predicted-absent CARET at each absent order.
                Shown when tracked (the migration line threads them) so the
                missing peak reads as missing. */}
            {isTracked && track.ghostRings.map((g, gi) => {
              const px = ctx.applyQ(g.q);
              const py = ctx.applyY(g.y);
              if (px === null || py === null) return null;
              const descriptor = peakGlyph({
                color: ACCENT, predictedAbsent: true, source: "auto", r: 4.2,
              });
              return (
                <PeakGlyph
                  key={`ghost-${gi}`}
                  descriptor={descriptor}
                  x={px}
                  y={py}
                  haloStroke={PLATE}
                />
              );
            })}
          </g>
        );
      })}

      {/* ── E-3: hover/focus anchor handles ───────────────────────────── */}
      {beads.map((b, bi) => {
        const px = ctx.applyQ(b.q);
        const py = ctx.applyY(b.y);
        if (px === null || py === null) return null;
        const tracked = trackedKey === b.key;
        return (
          <g key={`bead-${bi}`} data-role="series-anchor">
            <circle
              data-role="series-anchor-bead"
              data-key={b.key}
              cx={px}
              cy={py}
              r={tracked ? 4.2 : 2.9}
              fill={b.color}
              stroke={tracked ? ACCENT : PLATE}
              strokeWidth={tracked ? 1.8 : 1.4}
            />
            {/* Larger transparent hit target + keyboard focus path (a11y P1). */}
            <circle
              data-role="series-anchor-hit"
              data-key={b.key}
              cx={px}
              cy={py}
              r={9}
              fill="transparent"
              tabIndex={0}
              role="button"
              aria-label={`Track ${b.key.replace(":", " order ")} across the series`}
              style={{ cursor: "pointer", pointerEvents: "all" }}
              onMouseEnter={() => onTrack?.(b.key)}
              onMouseLeave={() => onTrack?.(null)}
              onFocus={() => onTrack?.(b.key)}
              onBlur={() => onTrack?.(null)}
            />
          </g>
        );
      })}
    </g>
  );
}
