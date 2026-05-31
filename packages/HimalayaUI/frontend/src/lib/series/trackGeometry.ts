/**
 * series/trackGeometry — turn the `(phase,order)` anchor map into renderable
 * migration-track geometry in (q, y) DATA space (Plan E, Task E-2 overlay half).
 *
 * The SVG overlay in `MultiTracePlot` projects these (q, y) coordinates through
 * PlotSurface's `applyQ`/`applyY` so the track + ghost rings live on the same
 * scales as the traces. Geometry-only here (no DOM) → unit-testable.
 *
 * Each track is one polyline per `(phase, order)`, threading the SAME
 * reflection across every carrying member at that member's BASELINE (the band
 * bottom, where peak ticks plant their root). Where a member carries the phase
 * but the order is predicted-but-absent, the vertex routes through the
 * predicted q (so the line doesn't kink to zero) AND a hollow ghost ring is
 * recorded there, which the overlay paints as a `<PeakGlyph>` caret
 * (predictedAbsent). A track needs ≥ 2 carrying members to read as a line.
 */
import type { SeriesMember } from "../../api";
import { phaseColor } from "../../phases";
import type { AnchorMap } from "./anchors";

/** A track vertex in (q, y) data space. `absent` vertices route through the
 *  predicted q and carry the ghost ring. */
export interface TrackPoint {
  q: number;
  y: number;
  absent: boolean;
}

/** A hollow ghost-ring marker (predicted-but-absent order) in (q, y) space. */
export interface GhostRing {
  q: number;
  y: number;
}

export interface SeriesTrack {
  /** `"phase:order"`. */
  key: string;
  phase: string;
  order: number;
  /** Resolved phase colour (the migration connector is a phase relation). */
  color: string;
  /** Polyline vertices in member order (top→bottom of the stack as passed). */
  points: TrackPoint[];
  /** Ghost rings at absent orders. */
  ghostRings: GhostRing[];
}

export interface SeriesTrackGeometry {
  tracks: SeriesTrack[];
}

/**
 * Build the migration-track geometry from the anchor map + per-member y-bands.
 * `yBands[pos]` is `[topPx, bottomPx]`; the vertex y is the band BOTTOM
 * (baseline) — matching the cross-trace polyline marks and where peak ticks
 * root in the waterfall layout.
 */
export function buildSeriesTrackGeometry(
  anchorMap: AnchorMap,
  members: SeriesMember[],
  yBands: Array<[number, number]>,
): SeriesTrackGeometry {
  // `members` is part of the contract (callers pass the same render-ordered
  // array used to build the anchor map + y-bands) but the geometry is fully
  // determined by the anchor map's member positions + the y-bands.
  void members;
  const tracks: SeriesTrack[] = [];

  for (const [key, vertices] of anchorMap) {
    // A polyline needs ≥ 2 carrying members.
    if (vertices.length < 2) continue;
    const sep = key.lastIndexOf(":");
    const phase = key.slice(0, sep);
    const order = Number(key.slice(sep + 1));
    const color = phaseColor(phase);

    const points: TrackPoint[] = [];
    const ghostRings: GhostRing[] = [];
    for (const v of vertices) {
      const band = yBands[v.memberPos];
      if (!band) continue;
      const y = band[1]; // baseline (band bottom)
      // Absent → route through predicted q; present → observed q.
      const q = v.absent ? v.predictedQ : v.q!;
      points.push({ q, y, absent: v.absent });
      if (v.absent) ghostRings.push({ q, y });
    }
    if (points.length < 2) continue;
    // Drop an all-absent track: with no observed reflection there's nothing to
    // migrate *between* — it would be a connector of pure predictions (noise).
    if (points.every((p) => p.absent)) continue;
    tracks.push({ key, phase, order, color, points, ghostRings });
  }

  return tracks.length > 0 ? { tracks } : { tracks: [] };
}
