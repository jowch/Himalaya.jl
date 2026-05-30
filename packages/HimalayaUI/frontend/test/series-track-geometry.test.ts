/**
 * series track geometry (Plan E, Task E-2 — the overlay half).
 *
 * `buildSeriesTrackGeometry` turns the `(phase,order)` anchor map + per-member
 * y-bands into renderable geometry in (q, y) DATA space (the overlay projects
 * via PlotSurface's applyQ/applyY). Each track is a polyline whose vertices
 * route through the predicted-q where an order is absent; absent vertices also
 * carry a ghost-ring marker so the overlay can plant a hollow caret there.
 */
import { describe, it, expect } from "vitest";
import type { SeriesMember } from "../src/api";
import { buildAnchorMap } from "../src/lib/series/anchors";
import { buildSeriesTrackGeometry } from "../src/lib/series/trackGeometry";

function member(
  id: number,
  phase: string | null,
  lattice_d: number,
  observedQs: number[],
): SeriesMember {
  const peakIds = observedQs.map((_, i) => id * 100 + i);
  return {
    id, series_id: 1, exposure_id: id * 10, display_order: id - 1,
    band_height: 1, y_offset: 0, normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: phase
      ? {
          effective_peaks: observedQs.map((q, i) => ({
            id: peakIds[i]!, q, intensity: 50, sharpness: 1, source: "auto" as const,
          })),
          confirmed_index: {
            id: id * 1000, phase, lattice_d, r_squared: 0.99, ngc: -1.5,
            peak_ids: peakIds,
          },
          analysis_inputs_hash: "h",
        }
      : { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
    is_stale: false, created_by: null, created_at: null,
  };
}

const a = 200;
const qOf = (rad: number) => (2 * Math.PI * Math.sqrt(rad)) / a;

describe("buildSeriesTrackGeometry", () => {
  it("emits one track per (phase,order) with vertices at member baselines", () => {
    const m1 = member(1, "Pn3m", a, [qOf(2), qOf(3)]);
    const m2 = member(2, "Pn3m", a, [qOf(2), qOf(3)]);
    const yBands: Array<[number, number]> = [[0, 50], [50, 100]];
    const map = buildAnchorMap([m1, m2]);
    const geom = buildSeriesTrackGeometry(map, [m1, m2], yBands);
    const t0 = geom.tracks.find((t) => t.key === "Pn3m:0")!;
    expect(t0.points.map((p) => p.y)).toEqual([50, 100]); // baselines (band bottoms)
    expect(t0.points.every((p) => !p.absent)).toBe(true);
    expect(t0.ghostRings).toHaveLength(0);
    expect(t0.color).toBeTruthy();
  });

  it("routes a track through the predicted-q at an absent order and records a ghost ring there", () => {
    // m1 absents order 1 (√3); m2 observes it.
    const m1 = member(1, "Pn3m", a, [qOf(2), qOf(4)]);
    const m2 = member(2, "Pn3m", a, [qOf(2), qOf(3), qOf(4)]);
    const yBands: Array<[number, number]> = [[0, 50], [50, 100]];
    const map = buildAnchorMap([m1, m2]);
    const geom = buildSeriesTrackGeometry(map, [m1, m2], yBands);
    const t1 = geom.tracks.find((t) => t.key === "Pn3m:1")!;
    // Both members contribute a vertex; the absent one routes through predicted-q.
    const v0 = t1.points.find((p) => p.y === 50)!;
    expect(v0.absent).toBe(true);
    expect(v0.q).toBeCloseTo(qOf(3), 6); // routed through predicted √3
    // One ghost ring planted at the absent vertex.
    expect(t1.ghostRings).toHaveLength(1);
    expect(t1.ghostRings[0]!.q).toBeCloseTo(qOf(3), 6);
    expect(t1.ghostRings[0]!.y).toBe(50);
  });

  it("drops a track with fewer than two carrying members", () => {
    const m1 = member(1, "Pn3m", a, [qOf(2)]);
    const m2 = member(2, "Lamellar", 60, [(2 * Math.PI) / 60]);
    const yBands: Array<[number, number]> = [[0, 50], [50, 100]];
    const map = buildAnchorMap([m1, m2]);
    const geom = buildSeriesTrackGeometry(map, [m1, m2], yBands);
    // Only single-carrier tracks → no polylines.
    expect(geom.tracks).toHaveLength(0);
  });
});
