/**
 * CrossTraceTrackingLayer tests (#208 — render-core finish).
 *
 * The cross-trace peak-tracking annotation links the same Miller-order
 * reflection across stacked traces. Tests pin:
 *   - no marks when fewer than two members share a phase
 *   - one Plot.line call grouped by (phase, order) when ≥ 2 carriers exist
 *   - each polyline carries the phase hue + an opacity ramp on order
 *   - vertex y is the *baseline* of each member's y-band
 *   - peak_ids ordering (Miller-order from the Julia indexer) is preserved
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import * as Plot from "@observablehq/plot";
import type { SeriesMember } from "../src/api";
import { buildCrossTraceTrackingMarks } from "../src/print/export/CrossTraceTrackingLayer";
import { phaseColor } from "../src/phases";

vi.mock("@observablehq/plot", () => ({
  line: vi.fn((data: unknown, opts: unknown) => ({ _kind: "line", data, opts })),
}));

beforeEach(() => {
  (Plot.line as unknown as { mockClear: () => void }).mockClear();
});

function member(
  id: number,
  phase: string | null,
  peakIds: number[],
  qs: number[],
): SeriesMember {
  const effective_peaks = peakIds.map((pid, i) => ({
    id: pid, q: qs[i]!, intensity: 50, sharpness: 1, source: "auto" as const,
  }));
  return {
    id, series_id: 100, exposure_id: id * 10, display_order: id - 1,
    band_height: 1, y_offset: 0, normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: phase
      ? {
          effective_peaks,
          confirmed_index: {
            id: id * 100, phase, lattice_d: 12, r_squared: 0.99, ngc: -1.5,
            peak_ids: peakIds,
          },
          analysis_inputs_hash: "h",
        }
      : null,
    is_stale: false, created_by: null, created_at: null,
  };
}

describe("buildCrossTraceTrackingMarks", () => {
  it("emits no marks when no phase is shared by two-or-more members", () => {
    const m1 = member(1, "Pn3m", [11, 12], [0.10, 0.20]);
    const m2 = member(2, "Hexagonal", [21, 22], [0.05, 0.10]);
    expect(buildCrossTraceTrackingMarks({
      members: [m1, m2],
      yBands: [[0, 50], [50, 100]],
    })).toEqual([]);
  });

  it("emits no marks when every member is unindexed", () => {
    const m1 = member(1, null, [], []);
    const m2 = member(2, null, [], []);
    expect(buildCrossTraceTrackingMarks({
      members: [m1, m2],
      yBands: [[0, 50], [50, 100]],
    })).toEqual([]);
  });

  it("emits one Plot.line mark with z-grouped polylines for shared-phase members", () => {
    const m1 = member(1, "Pn3m", [11, 12], [0.10, 0.20]);
    const m2 = member(2, "Pn3m", [21, 22], [0.11, 0.22]);
    const marks = buildCrossTraceTrackingMarks({
      members: [m1, m2],
      yBands: [[0, 50], [50, 100]],
    });
    expect(marks).toHaveLength(1);
    expect(Plot.line).toHaveBeenCalledTimes(1);

    const [data, opts] = (Plot.line as unknown as { mock: { calls: unknown[][] } }).mock.calls[0]!;
    const verts = data as Array<{ key: string; q: number; y: number; phase: string; order: number }>;
    // Two polylines (k=0, k=1) × two members each → 4 vertices.
    expect(verts).toHaveLength(4);
    // Polyline grouping uses `z: "key"` so each (phase, order) is its own line.
    expect((opts as { z: string }).z).toBe("key");
    // Both members carry the phase — the connector is rooted in their bands.
    expect(verts.every((v) => v.phase === "Pn3m")).toBe(true);
    // Each member's vertex y is the baseline (bottom) of its band.
    const yByPos = verts
      .filter((v) => v.key === "Pn3m-0")
      .map((v) => v.y);
    expect(yByPos).toEqual([50, 100]);
  });

  it("respects peak_ids ordering (Miller-order from the Julia indexer)", () => {
    // peak_ids stored OUT of q-order: id 12 at q=0.10 (first reflection),
    // id 11 at q=0.20 (second). The layer must connect by Miller-order
    // *index* into peak_ids, not by q-sort.
    const m1 = member(1, "Im3m", [12, 11], [0.10, 0.20]);
    const m2 = member(2, "Im3m", [22, 21], [0.11, 0.22]);
    // Reorder effective_peaks so the resolver has to look up by id, not by
    // peak position. (`member` already builds peaks in peakIds order; this
    // exercise just confirms the look-up.)
    const marks = buildCrossTraceTrackingMarks({
      members: [m1, m2],
      yBands: [[0, 50], [50, 100]],
    });
    expect(marks).toHaveLength(1);
    const verts = (Plot.line as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![0] as Array<{ key: string; q: number; order: number }>;
    const k0 = verts.filter((v) => v.order === 0).map((v) => v.q);
    const k1 = verts.filter((v) => v.order === 1).map((v) => v.q);
    expect(k0).toEqual([0.10, 0.11]); // first Miller order in each member
    expect(k1).toEqual([0.20, 0.22]); // second Miller order in each member
  });

  it("colours each polyline by phase and ramps opacity down on higher Miller orders", () => {
    const m1 = member(1, "Pn3m", [11, 12, 13], [0.10, 0.12, 0.14]);
    const m2 = member(2, "Pn3m", [21, 22, 23], [0.11, 0.13, 0.15]);
    buildCrossTraceTrackingMarks({
      members: [m1, m2],
      yBands: [[0, 50], [50, 100]],
    });
    const verts = (Plot.line as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![0] as Array<{ stroke: string; opacity: number; order: number }>;
    expect(verts.every((v) => v.stroke === phaseColor("Pn3m"))).toBe(true);
    // Opacity ramps down: order 0 > order 1 > order 2.
    const op = (k: number) => verts.find((v) => v.order === k)!.opacity;
    expect(op(0)).toBeGreaterThan(op(1));
    expect(op(1)).toBeGreaterThan(op(2));
  });

  it("connects two phases independently when both have ≥ 2 carriers", () => {
    const m1 = member(1, "Pn3m", [11, 12], [0.10, 0.20]);
    const m2 = member(2, "Pn3m", [21, 22], [0.11, 0.22]);
    const m3 = member(3, "Hexagonal", [31, 32], [0.05, 0.10]);
    const m4 = member(4, "Hexagonal", [41, 42], [0.06, 0.12]);
    buildCrossTraceTrackingMarks({
      members: [m1, m2, m3, m4],
      yBands: [[0, 25], [25, 50], [50, 75], [75, 100]],
    });
    const verts = (Plot.line as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![0] as Array<{ phase: string; key: string }>;
    const phases = new Set(verts.map((v) => v.phase));
    expect(phases).toEqual(new Set(["Pn3m", "Hexagonal"]));
    // Per phase, two Miller orders (k=0, k=1) × two carriers = 4 verts each.
    expect(verts.filter((v) => v.phase === "Pn3m")).toHaveLength(4);
    expect(verts.filter((v) => v.phase === "Hexagonal")).toHaveLength(4);
  });

  it("does NOT drop a phase when a third member is unindexed in between (only carriers connect)", () => {
    // 1: Pn3m, 2: <none>, 3: Pn3m — the polyline jumps over m2's row, which
    // is the visual intent: the tracking line is between carrying members.
    const m1 = member(1, "Pn3m", [11, 12], [0.10, 0.20]);
    const m2 = member(2, null, [], []);
    const m3 = member(3, "Pn3m", [31, 32], [0.12, 0.24]);
    buildCrossTraceTrackingMarks({
      members: [m1, m2, m3],
      yBands: [[0, 33], [33, 66], [66, 100]],
    });
    const verts = (Plot.line as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![0] as Array<{ q: number; y: number; order: number }>;
    // Two carriers × two Miller orders = 4 verts. The y values are the
    // baselines of m1 (33) and m3 (100). m2's baseline (66) is absent.
    const k0Ys = verts.filter((v) => v.order === 0).map((v) => v.y).sort((a, b) => a - b);
    expect(k0Ys).toEqual([33, 100]);
  });
});
