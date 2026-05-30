/**
 * series/anchors tests (Plan E, Task E-2).
 *
 * `buildAnchorMap(members)` keys reflections by `"phase:order"` and yields one
 * vertex per carrying member. A member that carries the phase but does NOT
 * observe order k yields a `q = null` (predicted-absent) vertex at the
 * predicted q (derived from `lattice_d` × phase ratios — the same physics as
 * Focus combs). A member without the phase produces no vertex (the track ends).
 */
import { describe, it, expect } from "vitest";
import type { SeriesMember } from "../src/api";
import { buildAnchorMap, predictedQForOrder } from "../src/lib/series/anchors";

// Build a member whose confirmed_index claims `phase` with the given
// `lattice_d`, observing exactly the peaks at `observedQs` (in Miller order).
function member(
  id: number,
  phase: string | null,
  lattice_d: number,
  observedQs: number[],
): SeriesMember {
  const peakIds = observedQs.map((_, i) => id * 100 + i);
  const effective_peaks = observedQs.map((q, i) => ({
    id: peakIds[i]!, q, intensity: 50, sharpness: 1, source: "auto" as const,
  }));
  return {
    id, series_id: 1, exposure_id: id * 10, display_order: id - 1,
    band_height: 1, y_offset: 0, normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: phase
      ? {
          effective_peaks,
          confirmed_index: {
            id: id * 1000, phase, lattice_d, r_squared: 0.99, ngc: -1.5,
            peak_ids: peakIds,
          },
          analysis_inputs_hash: "h",
        }
      : { effective_peaks, confirmed_index: null, analysis_inputs_hash: "h" },
    is_stale: false, created_by: null, created_at: null,
  };
}

describe("predictedQForOrder", () => {
  it("derives predicted q from lattice_d × phase ratio (cubic Pn3m √2)", () => {
    // q(cubic, order k) = 2π·√(radicand[k]) / a. Pn3m radicands = [2,3,4,…].
    const a = 205;
    expect(predictedQForOrder("Pn3m", a, 0)).toBeCloseTo((2 * Math.PI * Math.sqrt(2)) / a, 6);
    expect(predictedQForOrder("Pn3m", a, 1)).toBeCloseTo((2 * Math.PI * Math.sqrt(3)) / a, 6);
  });
  it("derives lamellar orders as 2π·n / d", () => {
    const d = 60;
    expect(predictedQForOrder("Lamellar", d, 0)).toBeCloseTo((2 * Math.PI * 1) / d, 6);
    expect(predictedQForOrder("Lamellar", d, 2)).toBeCloseTo((2 * Math.PI * 3) / d, 6);
  });
  it("returns null for an unknown phase", () => {
    expect(predictedQForOrder("Bogus", 100, 0)).toBeNull();
  });
});

describe("buildAnchorMap", () => {
  it("keys by phase:order and records one vertex per carrying member", () => {
    const a = 205;
    const q0 = (2 * Math.PI * Math.sqrt(2)) / a;
    const q1 = (2 * Math.PI * Math.sqrt(3)) / a;
    const m1 = member(1, "Pn3m", a, [q0, q1]);
    const m2 = member(2, "Pn3m", a, [q0, q1]);
    const map = buildAnchorMap([m1, m2]);
    const o0 = map.get("Pn3m:0")!;
    expect(o0).toHaveLength(2);
    expect(o0.map((v) => v.memberPos)).toEqual([0, 1]);
    expect(o0.every((v) => v.q !== null)).toBe(true);
  });

  it("yields a q=null absent vertex when the phase is present but order k is unobserved", () => {
    // The mockup's 1:0.25 case: Pn3m present, but the √3 (order 1) reflection
    // is predicted-but-not-observed in this member.
    const a = 198;
    const q0 = (2 * Math.PI * Math.sqrt(2)) / a;  // observed
    const q2 = (2 * Math.PI * Math.sqrt(4)) / a;  // observed
    // order 1 (√3) is NOT observed
    const present = member(1, "Pn3m", a, [q0, q2]);
    const full = member(2, "Pn3m", a, [
      (2 * Math.PI * Math.sqrt(2)) / a,
      (2 * Math.PI * Math.sqrt(3)) / a,
      (2 * Math.PI * Math.sqrt(4)) / a,
    ]);
    const map = buildAnchorMap([present, full]);
    const o1 = map.get("Pn3m:1")!;
    // Member 0 carries the phase but absents order 1 → q === null with a
    // predictedQ at the √3 position.
    const v0 = o1.find((v) => v.memberPos === 0)!;
    expect(v0.q).toBeNull();
    expect(v0.absent).toBe(true);
    expect(v0.predictedQ).toBeCloseTo((2 * Math.PI * Math.sqrt(3)) / a, 6);
    // Member 1 observes order 1 → q present.
    const v1 = o1.find((v) => v.memberPos === 1)!;
    expect(v1.q).toBeCloseTo((2 * Math.PI * Math.sqrt(3)) / a, 6);
    expect(v1.absent).toBe(false);
  });

  it("produces NO vertex for a member that does not carry the phase (track ends)", () => {
    const a = 205;
    const q0 = (2 * Math.PI * Math.sqrt(2)) / a;
    const carrier = member(1, "Pn3m", a, [q0]);
    const other = member(2, "Lamellar", 60, [(2 * Math.PI) / 60]);
    const map = buildAnchorMap([carrier, other]);
    const o0 = map.get("Pn3m:0")!;
    // Only the carrier appears; the lamellar member contributes no Pn3m vertex.
    expect(o0.map((v) => v.memberPos)).toEqual([0]);
  });

  it("does not emit anchor entries for form-factor / null members", () => {
    const ff = member(1, null, 0, []);
    const map = buildAnchorMap([ff]);
    expect(map.size).toBe(0);
  });
});
