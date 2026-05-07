import { describe, it, expect } from "vitest";
import {
  computeReference,
  applyNormalization,
} from "../src/lib/comparison/normalization";
import type { Peak } from "../src/lib/comparison/normalization";

// ── computeReference ──────────────────────────────────────────────────────────

describe("computeReference", () => {
  const trace = {
    q: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
    I: [10, 20, 30, 50, 80, 60, 40, 25, 15, 12],
  };

  it("'none' returns the global signal max (fills the band without additional scaling)", () => {
    // "none" is no longer a literal passthrough (reference=1) because raw
    // SAXS intensities are in the hundreds–thousands, which would clip every
    // point to the top of the band. It now behaves like "max" with no q-window.
    expect(computeReference(trace, [], null, "none")).toBe(80);
    // With an explicit window the reference is the max inside that window.
    expect(computeReference(trace, [], [0.2, 0.5], "none")).toBe(80);
  });

  it("'max' returns the signal max in the q_window (full trace if window null)", () => {
    expect(computeReference(trace, [], null, "max")).toBe(80);
    // window [0.1, 0.4] only includes I = [10, 20, 30, 50] → max 50.
    expect(computeReference(trace, [], [0.1, 0.4], "max")).toBe(50);
  });

  it("'area' returns the trapezoidal integral of I dq over q_window", () => {
    const flat = { q: [0, 1, 2, 3], I: [1, 1, 1, 1] };
    // trapezoidal of constant 1 from q=0..3 = 3
    expect(computeReference(flat, [], null, "area")).toBeCloseTo(3, 6);
    // restricted to [1, 2] → integral = 1
    expect(computeReference(flat, [], [1, 2], "area")).toBeCloseTo(1, 6);
  });

  it("'qwindow' uses peak-fit when peaks fall in window with non-null intensity", () => {
    const peaks: Peak[] = [
      { q: 0.4, intensity: 50 }, // in window
      { q: 0.5, intensity: 80 }, // in window — max
      { q: 0.9, intensity: 15 }, // out of window
    ];
    expect(computeReference(trace, peaks, [0.3, 0.6], "qwindow")).toBe(80);
  });

  it("'qwindow' falls back to signal max in window when no peaks present", () => {
    expect(computeReference(trace, [], [0.3, 0.7], "qwindow")).toBe(80);
  });

  it("'qwindow' falls back to signal-fit when ALL peaks in window are manual (null intensity)", () => {
    const peaks: Peak[] = [
      { q: 0.4, intensity: null }, // manual, in window
      { q: 0.5, intensity: null }, // manual, in window
    ];
    // Should fall through to signal max in [0.3, 0.6] which is 80.
    expect(computeReference(trace, peaks, [0.3, 0.6], "qwindow")).toBe(80);
  });

  it("'qwindow' uses peak-fit when at least one peak has intensity", () => {
    const peaks: Peak[] = [
      { q: 0.4, intensity: null }, // manual
      { q: 0.5, intensity: 35 },   // auto with intensity → reference = 35
    ];
    expect(computeReference(trace, peaks, [0.3, 0.6], "qwindow")).toBe(35);
  });

  it("does not crash on empty trace", () => {
    expect(computeReference({ q: [], I: [] }, [], null, "max")).toBeGreaterThan(0);
    expect(computeReference({ q: [], I: [] }, [], null, "area")).toBeGreaterThan(0);
    expect(computeReference({ q: [], I: [] }, [], null, "qwindow")).toBeGreaterThan(0);
  });

  it("does not crash with empty q_window (min == max)", () => {
    // empty selection → still safe (positive reference avoids division-by-zero).
    expect(computeReference(trace, [], [0.5, 0.5], "max")).toBeGreaterThan(0);
    expect(computeReference(trace, [], [0.5, 0.5], "area")).toBeGreaterThan(0);
  });
});

// ── applyNormalization ────────────────────────────────────────────────────────

describe("applyNormalization", () => {
  const trace = {
    q: [0.1, 0.2, 0.3, 0.4, 0.5],
    I: [10, 20, 50, 30, 15],
  };

  it("maps reference value to top of working band (default 0.7)", () => {
    // yBand = [0, 100] (top-down y; we treat the band as a numeric envelope).
    // working_band_fraction = 0.7 means inner 70% of yBand.
    // For yBand [0, 100], working band = [15, 85].
    // The TOP of the working band (in plot coordinates where y axis points up)
    // corresponds to the smaller numeric value of yBand if yBand is a vertical
    // interval [top_px, bottom_px]; we let yBand simply describe the numeric
    // range and map intensity → y by linear scale.
    const out = applyNormalization(trace, 50, [0, 100], 0.7);
    expect(out).toHaveLength(5);
    // Reference (50) should map to the top edge of the working band.
    // Working band is centered: top = 15, bottom = 85.
    const refPoint = out.find((p) => p.q === 0.3);
    expect(refPoint).toBeDefined();
    expect(refPoint!.y).toBeCloseTo(15, 5);
  });

  it("maps zero intensity to the bottom of the working band", () => {
    const flat = { q: [0, 1], I: [0, 0] };
    const out = applyNormalization(flat, 100, [0, 100], 0.7);
    // Bottom of working band (in numeric envelope [0,100]) = 85.
    expect(out[0]!.y).toBeCloseTo(85, 5);
    expect(out[1]!.y).toBeCloseTo(85, 5);
  });

  it("clips tail at total band envelope", () => {
    const tall = { q: [0, 1, 2], I: [1, 5, 1] }; // peak (5) far above reference (1)
    const out = applyNormalization(tall, 1, [0, 100], 0.7);
    // y = top - (I/ref)*(workingBandHeight). For I=5 the un-clipped y would be
    // top - 5*70 = 15 - 350 = -335 → clipped to total band top (0).
    const top = out.find((p) => p.q === 1)!;
    expect(top.y).toBe(0);
    // Bottom edge: I=1 maps inside band, no clip.
    expect(out[0]!.y).toBeGreaterThanOrEqual(0);
    expect(out[0]!.y).toBeLessThanOrEqual(100);
  });

  it("clamps negative intensities to the bottom of the working band", () => {
    // Negative-intensity points (rare but possible after baseline subtraction)
    // are undefined in log space; the y-mapping clamps I to 0 before log1p,
    // so a negative point lands at the working band's bottom edge — never
    // strays into the adjacent member's band.
    const t = { q: [0, 1], I: [-50, 1] };
    const out = applyNormalization(t, 1, [0, 100], 0.7);
    // workBottom = bandBottom - padding = 100 - 15 = 85.
    expect(out[0]!.y).toBeCloseTo(85, 5);
  });

  it("respects a custom working_band_fraction", () => {
    // workingBandFraction=1.0 → no headroom; reference maps to top of total band.
    const out = applyNormalization(trace, 50, [0, 100], 1.0);
    const refPoint = out.find((p) => p.q === 0.3)!;
    expect(refPoint.y).toBeCloseTo(0, 5);
  });

  it("defaults workingBandFraction to 0.7 when omitted", () => {
    const out1 = applyNormalization(trace, 50, [0, 100]);
    const out2 = applyNormalization(trace, 50, [0, 100], 0.7);
    for (let i = 0; i < out1.length; i++) {
      expect(out1[i]!.y).toBeCloseTo(out2[i]!.y, 6);
    }
  });

  it("does not crash on a zero reference (avoids division-by-zero)", () => {
    // Treat zero reference as a degenerate trace; should not produce NaN/Infinity.
    const out = applyNormalization(trace, 0, [0, 100], 0.7);
    expect(out).toHaveLength(5);
    for (const p of out) expect(Number.isFinite(p.y)).toBe(true);
  });

  it("does not crash on an empty trace", () => {
    expect(applyNormalization({ q: [], I: [] }, 1, [0, 100])).toEqual([]);
  });

  // Issue #56 regression: a SAXS trace spanning ~3 decades inside one band
  // must not collapse to a smooth roll-off near the floor. Linear y-mapping
  // squashes everything except the highest decade; log mapping keeps the
  // low- and mid-q tail visible. The shape-level assertion: the y-values at
  // ~mid-q (q=0.10) and ~high-q (q=0.20) are well separated relative to the
  // working-band height, even though I drops 100x between those samples.
  it("issue #56: 3-decade SAXS dynamic range stays visible inside the band", () => {
    // Synthesised trace: I drops from 10000 at low q to 10 at high q
    // (3 decades) — typical SAXS shape. Reference = max in window = 10000.
    const trace = {
      q: [0.05, 0.10, 0.15, 0.20, 0.25],
      I:  [10000, 1000, 100, 10, 5],
    };
    const yBand: [number, number] = [0, 100]; // working band = [15, 85], height 70
    const out = applyNormalization(trace, 10000, yBand, 0.7);

    // Reference (q=0.05, I=10000) lands at top of working band (~15).
    expect(out[0]!.y).toBeCloseTo(15, 0);

    // The point at I=10 (1000x below the reference) must be visibly above
    // the working-band floor — i.e. NOT collapsed to the floor (~85). With
    // a linear map, y(I=10) = 85 - (10/10000)*70 ≈ 84.93, indistinguishable
    // from the floor (85) at any pixel resolution. With a log map this is
    // pushed up into the band (around y ≈ 32, well above floor).
    const yHighQ = out[3]!.y;
    expect(yHighQ).toBeLessThan(70);

    // Mid-q (I=1000, 10x below reference) and high-q (I=10, 1000x below
    // reference) must be SEPARATED by a meaningful fraction of the band.
    // Linear mapping gives ∆y ≈ 6.93 px (10% of band). Log mapping gives
    // ∆y ≈ 1/3 of the band (one decade in three).
    const yMidQ = out[1]!.y;
    const spread = Math.abs(yHighQ - yMidQ);
    const workHeight = 70;
    expect(spread).toBeGreaterThan(workHeight * 0.2);
  });
});
