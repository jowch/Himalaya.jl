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

  it("clips overflow at the bottom of the total band envelope", () => {
    // Negative-intensity points (rare but possible after baseline subtraction)
    // would map BELOW the working band; should clip at total envelope bottom.
    const t = { q: [0, 1], I: [-50, 1] };
    const out = applyNormalization(t, 1, [0, 100], 0.7);
    // The negative point's mapped y exceeds the total band bottom (100) → clip.
    expect(out[0]!.y).toBe(100);
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
});
