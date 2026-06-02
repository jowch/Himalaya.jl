import { describe, it, expect } from "vitest";
import {
  buildPrintDetectorLut, PRINT_DETECTOR_STOPS,
} from "../../src/print/detector/detectorLut";

const luma = (lut: Uint8ClampedArray, i: number) =>
  0.2126 * lut[i*3] + 0.7152 * lut[i*3+1] + 0.0722 * lut[i*3+2];

describe("buildPrintDetectorLut (default = neutral)", () => {
  const lut = buildPrintDetectorLut();

  it("has 256 RGB triples (768 bytes)", () => {
    expect(lut.length).toBe(256 * 3);
  });

  it("is non-inverting — luminance is monotonically non-decreasing", () => {
    for (let i = 1; i < 256; i++) {
      // Allow a tiny epsilon for rounding; the ramp must never step down.
      expect(luma(lut, i)).toBeGreaterThanOrEqual(luma(lut, i - 1) - 1.5);
    }
  });

  it("endpoints: index 0 is dark, index 255 is light", () => {
    const sum0 = lut[0] + lut[1] + lut[2];
    const sum255 = lut[255*3] + lut[255*3+1] + lut[255*3+2];
    expect(sum0).toBeLessThan(180);          // near-black window backing
    expect(sum255).toBeGreaterThan(600);     // paper-cream highlight
  });

  it("is NEAR-NEUTRAL — barely any chroma, so the phase rings own all the colour", () => {
    // R, G, B stay within a tight band at every index (a warm gray, not a hue).
    for (let i = 0; i < 256; i++) {
      const r = lut[i*3], g = lut[i*3+1], b = lut[i*3+2];
      const spread = Math.max(r, g, b) - Math.min(r, g, b);
      expect(spread).toBeLessThan(18);       // saturated terracotta would be ~50+
    }
  });

  it("is faintly WARM (R ≥ B), keeping the Print paper feel", () => {
    // Sampled across the ramp: red channel never below blue (warm-gray, not cool).
    for (const t of [0.1, 0.4, 0.7, 1.0]) {
      const i = Math.round(t * 255);
      expect(lut[i*3]).toBeGreaterThanOrEqual(lut[i*3+2]);
    }
  });

  it("exposes its stops for the Storybook swatch", () => {
    expect(PRINT_DETECTOR_STOPS.length).toBeGreaterThanOrEqual(4);
    expect(PRINT_DETECTOR_STOPS[0].t).toBe(0);
    expect(PRINT_DETECTOR_STOPS[PRINT_DETECTOR_STOPS.length - 1].t).toBe(1);
  });
});

describe("buildPrintDetectorLut('warm') — the optional beauty ramp", () => {
  const lut = buildPrintDetectorLut("warm");

  it("is non-inverting and ends light + warm", () => {
    for (let i = 1; i < 256; i++) {
      expect(luma(lut, i)).toBeGreaterThanOrEqual(luma(lut, i - 1) - 1.5);
    }
    expect(lut[255*3] + lut[255*3+1] + lut[255*3+2]).toBeGreaterThan(600);
    expect(lut[255*3]).toBeGreaterThanOrEqual(lut[255*3+2]); // warm: R >= B
  });

  it("blooms warm in the mid-highs (terracotta: R clearly exceeds B near t~0.62)", () => {
    const i = Math.round(0.62 * 255);
    expect(lut[i*3] - lut[i*3+2]).toBeGreaterThan(20);
  });

  it("is a DISTINCT, more saturated ramp than the neutral default", () => {
    const neutral = buildPrintDetectorLut();
    const i = Math.round(0.62 * 255);
    const warmSpread = Math.max(lut[i*3], lut[i*3+1], lut[i*3+2]) - Math.min(lut[i*3], lut[i*3+1], lut[i*3+2]);
    const neutralSpread = Math.max(neutral[i*3], neutral[i*3+1], neutral[i*3+2]) - Math.min(neutral[i*3], neutral[i*3+1], neutral[i*3+2]);
    expect(warmSpread).toBeGreaterThan(neutralSpread + 20);
  });
});
