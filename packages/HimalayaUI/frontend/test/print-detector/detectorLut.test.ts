import { describe, it, expect } from "vitest";
import { buildPrintDetectorLut, PRINT_DETECTOR_STOPS } from "../../src/print/detector/detectorLut";

describe("buildPrintDetectorLut", () => {
  const lut = buildPrintDetectorLut();

  it("has 256 RGB triples (768 bytes)", () => {
    expect(lut.length).toBe(256 * 3);
  });

  it("is non-inverting — luminance is monotonically non-decreasing", () => {
    const luma = (i: number) => 0.2126 * lut[i*3] + 0.7152 * lut[i*3+1] + 0.0722 * lut[i*3+2];
    for (let i = 1; i < 256; i++) {
      // Allow a tiny epsilon for rounding; the ramp must never step down.
      expect(luma(i)).toBeGreaterThanOrEqual(luma(i - 1) - 1.5);
    }
  });

  it("endpoints: index 0 is dark, index 255 is light and warm", () => {
    const sum0 = lut[0] + lut[1] + lut[2];
    const sum255 = lut[255*3] + lut[255*3+1] + lut[255*3+2];
    expect(sum0).toBeLessThan(180);          // near-black window backing
    expect(sum255).toBeGreaterThan(600);     // cream highlight
    expect(lut[255*3]).toBeGreaterThanOrEqual(lut[255*3+2]); // warm: R >= B
  });

  it("blooms warm in the mid-highs (terracotta: R clearly exceeds B near t~0.62)", () => {
    const i = Math.round(0.62 * 255);
    expect(lut[i*3] - lut[i*3+2]).toBeGreaterThan(20);
  });

  it("exposes its stops for the Storybook swatch", () => {
    expect(PRINT_DETECTOR_STOPS.length).toBeGreaterThanOrEqual(4);
    expect(PRINT_DETECTOR_STOPS[0].t).toBe(0);
    expect(PRINT_DETECTOR_STOPS[PRINT_DETECTOR_STOPS.length - 1].t).toBe(1);
  });
});
