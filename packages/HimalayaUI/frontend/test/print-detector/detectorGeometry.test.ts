import { describe, it, expect } from "vitest";
import {
  qToImageRadius, buildRingPlacements, type DetectorCalibration,
} from "../../src/print/detector/detectorGeometry";

// lambda = 12.398/12.398 = 1.0 Angstrom; 1000 mm distance; 0.1 mm pixels; 1000x1000 image.
const CAL: DetectorCalibration = {
  beamCenterPx: { x: 500, y: 500 },
  imageSizePx: { w: 1000, h: 1000 },
  sampleDistanceMm: 1000,
  pixelSizeMm: 0.1,
  energyKeV: 12.398,
};

describe("qToImageRadius", () => {
  it("matches the hand-computed small-angle radius for q=0.1", () => {
    // theta = asin(0.1*1/(4pi)) = 0.0079577 rad; 2theta -> tan ~ 0.015916;
    // r_mm = 1000*0.015916 = 15.916; r_px = 159.16; /1000 = 0.15916.
    expect(qToImageRadius(0.1, CAL)).toBeCloseTo(0.1592, 3);
  });
  it("grows monotonically with q", () => {
    expect(qToImageRadius(0.2, CAL)).toBeGreaterThan(qToImageRadius(0.1, CAL));
  });
});

describe("buildRingPlacements", () => {
  it("with calibration: beam center is the normalized px center; radii from geometry", () => {
    const { beamCenter, rings } = buildRingPlacements([0.1, 0.2], CAL);
    expect(beamCenter).toEqual({ x: 0.5, y: 0.5 });
    expect(rings.map((r) => r.q)).toEqual([0.1, 0.2]);
    expect(rings[0].r).toBeCloseTo(0.1592, 3);
  });
  it("off-center beam normalizes x directly and FLIPS y (bottom-left cal → top-left screen)", () => {
    const { beamCenter } = buildRingPlacements([0.1], { ...CAL, beamCenterPx: { x: 950, y: 60 } });
    expect(beamCenter.x).toBeCloseTo(0.95, 3);
    // y=60 from the bottom of a 1000px image → 0.94 from the top.
    expect(beamCenter.y).toBeCloseTo(0.94, 3);
  });
  it("null calibration -> presentational fallback: centered, radius prop to q-range", () => {
    const { beamCenter, rings } = buildRingPlacements([0.1, 0.2, 0.3], null);
    expect(beamCenter).toEqual({ x: 0.5, y: 0.5 });
    // innermost ~ RING_R_MIN/VIEWBOX, span up with q; monotonic + within (0, 0.5].
    expect(rings[0].r).toBeGreaterThan(0);
    expect(rings[2].r).toBeGreaterThan(rings[0].r);
    expect(rings[2].r).toBeLessThanOrEqual(0.5);
  });
  it("carries color + ghost through from the q-set entries", () => {
    const { rings } = buildRingPlacements(
      [{ q: 0.1, color: "var(--color-success)" }, { q: 0.2, ghost: true }], CAL,
    );
    expect(rings[0].color).toBe("var(--color-success)");
    expect(rings[1].ghost).toBe(true);
  });
});
