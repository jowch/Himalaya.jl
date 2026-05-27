import { describe, it, expect } from "vitest";
import { qToRadius, nearestRingQ, RING_VIEWBOX } from "../src/lib/qRing";

describe("qToRadius (normalized, mockup parity)", () => {
  it("maps qLo to the beamstop offset and qHi to the outer radius", () => {
    const lo = qToRadius(0.04, 0.04, 0.30);
    const hi = qToRadius(0.30, 0.04, 0.30);
    expect(lo).toBeCloseTo(12, 5);        // RING_R_MIN
    expect(hi).toBeCloseTo(12 + 33, 5);   // RING_R_MIN + RING_R_SPAN
    expect(hi).toBeLessThan(RING_VIEWBOX / 2); // stays inside the 100x100 box
  });

  it("is monotonic in q", () => {
    const a = qToRadius(0.1, 0.04, 0.30);
    const b = qToRadius(0.2, 0.04, 0.30);
    expect(b).toBeGreaterThan(a);
  });

  it("degenerate qLo==qHi returns the min radius (no divide-by-zero)", () => {
    expect(qToRadius(0.1, 0.1, 0.1)).toBeCloseTo(12, 5);
  });
});

describe("nearestRingQ", () => {
  const peaks = [0.045, 0.103, 0.206];
  it("returns the nearest peak q within tolerance", () => {
    expect(nearestRingQ(0.104, peaks, 0.005)).toBe(0.103);
  });
  it("returns undefined when nothing is within tolerance", () => {
    expect(nearestRingQ(0.5, peaks, 0.005)).toBeUndefined();
  });
  it("returns undefined for an undefined hover", () => {
    expect(nearestRingQ(undefined, peaks, 0.005)).toBeUndefined();
  });
});
