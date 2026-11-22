import { describe, expect, it } from "vitest";
import { formatCurvature, gaussianCurvature } from "../../src/workspace/curvature";

describe("gaussianCurvature", () => {
  it("gives the mean Gaussian curvature in nm⁻² from a lattice in Å", () => {
    expect(gaussianCurvature("Pn3m", 100)).toBeCloseTo((2 * Math.PI * -2) / (1.919 * 100), 10);
    expect(gaussianCurvature("Im3m", 150)).toBeCloseTo((2 * Math.PI * -4) / (2.345 * 225), 10);
    expect(gaussianCurvature("Ia3d", 200)).toBeCloseTo((2 * Math.PI * -8) / (3.091 * 400), 10);
  });

  it("is always negative for the bicontinuous cubics", () => {
    for (const phase of ["Pn3m", "Im3m", "Ia3d"]) expect(gaussianCurvature(phase, 120)).toBeLessThan(0);
  });

  it("agrees across Bonnet-related cubics, whose lattices keep equal curvature", () => {
    const d = gaussianCurvature("Pn3m", 100)!;
    expect(gaussianCurvature("Im3m", 127.9)!).toBeCloseTo(d, 3);
    expect(gaussianCurvature("Ia3d", 157.6)!).toBeCloseTo(d, 3);
  });

  it("is null for phases without a minimal surface and for a missing lattice", () => {
    for (const phase of ["Fm3m", "Fd3m", "Hexagonal", "Lamellar", "Square"]) expect(gaussianCurvature(phase, 100)).toBeNull();
    expect(gaussianCurvature("Pn3m", null)).toBeNull();
    expect(gaussianCurvature("Pn3m", 0)).toBeNull();
  });
});

describe("formatCurvature", () => {
  it("writes κ with a true minus sign and nm⁻²", () => {
    expect(formatCurvature(-0.06548)).toBe("κ −0.0655 nm⁻²");
  });
});
