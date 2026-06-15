import { describe, it, expect } from "vitest";
import { kappaForPhase } from "../src/lib/curvature";
import { formatKappaPretty } from "../src/lib/units";

describe("kappaForPhase", () => {
  it("matches the backend formula κ = -2π·(χ/A₀)/a² for the three cubics", () => {
    // Im3m: χ=-4, A₀=2.345, a=252 Å → +2π·(4/2.345)/252² ≈ 1.69e-4
    const k = kappaForPhase("Im3m", 252)!;
    expect(k).toBeGreaterThan(0); // χ<0 ⇒ κ reported positive
    expect(k).toBeCloseTo((-2 * Math.PI * (-4 / 2.345)) / (252 * 252), 12);
  });

  it("scales as 1/a² (doubling the lattice quarters κ)", () => {
    const a = kappaForPhase("Pn3m", 100)!;
    const b = kappaForPhase("Pn3m", 200)!;
    expect(b / a).toBeCloseTo(0.25, 6);
  });

  it("returns null for non-bicontinuous-cubic phases", () => {
    expect(kappaForPhase("Fm3m", 200)).toBeNull();
    expect(kappaForPhase("Lamellar", 65)).toBeNull();
    expect(kappaForPhase("Hexagonal", 70)).toBeNull();
  });

  it("returns null for missing / non-positive lattice", () => {
    expect(kappaForPhase("Im3m", null)).toBeNull();
    expect(kappaForPhase("Im3m", 0)).toBeNull();
    expect(kappaForPhase("Im3m", -5)).toBeNull();
  });
});

describe("formatKappaPretty", () => {
  it("renders the small-magnitude branch with a ×10ⁿ superscript", () => {
    expect(formatKappaPretty(1.7e-4)).toBe("1.70×10⁻⁴");
    expect(formatKappaPretty(9.31e-4)).toBe("9.31×10⁻⁴");
  });

  it("passes through fixed-point and zero unchanged", () => {
    expect(formatKappaPretty(0)).toBe("0");
    expect(formatKappaPretty(0.5)).toBe("0.500");
  });
});
