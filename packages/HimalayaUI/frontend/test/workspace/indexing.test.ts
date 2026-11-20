import { describe, it, expect } from "vitest";
import { clampLattice, draftCommit, indexFit, initialDraft } from "../../src/workspace/indexing";

describe("indexing", () => {
  it("fits q against the reflection ratio through the origin", () => {
    const predicted = [0.1, 0.1 * Math.sqrt(1.5), 0.1 * Math.sqrt(2), 0.1 * Math.sqrt(3)];
    const peaks = [0.1, 0.1 * Math.sqrt(1.5) * 1.01, 0.1 * Math.sqrt(3) * 0.995, 0.05];
    const fit = indexFit(predicted, peaks);
    expect(fit.landed).toBe(3);
    expect(fit.points.map((p) => p.observed !== null)).toEqual([true, true, false, true]);
    expect(fit.slope).toBeCloseTo(0.1, 2);
    expect(fit.r2).toBeGreaterThan(0.99);
  });

  it("uses each peak once and reports no R² below two matches", () => {
    const fit = indexFit([0.1, 0.1005], [0.1]);
    expect(fit.landed).toBe(1);
    expect(fit.r2).toBeNull();
  });

  it("commits the basis and the ratios the comb drew", () => {
    const c = draftCommit({ phase: "Lamellar", lattice: 62.8 });
    expect(c.basis).toBeCloseTo((2 * Math.PI) / 62.8, 9);
    c.ratios.forEach((r, i) => expect(r).toBeCloseTo(i + 1, 9));
    expect(c.ratios).toHaveLength(5);
  });

  it("clamps the lattice to the phase's range and rounds to 0.1 Å", () => {
    expect(clampLattice("Pn3m", 5, null)).toBe(10);
    expect(clampLattice("Pn3m", 123.456, null)).toBe(123.5);
  });

  it("starts from the top candidate, or Pn3m when there is none", () => {
    expect(initialDraft({ phase: "Im3m", lattice_d: 140.04 })).toEqual({ phase: "Im3m", lattice: 140 });
    expect(initialDraft(undefined).phase).toBe("Pn3m");
  });
});
