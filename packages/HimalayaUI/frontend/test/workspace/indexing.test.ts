import { describe, it, expect } from "vitest";
import { claimPeaks, clampLattice, draftCommit, indexFit, initialDraft, reviseDraft, storedFit } from "../../src/workspace/indexing";
import type { IndexPeakRef } from "../../src/api";

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

  it("gives a contested peak to the closer order, as the backend claims it", () => {
    expect(claimPeaks([0.1, 0.103], [0.1025])).toEqual([null, 0.1025]);
    expect(claimPeaks([0.1, 0.103], [0.1005])).toEqual([0.1005, null]);
  });

  it("an order whose best peak went to a closer order stays unmatched", () => {
    expect(claimPeaks([0.1, 0.103], [0.1028, 0.0985])).toEqual([0.0985, 0.1028]);
    expect(claimPeaks([0.1, 0.1025], [0.102, 0.0979])).toEqual([null, 0.102]);
  });

  it("pins the slope and R² when the first peak sits off the comb", () => {
    const fit = indexFit([0.1, 0.2, 0.3], [0.102, 0.2, 0.3]);
    const slope = (0.102 * 1 + 0.2 * 2 + 0.3 * 3) / (1 + 4 + 9);
    expect(fit.slope).toBeCloseTo(slope, 12);
    const obs = [0.102, 0.2, 0.3];
    const mean = obs.reduce((a, b) => a + b, 0) / 3;
    const ssr = obs.reduce((a, o, i) => a + (o - slope * (i + 1)) ** 2, 0);
    const sst = obs.reduce((a, o) => a + (o - mean) ** 2, 0);
    expect(fit.r2).toBeCloseTo(1 - ssr / sst, 12);
    expect(fit.r2!).toBeLessThan(0.9999);
  });

  it("counts only the orders inside the trace's q range", () => {
    const fit = indexFit([0.1, 0.2, 0.3, 0.4], [0.1, 0.2], [0.05, 0.25]);
    expect(fit.points.map((p) => p.predicted)).toEqual([0.1, 0.2]);
    expect(fit.landed).toBe(2);
    expect(indexFit([0.1, 0.2, 0.3, 0.4], [0.1, 0.2]).points).toHaveLength(4);
  });

  it("plots a stored index from the peaks the backend claimed, not a fresh match", () => {
    const ref = (ratio_position: unknown, q_observed: number): IndexPeakRef =>
      ({ peak_id: 1, ratio_position, residual: 0, q_observed }) as IndexPeakRef;
    const ix = { basis: 0.1, predicted_q: [0.1, 0.2, 0.3, 0.4], peaks: [ref(3, 0.31), ref("\u0002", 0.205)] };
    const fit = storedFit(ix, [0.05, 0.35]);
    expect(fit.points.map((p) => p.observed)).toEqual([null, 0.205, 0.31]);
    expect(fit.landed).toBe(2);
  });

  it("commits the basis and the ratios the comb drew", () => {
    const c = draftCommit({ phase: "Lamellar", lattice: 62.8 });
    expect(c.basis).toBeCloseTo((2 * Math.PI) / 62.8, 9);
    c.ratios.forEach((r, i) => expect(r).toBeCloseTo(i + 1, 9));
    expect(c.ratios).toHaveLength(11);
  });

  it("clamps the lattice to the phase's range and rounds to 0.1 Å", () => {
    expect(clampLattice("Pn3m", 5, null)).toBe(10);
    expect(clampLattice("Pn3m", 123.456, null)).toBe(123.5);
  });

  it("starts from the top candidate, or Pn3m when there is none", () => {
    expect(initialDraft({ phase: "Im3m", lattice_d: 140.04 })).toEqual({ phase: "Im3m", lattice: 140 });
    expect(initialDraft(undefined).phase).toBe("Pn3m");
  });

  it("lets go of a snapped peak when the lattice moves by hand or the peak is gone, and keeps it when only the symmetry changes", () => {
    const snap = { id: 3, source: "manual" as const };
    const snapped = { phase: "Pn3m", lattice: (2 * Math.PI * Math.SQRT2) / 0.12, snap };
    const at = (q: number | undefined) => () => q;
    expect(reviseDraft(snapped, { ...snapped, lattice: snapped.lattice + 0.5 }, at(0.12))).toEqual({ phase: "Pn3m", lattice: snapped.lattice + 0.5 });
    expect(reviseDraft(snapped, { ...snapped, phase: "Lamellar" }, at(undefined))).toEqual({ phase: "Lamellar", lattice: snapped.lattice });
    const moved = reviseDraft(snapped, { ...snapped, phase: "Lamellar" }, at(0.12));
    expect(moved.snap).toEqual(snap);
    expect(moved.lattice).toBeCloseTo((2 * Math.PI) / 0.12, 9);
    expect(reviseDraft(snapped, { ...snapped }, at(0.2))).toEqual(snapped);
  });
});
