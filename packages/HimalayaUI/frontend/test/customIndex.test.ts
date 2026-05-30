import { describe, it, expect } from "vitest";
import {
  customRefls, basisFor, landsOn, snapValues, snapTo, latticeForFirstOrderOnPeak, SYMS,
} from "../src/lib/customIndex";

describe("customIndex physics (Plan D-9)", () => {
  it("cubic refls follow 2π√N/a", () => {
    const r = customRefls("Pn3m", 197);
    expect(r[0]!.N).toBe(2);
    expect(r[0]!.q).toBeCloseTo((2 * Math.PI * Math.sqrt(2)) / 197, 6);
  });

  it("lamellar refls follow 2πn/d", () => {
    const r = customRefls("Lamellar", 60);
    expect(r[0]!.q).toBeCloseTo((2 * Math.PI * 1) / 60, 6);
  });

  it("hex refls follow q1 = 4π/(√3·a)", () => {
    const r = customRefls("Hexagonal", 70);
    const q1 = (4 * Math.PI) / (Math.sqrt(3) * 70);
    expect(r.find((x) => x.N === 1)!.q).toBeCloseTo(q1, 6);
  });

  it("basisFor returns the FIRST reflection q (the q₁ slope) — Ia3d √6", () => {
    const a = 100;
    expect(basisFor("Ia3d", a)).toBeCloseTo((2 * Math.PI * Math.sqrt(6)) / a, 8);
  });

  it("landsOn counts orders that hit observed peaks", () => {
    const a = 197;
    const peakQs = customRefls("Pn3m", a).slice(0, 2).map((r) => r.q);
    expect(landsOn("Pn3m", a, peakQs)).toBe(2);
  });

  it("snapTo pulls 200→198 / 196→197 but leaves 215 free (mockup parity)", () => {
    // A Pn3m peak that 'wants' a=197 (first order at q = 2π√2/197).
    const peakQ = (2 * Math.PI * Math.sqrt(2)) / 197;
    // 197 is the snap target; a window of 0.75 absorbs nearby drags.
    expect(snapTo("Pn3m", 197.3, [peakQ])).toBeCloseTo(197, 3);
    expect(snapTo("Pn3m", 196.7, [peakQ])).toBeCloseTo(197, 3);
    // 215 is well outside the window → stays free.
    expect(snapTo("Pn3m", 215, [peakQ])).toBe(215);
  });

  it("snapValues stay within the symmetry's [min,max] band", () => {
    const vals = snapValues("Pn3m", [0.045, 0.055]);
    for (const v of vals) {
      expect(v).toBeGreaterThanOrEqual(SYMS.Pn3m!.min);
      expect(v).toBeLessThanOrEqual(SYMS.Pn3m!.max);
    }
  });

  it("latticeForFirstOrderOnPeak puts the first order on a clicked peak", () => {
    const pq = 0.045;
    const a = latticeForFirstOrderOnPeak("Pn3m", pq);
    // first order back-computes to pq
    expect(customRefls("Pn3m", a)[0]!.q).toBeCloseTo(pq, 8);
  });
});
