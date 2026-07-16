import { describe, it, expect } from "vitest";
import {
  customRefls, basisFor, landsOn, snapValues, snapTo, latticeForFirstOrderOnPeak,
  latticeBounds, SYMS,
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

  // ── F13: the three phases the custom-index modal was missing. The allowed
  // reflection indices (Ms) are the squared magnitudes h²+k²+l² behind core
  // Himalaya's `phaseratios` (src/phase.jl) — the single source of truth — so a
  // user-fit comb in the UI matches the indexer's own ratio series exactly.
  it("Fm3m refls follow the √3,√4,√8,√11,√12 series (cubic 2π√N/a)", () => {
    const a = 200;
    const r = customRefls("Fm3m", a);
    expect(r.map((x) => x.N)).toEqual([3, 4, 8, 11, 12]);
    expect(r[0]!.q).toBeCloseTo((2 * Math.PI * Math.sqrt(3)) / a, 6);
  });

  it("Fd3m refls follow the √3,√8,√11,√12,√16,√19,√24,√27,√32,√35,√36 series", () => {
    const a = 250;
    const r = customRefls("Fd3m", a);
    expect(r.map((x) => x.N)).toEqual([3, 8, 11, 12, 16, 19, 24, 27, 32, 35, 36]);
    expect(r[0]!.q).toBeCloseTo((2 * Math.PI * Math.sqrt(3)) / a, 6);
  });

  it("Square (2D) refls follow √(h²+k²): 1,2,4,5,8,9,10,13,16,17,18,20 at 2π√N/a", () => {
    const a = 80;
    const r = customRefls("Square", a);
    expect(r.map((x) => x.N)).toEqual([1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20]);
    // first order q = 2π/a (N = 1)
    expect(r[0]!.q).toBeCloseTo((2 * Math.PI) / a, 6);
    // √2 reflection
    expect(r[1]!.q).toBeCloseTo((2 * Math.PI * Math.sqrt(2)) / a, 6);
  });

  it("SYMS carries all eight canonical phases", () => {
    expect(Object.keys(SYMS).sort()).toEqual(
      ["Fd3m", "Fm3m", "Hexagonal", "Ia3d", "Im3m", "Lamellar", "Pn3m", "Square"],
    );
  });

  it("every phase's baseline envelope spans at least 1 nm … 50 nm (deliberately wide)", () => {
    // Regression floor for the user's requirement: 10 Å (1 nm) reachable for all
    // phases (Square down to 1 nm was the driving case) and ≥ 500 Å (50 nm) up top.
    for (const [name, s] of Object.entries(SYMS)) {
      expect(s.min, `${name} floor`).toBeLessThanOrEqual(10);
      expect(s.max, `${name} ceiling`).toBeGreaterThanOrEqual(500);
    }
  });
});

// ── latticeBounds: widen-only, q-window-aware slider bounds ────────────────────
// The slider is the static envelope (SYMS.min/max) UNIONED with what the trace's
// own q-window can reach — it errs wide: never narrower than the envelope, and
// wider when the trace reaches beyond it. a = 2π√(Ms[0])/q (cubic/square/lamellar
// first order), a = 4π/(√3·q) (hex). q_max → smallest a, q_min → largest a.
describe("latticeBounds — widen-only q-window bounds", () => {
  it("returns the bare envelope when no q-window is available (trace not loaded)", () => {
    expect(latticeBounds("Pn3m", null)).toEqual({ min: SYMS.Pn3m!.min, max: SYMS.Pn3m!.max });
  });

  it("is never narrower than the envelope, whatever the q-window", () => {
    for (const win of [[0.1, 0.5], [0.02, 0.9], [0.8, 2.7], [0.002, 0.3]] as [number, number][]) {
      for (const name of Object.keys(SYMS)) {
        const b = latticeBounds(name, win);
        expect(b.min).toBeLessThanOrEqual(SYMS[name]!.min);
        expect(b.max).toBeGreaterThanOrEqual(SYMS[name]!.max);
      }
    }
  });

  it("widens BELOW the floor when a trace reaches high q (sub-1 nm becomes selectable)", () => {
    // Square, WAXS-wide q_max 2.7 → a = 2π/2.7 ≈ 2.33 Å < the 10 Å envelope floor.
    const b = latticeBounds("Square", [0.05, 2.7]);
    expect(b.min).toBeCloseTo((2 * Math.PI) / 2.7, 4);
    expect(b.min).toBeLessThan(SYMS.Square!.min);
  });

  it("widens ABOVE the ceiling when a trace reaches tiny q (super-swollen becomes selectable)", () => {
    // Lamellar, USAXS q_min 0.002 → d = 2π/0.002 ≈ 3142 Å > the 1000 Å envelope ceiling.
    const b = latticeBounds("Lamellar", [0.002, 0.3]);
    expect(b.max).toBeCloseTo((2 * Math.PI) / 0.002, 3);
    expect(b.max).toBeGreaterThan(SYMS.Lamellar!.max);
  });

  it("uses the hexagonal law a = 4π/(√3·q) when widening a hex phase", () => {
    // q_max 2.5 → a = 4π/(√3·2.5) ≈ 2.90 Å, below the hex floor → widened.
    const b = latticeBounds("Hexagonal", [0.05, 2.5]);
    expect(b.min).toBeCloseTo((4 * Math.PI) / (Math.sqrt(3) * 2.5), 4);
  });

  it("a q-window well inside the envelope leaves the slider at the envelope", () => {
    // Square [0.1, 0.5] → derived a ∈ [12.6, 62.8], both inside [10, 500] → no widening.
    expect(latticeBounds("Square", [0.1, 0.5])).toEqual({ min: SYMS.Square!.min, max: SYMS.Square!.max });
  });
});
