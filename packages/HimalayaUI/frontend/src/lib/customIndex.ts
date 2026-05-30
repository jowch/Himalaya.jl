// customIndex.ts (Plan D-9) — pure physics for the custom-index modal.
//
// Symmetry → predicted reflections from a user-chosen lattice, the basis the
// backend stores (the q₁ slope), and snap-to-peaks helpers. Mirrors the mockup
// `SYMS` / `customRefls` / `snapValues` (docs/redesign-mockups/2026-05-29-focus-plot.html).

const TWO_PI = 2 * Math.PI;

export type SymmetryKind = "cubic" | "lamellar" | "hex";

export interface SymmetrySpec {
  kind: SymmetryKind;
  /** Allowed reflection indices N (√N orders for cubic; n for lamellar; M for hex). */
  Ms: number[];
  /** The lattice parameter label (a or d). */
  param: "a" | "d";
  def: number;
  min: number;
  max: number;
}

export const SYMS: Record<string, SymmetrySpec> = {
  Pn3m:      { kind: "cubic",    Ms: [2, 3, 4, 6, 8, 9],      param: "a", def: 197, min: 120, max: 320 },
  Im3m:      { kind: "cubic",    Ms: [2, 4, 6, 8, 10, 12],    param: "a", def: 252, min: 120, max: 360 },
  Ia3d:      { kind: "cubic",    Ms: [6, 8, 14, 16, 20, 22],  param: "a", def: 218, min: 120, max: 360 },
  Lamellar:  { kind: "lamellar", Ms: [1, 2, 3, 4, 5],         param: "d", def: 60,  min: 30,  max: 130 },
  Hexagonal: { kind: "hex",      Ms: [1, 3, 4, 7, 9, 12],     param: "a", def: 70,  min: 40,  max: 160 },
};

export interface CustomReflection { N: number; q: number }

/** Predicted reflections for a symmetry at lattice value `val`. */
export function customRefls(sym: string, val: number): CustomReflection[] {
  const s = SYMS[sym];
  if (!s || val <= 0) return [];
  if (s.kind === "lamellar") return s.Ms.map((n) => ({ N: n, q: (TWO_PI * n) / val }));
  if (s.kind === "hex") {
    const q1 = (4 * Math.PI) / (Math.sqrt(3) * val);
    return s.Ms.map((M) => ({ N: M, q: q1 * Math.sqrt(M) }));
  }
  // cubic: q = 2π√N / a
  return s.Ms.map((N) => ({ N, q: (TWO_PI * Math.sqrt(N)) / val }));
}

/**
 * The basis the backend stores — the FIRST reflection's q (the q₁ slope).
 * This is the value such that `predicted_q_for_phase(phase, basis)` (which
 * multiplies by the NORMALIZED ratios, first entry 1.0) reproduces this comb.
 * Equivalent to `2π/a × first(phaseratios(P))` (review finding #4).
 */
export function basisFor(sym: string, val: number): number {
  const refls = customRefls(sym, val);
  return refls.length ? refls[0]!.q : 0;
}

/** How many predicted orders land on an observed peak (within relTol). */
export function landsOn(sym: string, val: number, peakQs: number[], relTol = 0.022): number {
  let lands = 0;
  for (const rf of customRefls(sym, val)) {
    if (peakQs.some((pq) => Math.abs(pq - rf.q) / rf.q < relTol)) lands++;
  }
  return lands;
}

/**
 * Lattice values at which ANY allowed order lands exactly on an observed peak
 * (the magnetic snap targets). For each (peak, allowed order) solve for the
 * lattice that places that order on that peak.
 */
export function snapValues(sym: string, peakQs: number[]): number[] {
  const s = SYMS[sym];
  if (!s) return [];
  const out: number[] = [];
  for (const pq of peakQs) {
    if (pq <= 0) continue;
    for (const M of s.Ms) {
      const v = s.kind === "lamellar"
        ? (TWO_PI * M) / pq
        : s.kind === "hex"
          ? (4 * Math.PI * Math.sqrt(M)) / (Math.sqrt(3) * pq)
          : (TWO_PI * Math.sqrt(M)) / pq;
      if (v >= s.min && v <= s.max) out.push(v);
    }
  }
  return out;
}

/** Magnetic-drag snap: pull a raw lattice value to the nearest snap target. */
export function snapTo(sym: string, raw: number, peakQs: number[], window = 0.75): number {
  const targets = snapValues(sym, peakQs);
  let best = raw; let bestD = window;
  for (const v of targets) {
    const d = Math.abs(v - raw);
    if (d < bestD) { bestD = d; best = v; }
  }
  return best;
}

/** The lattice that puts the first allowed order on a clicked peak q. */
export function latticeForFirstOrderOnPeak(sym: string, pq: number): number {
  const s = SYMS[sym];
  if (!s || pq <= 0) return 0;
  const M0 = s.Ms[0]!;
  return s.kind === "lamellar"
    ? (TWO_PI * M0) / pq
    : s.kind === "hex"
      ? (4 * Math.PI * Math.sqrt(M0)) / (Math.sqrt(3) * pq)
      : (TWO_PI * Math.sqrt(M0)) / pq;
}
