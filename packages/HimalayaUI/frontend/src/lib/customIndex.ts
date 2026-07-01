// customIndex.ts (Plan D-9) — pure physics for the custom-index modal.
//
// Symmetry → predicted reflections from a user-chosen lattice, the basis the
// backend stores (the q₁ slope), and snap-to-peaks helpers. Mirrors the mockup
// `SYMS` / `customRefls` / `snapValues`.

const TWO_PI = 2 * Math.PI;

// "square" is the 2D square lattice. Its q-law is q = 2π√N / a with
// N = h²+k² — algebraically identical to the cubic case (N = h²+k²+l²), so it
// falls through to the cubic branch everywhere kind is switched on. It is named
// distinctly only so the picker reads honestly (a 2D phase, not a cubic one).
export type SymmetryKind = "cubic" | "lamellar" | "hex" | "square";

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

// Allowed reflection indices (Ms) are the squared magnitudes h²+k²+l² (3D) or
// h²+k² (2D) behind core Himalaya's `phaseratios` (src/phase.jl) — the single
// source of truth — so a user-fit comb here reproduces the indexer's own ratio
// series exactly. Fm3m / Fd3m / Square added in F13; their Ms mirror phase.jl
// verbatim.
//
// [min, max] is the BASELINE slider envelope (Å): deliberately wide — a 10 Å
// (1 nm) floor for every phase and generous ceilings covering swollen extremes
// (lamellar super-swells past 100 nm; bicontinuous/micellar cubics reach ~65 nm),
// sized from a SAXS literature review (2026-07). `latticeBounds()` only ever
// WIDENS this per-trace (never narrows), so the slider is never stricter than the
// envelope. `def` is the initial slider value (used by stories); it stays inside
// the envelope.
export const SYMS: Record<string, SymmetrySpec> = {
  Pn3m:      { kind: "cubic",    Ms: [2, 3, 4, 6, 8, 9],      param: "a", def: 197, min: 10, max: 500 },
  Im3m:      { kind: "cubic",    Ms: [2, 4, 6, 8, 10, 12],    param: "a", def: 252, min: 10, max: 700 },
  Ia3d:      { kind: "cubic",    Ms: [6, 8, 14, 16, 20, 22],  param: "a", def: 218, min: 10, max: 650 },
  Fm3m:      { kind: "cubic",    Ms: [3, 4, 8, 11, 12],       param: "a", def: 200, min: 10, max: 500 },
  Fd3m:      { kind: "cubic",    Ms: [3, 8, 11, 12, 16, 19, 24, 27, 32, 35, 36], param: "a", def: 250, min: 10, max: 600 },
  Lamellar:  { kind: "lamellar", Ms: [1, 2, 3, 4, 5],         param: "d", def: 60,  min: 10, max: 1000 },
  Hexagonal: { kind: "hex",      Ms: [1, 3, 4, 7, 9, 12],     param: "a", def: 70,  min: 10, max: 500 },
  Square:    { kind: "square",   Ms: [1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20], param: "a", def: 80, min: 10, max: 500 },
};

/**
 * Slider bounds for the lattice parameter: the static per-phase envelope
 * (SYMS.min/max) WIDENED by whatever the trace's own q-window can reach. We err
 * wide — the range is never stricter than the envelope, and expands further when
 * a trace reaches beyond it (a high q_max opens sub-1 nm lattices; a tiny q_min
 * opens super-swollen ones), so no bound ever needs manual widening per-trace.
 * The first order inverts the same law `customRefls` draws with:
 *   a = 2π√(Ms[0])/q  (cubic / square / lamellar)   a = 4π/(√3·q)  (hex).
 * q_max → smallest a, q_min → largest a. `qWindow === null` (trace not loaded) or
 * a degenerate window → the bare envelope.
 */
export function latticeBounds(
  sym: string,
  qWindow: readonly [number, number] | null,
): { min: number; max: number } {
  const s = SYMS[sym];
  if (!s) return { min: 0, max: 0 };
  if (!qWindow) return { min: s.min, max: s.max };
  const [qmin, qmax] = qWindow;
  if (!(qmin > 0) || !(qmax > qmin)) return { min: s.min, max: s.max };
  const aAtQ = (q: number) =>
    s.kind === "hex"
      ? (4 * Math.PI) / (Math.sqrt(3) * q)
      : (TWO_PI * Math.sqrt(s.Ms[0]!)) / q;
  // union with the envelope: widen only, never narrow.
  return { min: Math.min(aAtQ(qmax), s.min), max: Math.max(aAtQ(qmin), s.max) };
}

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
