// customIndex.ts (Plan D-9) — pure physics for the custom-index modal.
//
// Symmetry → predicted reflections from a user-chosen lattice. Mirrors the
// mockup `SYMS` / `customRefls`.

import phaseRatioIndices from "./phaseRatioIndices.json";

const TWO_PI = 2 * Math.PI;

export const CUSTOM_SNAP_TOL = 0.022;

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
  Pn3m:      { kind: "cubic",    Ms: phaseRatioIndices.Pn3m,      param: "a", def: 197, min: 10, max: 500 },
  Im3m:      { kind: "cubic",    Ms: phaseRatioIndices.Im3m,    param: "a", def: 252, min: 10, max: 700 },
  Ia3d:      { kind: "cubic",    Ms: phaseRatioIndices.Ia3d,  param: "a", def: 218, min: 10, max: 650 },
  Fm3m:      { kind: "cubic",    Ms: phaseRatioIndices.Fm3m,       param: "a", def: 200, min: 10, max: 500 },
  Fd3m:      { kind: "cubic",    Ms: phaseRatioIndices.Fd3m, param: "a", def: 250, min: 10, max: 600 },
  Lamellar:  { kind: "lamellar", Ms: phaseRatioIndices.Lamellar,         param: "d", def: 60,  min: 10, max: 1000 },
  Hexagonal: { kind: "hex",      Ms: phaseRatioIndices.Hexagonal,     param: "a", def: 70,  min: 10, max: 500 },
  Square:    { kind: "square",   Ms: phaseRatioIndices.Square, param: "a", def: 80, min: 10, max: 500 },
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
