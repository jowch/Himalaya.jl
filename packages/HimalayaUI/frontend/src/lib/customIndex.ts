import phaseRatioIndices from "./phaseRatioIndices.json";

const TWO_PI = 2 * Math.PI;

export const CUSTOM_SNAP_TOL = 0.022;

// physics: the 2D square lattice obeys q = 2π√N / a with N = h²+k², the same
// physics: law as the cubics (N = h²+k²+l²), so it takes the cubic branch.
export type SymmetryKind = "cubic" | "lamellar" | "hex" | "square";

export interface SymmetrySpec {
  kind: SymmetryKind;
  Ms: number[];
  param: "a" | "d";
  def: number;
  min: number;
  max: number;
}

// physics: Ms are the allowed h²+k²+l² (3D) or h²+k² (2D) from core Himalaya's
// physics: phaseratios, so a custom comb matches the indexer's ratio series.
// physics: [min, max] is a wide slider envelope in Å: a 10 Å floor, lamellar
// physics: swelling past 100 nm and bicontinuous or micellar cubics up to ~65 nm.
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

// physics: the envelope only widens to what the trace's q window can reach. The
// physics: first order inverts customRefls: a = 2π√(Ms[0])/q for cubic, square and
// physics: lamellar, and a = 4π/(√3·q) for hex. q_max gives the smallest a.
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
  return { min: Math.min(aAtQ(qmax), s.min), max: Math.max(aAtQ(qmin), s.max) };
}

export interface CustomReflection { N: number; q: number }

export function customRefls(sym: string, val: number): CustomReflection[] {
  const s = SYMS[sym];
  if (!s || val <= 0) return [];
  if (s.kind === "lamellar") return s.Ms.map((n) => ({ N: n, q: (TWO_PI * n) / val }));
  if (s.kind === "hex") {
    const q1 = (4 * Math.PI) / (Math.sqrt(3) * val);
    return s.Ms.map((M) => ({ N: M, q: q1 * Math.sqrt(M) }));
  }
  return s.Ms.map((N) => ({ N, q: (TWO_PI * Math.sqrt(N)) / val }));
}

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
