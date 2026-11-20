import { customRefls, landsOn, latticeBounds, SYMS } from "../lib/customIndex";

export const MATCH_TOL = 0.022;

export interface IndexPoint {
  ratio: number;
  predicted: number;
  observed: number | null;
}

export interface IndexFit {
  points: IndexPoint[];
  slope: number;
  r2: number | null;
  landed: number;
}

export function indexFit(predicted: readonly number[], peakQs: readonly number[], tol = MATCH_TOL): IndexFit {
  const q1 = predicted[0] ?? 0;
  const used = new Set<number>();
  const points: IndexPoint[] = predicted.map((q) => {
    let best: number | null = null;
    let bestD = tol;
    peakQs.forEach((p, i) => {
      const d = Math.abs(p - q) / q;
      if (d < bestD && !used.has(i)) {
        bestD = d;
        best = i;
      }
    });
    if (best !== null) used.add(best);
    return { ratio: q1 > 0 ? q / q1 : 0, predicted: q, observed: best === null ? null : peakQs[best]! };
  });
  const hit = points.filter((p): p is IndexPoint & { observed: number } => p.observed !== null);
  const sxy = hit.reduce((a, p) => a + p.ratio * p.observed, 0);
  const sxx = hit.reduce((a, p) => a + p.ratio * p.ratio, 0);
  const slope = sxx > 0 ? sxy / sxx : q1;
  let r2: number | null = null;
  if (hit.length >= 2) {
    const mean = hit.reduce((a, p) => a + p.observed, 0) / hit.length;
    const ssr = hit.reduce((a, p) => a + (p.observed - slope * p.ratio) ** 2, 0);
    const sst = hit.reduce((a, p) => a + (p.observed - mean) ** 2, 0);
    r2 = sst > 0 ? 1 - ssr / sst : null;
  }
  return { points, slope, r2, landed: hit.length };
}

export const CUSTOM_PHASES = Object.keys(SYMS);

export interface CustomDraft {
  phase: string;
  lattice: number;
}

export function clampLattice(phase: string, value: number, qWindow: readonly [number, number] | null): number {
  const { min, max } = latticeBounds(phase, qWindow);
  return Math.round(Math.min(max, Math.max(min, value)) * 10) / 10;
}

export function draftPrediction(d: CustomDraft): number[] {
  return customRefls(d.phase, d.lattice).map((r) => r.q);
}

export function draftLanded(d: CustomDraft, peakQs: number[]): number {
  return landsOn(d.phase, d.lattice, peakQs);
}

export function draftCommit(d: CustomDraft): { phase: string; basis: number; ratios: number[] } {
  const qs = draftPrediction(d);
  const q1 = qs[0] ?? 0;
  return { phase: d.phase, basis: q1, ratios: q1 > 0 ? qs.map((q) => q / q1) : [] };
}

export function initialDraft(top: { phase: string; lattice_d: number | null } | undefined): CustomDraft {
  if (top && SYMS[top.phase] && top.lattice_d !== null) return { phase: top.phase, lattice: Math.round(top.lattice_d * 10) / 10 };
  return { phase: "Pn3m", lattice: SYMS.Pn3m!.def };
}
