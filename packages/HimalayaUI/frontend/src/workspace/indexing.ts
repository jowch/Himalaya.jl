import type { IndexEntry } from "../api";
import { CUSTOM_SNAP_TOL, customRefls, latticeBounds, SYMS } from "../lib/customIndex";

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

type QWindow = readonly [number, number] | null;

function inWindow(q: number, window: QWindow): boolean {
  return window === null || (q >= window[0] && q <= window[1]);
}

export function claimPeaks(predicted: readonly number[], peakQs: readonly number[], tol = CUSTOM_SNAP_TOL): (number | null)[] {
  const best = predicted.map((q) => {
    let pick: { peak: number; rel: number } | null = null;
    peakQs.forEach((p, i) => {
      const rel = Math.abs(p - q) / q;
      if (rel <= tol && (pick === null || rel < pick.rel)) pick = { peak: i, rel };
    });
    return pick as { peak: number; rel: number } | null;
  });
  const out: (number | null)[] = predicted.map(() => null);
  const claimed = new Set<number>();
  best
    .map((b, order) => ({ b, order }))
    .filter((x): x is { b: { peak: number; rel: number }; order: number } => x.b !== null)
    .sort((x, y) => x.b.rel - y.b.rel)
    .forEach(({ b, order }) => {
      if (claimed.has(b.peak)) return;
      claimed.add(b.peak);
      out[order] = peakQs[b.peak]!;
    });
  return out;
}

function fitPoints(points: IndexPoint[], fallbackSlope: number): IndexFit {
  const hit = points.filter((p): p is IndexPoint & { observed: number } => p.observed !== null);
  const sxy = hit.reduce((a, p) => a + p.ratio * p.observed, 0);
  const sxx = hit.reduce((a, p) => a + p.ratio * p.ratio, 0);
  const slope = sxx > 0 ? sxy / sxx : fallbackSlope;
  let r2: number | null = null;
  if (hit.length >= 2) {
    const mean = hit.reduce((a, p) => a + p.observed, 0) / hit.length;
    const ssr = hit.reduce((a, p) => a + (p.observed - slope * p.ratio) ** 2, 0);
    const sst = hit.reduce((a, p) => a + (p.observed - mean) ** 2, 0);
    r2 = sst > 0 ? 1 - ssr / sst : null;
  }
  return { points, slope, r2, landed: hit.length };
}

export function indexFit(predicted: readonly number[], peakQs: readonly number[], window: QWindow = null, tol = CUSTOM_SNAP_TOL): IndexFit {
  const q1 = predicted[0] ?? 0;
  const observed = claimPeaks(predicted, peakQs, tol);
  const points = predicted
    .map((q, i) => ({ ratio: q1 > 0 ? q / q1 : 0, predicted: q, observed: observed[i] ?? null }))
    .filter((p) => p.observed !== null || inWindow(p.predicted, window));
  return fitPoints(points, q1);
}

export function storedFit(ix: Pick<IndexEntry, "predicted_q" | "peaks" | "basis">, window: QWindow = null): IndexFit {
  const predicted = ix.predicted_q;
  const q1 = predicted[0] ?? 0;
  const observed: (number | null)[] = predicted.map(() => null);
  for (const ref of ix.peaks) {
    const pos = ref.ratio_position as unknown;
    let order = typeof pos === "number" && Number.isInteger(pos) && pos >= 1 && pos <= predicted.length ? pos - 1 : -1;
    if (order < 0) {
      let bestD = Infinity;
      predicted.forEach((q, i) => {
        const d = Math.abs(q - ref.q_observed);
        if (d < bestD) {
          bestD = d;
          order = i;
        }
      });
    }
    if (order >= 0) observed[order] = ref.q_observed;
  }
  const points = predicted
    .map((q, i) => ({ ratio: q1 > 0 ? q / q1 : 0, predicted: q, observed: observed[i] ?? null }))
    .filter((p) => p.observed !== null || inWindow(p.predicted, window));
  return fitPoints(points, ix.basis);
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

export function draftCommit(d: CustomDraft): { phase: string; basis: number; ratios: number[] } {
  const qs = draftPrediction(d);
  const q1 = qs[0] ?? 0;
  return { phase: d.phase, basis: q1, ratios: q1 > 0 ? qs.map((q) => q / q1) : [] };
}

export function initialDraft(top: { phase: string; lattice_d: number | null } | undefined): CustomDraft {
  if (top && SYMS[top.phase] && top.lattice_d !== null) return { phase: top.phase, lattice: Math.round(top.lattice_d * 10) / 10 };
  return { phase: "Pn3m", lattice: SYMS.Pn3m!.def };
}
