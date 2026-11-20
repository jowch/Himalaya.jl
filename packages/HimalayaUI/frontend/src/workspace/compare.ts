import type { Trace } from "../api";

export const MAX_COMPARES = 3;

export function interpolateAt(trace: Pick<Trace, "q" | "I">, q: number): number {
  const xs = trace.q;
  const n = xs.length;
  if (n === 0 || !(q >= xs[0]!) || !(q <= xs[n - 1]!)) return Number.NaN;
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid]! <= q) lo = mid;
    else hi = mid;
  }
  const x0 = xs[lo]!;
  const x1 = xs[hi]!;
  const y0 = trace.I[lo]!;
  const y1 = trace.I[hi]!;
  return x1 === x0 ? y0 : y0 + ((y1 - y0) * (q - x0)) / (x1 - x0);
}

// physics: the comparison is scaled to the current trace by the median intensity
// physics: ratio over the top 15% of their shared q range, on the assumption that
// physics: by high q both curves are mostly solvent and capillary background.
// physics: Bragg reflections often still sit in that window (HII and higher Lα
// physics: orders), so q within PEAK_MASK of the current trace's own peaks is left
// physics: out. A sample's own smooth high-q scattering (bilayer form factor, more
// physics: lipid) still reads as extra background and over-scales the reference,
// physics: so the buffer is over-subtracted. Scaling by transmitted flux and
// physics: exposure time from the beam monitors is the right fix.
export const SCALE_WINDOW = 0.15;
export const PEAK_MASK = 0.008;

export function matchScale(
  main: Pick<Trace, "q" | "I">,
  ref: Pick<Trace, "q" | "I">,
  peakQs: readonly number[] = [],
): number {
  if (main.q.length < 2 || ref.q.length < 2) return 1;
  const lo = Math.max(main.q[0]!, ref.q[0]!);
  const hi = Math.min(main.q[main.q.length - 1]!, ref.q[ref.q.length - 1]!);
  if (!(hi > lo)) return 1;
  const from = hi - (hi - lo) * SCALE_WINDOW;
  const ratios: number[] = [];
  main.q.forEach((q, i) => {
    if (q < from || q > hi) return;
    if (peakQs.some((p) => Math.abs(q - p) < PEAK_MASK)) return;
    const r = interpolateAt(ref, q);
    const m = main.I[i]!;
    if (r > 0 && m > 0) ratios.push(m / r);
  });
  if (ratios.length === 0) return 1;
  ratios.sort((a, b) => a - b);
  return ratios[Math.floor(ratios.length / 2)]!;
}

export function scaledTrace(ref: Pick<Trace, "q" | "I">, scale: number): { q: number[]; I: number[] } {
  return { q: [...ref.q], I: ref.I.map((v) => v * scale) };
}

export function subtractTrace(
  main: Pick<Trace, "q" | "I">,
  ref: Pick<Trace, "q" | "I">,
  scale: number,
): { q: number[]; d: number[] } {
  const q: number[] = [];
  const d: number[] = [];
  main.q.forEach((x, i) => {
    const r = interpolateAt(ref, x);
    if (Number.isFinite(r)) {
      q.push(x);
      d.push(main.I[i]! - scale * r);
    }
  });
  return { q, d };
}

export function addCompare(refs: readonly number[], id: number): number[] {
  if (refs.includes(id)) return [...refs];
  return [...refs, id].slice(-MAX_COMPARES);
}

export function robustRange(values: readonly number[], lowFrac = 0.01, highFrac = 0.99): [number, number] {
  const finite = values.filter(Number.isFinite).sort((a, b) => a - b);
  if (finite.length === 0) return [-1, 1];
  const at = (f: number): number => finite[Math.min(finite.length - 1, Math.max(0, Math.round(f * (finite.length - 1))))]!;
  const lo = Math.min(0, at(lowFrac));
  const hi = Math.max(0, at(highFrac));
  return hi > lo ? [lo, hi] : [lo - 1, hi + 1];
}

export function clampDomain(
  domain: [number, number] | null,
  extent: [number, number] | null,
): [number, number] | null {
  if (!domain || !extent) return null;
  const lo = Math.max(domain[0], extent[0]);
  const hi = Math.min(domain[1], extent[1]);
  return hi > lo ? [lo, hi] : null;
}
