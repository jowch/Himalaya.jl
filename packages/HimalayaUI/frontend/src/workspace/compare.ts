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

/* physics: scale the comparison to the current trace over the top 15% of their
   shared q range. By high q the Bragg peaks have died away and both curves are
   mostly solvent and capillary background, so matching there lines a buffer up
   with the sample's background before it is subtracted. */
export const SCALE_WINDOW = 0.15;

export function matchScale(main: Pick<Trace, "q" | "I">, ref: Pick<Trace, "q" | "I">): number {
  if (main.q.length < 2 || ref.q.length < 2) return 1;
  const lo = Math.max(main.q[0]!, ref.q[0]!);
  const hi = Math.min(main.q[main.q.length - 1]!, ref.q[ref.q.length - 1]!);
  if (!(hi > lo)) return 1;
  const from = hi - (hi - lo) * SCALE_WINDOW;
  const ratios: number[] = [];
  main.q.forEach((q, i) => {
    if (q < from || q > hi) return;
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

export function toggleCompare(refs: readonly number[], id: number): number[] {
  if (refs.includes(id)) return refs.filter((r) => r !== id);
  return [...refs, id].slice(-MAX_COMPARES);
}
