import { stripRange } from "./compare";
import { fitYDomain, linearTicks, type AxisTick } from "./plotScale";
import { viewAxes, viewIntensity, type PlotView } from "./plotView";
import type { ScaleType } from "../plot/projection";
import { formatTicks } from "../lib/plot/formatAxis";

export interface DifferenceAxis {
  kind: "linear" | "pseudolog";
  to: (v: number) => number;
  domain: [number, number];
  ticks: AxisTick[];
}

// physics: a difference crosses zero, so a log axis cannot hold it. asinh(d / c) is
// physics: linear for |d| below c and logarithmic above, so it keeps a log axis's
// physics: reach for peaks without hiding the sign. c is the noise of the
// physics: difference (1.4826 × the median absolute deviation, the Gaussian σ) over
// physics: the upper half of the visible q range, where the residual has died away,
// physics: so noise sits in the linear band and real peaks in the logarithmic part.
const MAD_TO_SIGMA = 1.4826;

export function noiseScale(q: readonly number[], d: readonly number[], domain: [number, number]): number {
  const from = domain[0] + (domain[1] - domain[0]) / 2;
  const vals = d.filter((v, i) => q[i]! >= from && q[i]! <= domain[1] && Number.isFinite(v));
  const median = (xs: number[]): number => {
    const s = [...xs].sort((a, b) => a - b);
    return s.length === 0 ? 0 : s[Math.floor(s.length / 2)]!;
  };
  const m = median(vals);
  const sigma = MAD_TO_SIGMA * median(vals.map((v) => Math.abs(v - m)));
  if (sigma > 0) return sigma;
  const span = Math.max(...d.filter(Number.isFinite).map(Math.abs), 0);
  return span > 0 ? span / 100 : 1;
}

function padded(lo0: number, hi0: number): [number, number] {
  const lo = Math.min(lo0, 0);
  const hi = Math.max(hi0, 0);
  const pad = (hi - lo) * 0.1 || 1;
  return [lo - pad, hi + pad];
}

function pseudologTicks(domain: [number, number], c: number, to: (v: number) => number): AxisTick[] {
  const values: number[] = [0];
  const maxAbs = Math.max(Math.abs(Math.sinh(domain[0]) * c), Math.abs(Math.sinh(domain[1]) * c));
  for (let e = Math.ceil(Math.log10(c)); e <= Math.floor(Math.log10(maxAbs)); e++) values.push(10 ** e, -(10 ** e));
  const inside = values.filter((v) => to(v) >= domain[0] && to(v) <= domain[1]).sort((a, b) => a - b);
  const pos = inside.filter((v) => v > 0);
  const neg = inside.filter((v) => v < 0);
  const thin = (xs: number[]): number[] => {
    const stride = Math.ceil(xs.length / 3);
    return xs.filter((_, i) => (xs.length - 1 - i) % stride === 0);
  };
  const shown = [...thin(neg.reverse()).reverse(), 0, ...thin(pos)].filter((v) => inside.includes(v));
  const labels = formatTicks(shown);
  return shown.map((v, i) => ({ value: to(v), label: labels[i]! }));
}

export function differenceAxis(
  view: PlotView,
  q: readonly number[],
  d: readonly number[],
  shownX: [number, number],
  peakQs: readonly number[],
): DifferenceAxis {
  if (viewAxes(view).yType !== "log") {
    const domain = padded(...stripRange(q, d, shownX, peakQs));
    return { kind: "linear", to: (v) => v, domain, ticks: linearTicks(domain, 5) };
  }
  const c = noiseScale(q, d, shownX);
  const to = (v: number): number => Math.asinh(v / c);
  const domain = padded(...stripRange(q, d.map(to), shownX, peakQs));
  return { kind: "pseudolog", to, domain, ticks: pseudologTicks(domain, c, to) };
}

interface Subtracted {
  q: readonly number[];
  d: readonly number[];
}

export function subtractedSeries(
  view: PlotView,
  q: readonly number[],
  d: readonly number[],
  shownX: [number, number],
  others: readonly Subtracted[] = [],
): { I: number[]; yType: ScaleType; domain: [number, number]; others: number[][] } {
  const yType = viewAxes(view).yType;
  const viewed = (s: Subtracted): number[] => s.d.map((v, i) => viewIntensity(view, s.q[i]!, v));
  const I = viewed({ q, d });
  const otherI = others.map(viewed);
  const extra = others.map((s, i) => ({ q: s.q, I: otherI[i]! }));
  return { I, yType, domain: fitYDomain({ q, I }, extra, shownX, yType, true), others: otherI };
}
