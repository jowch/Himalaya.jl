import { scaleLinear } from "d3-scale";
import { axisTicks, makeAxis, type ScaleType } from "../print/plot/projection";
import { formatTicks } from "../lib/plot/formatAxis";

export interface Series {
  q: readonly number[];
  I: readonly number[];
}

export interface AxisTick {
  value: number;
  label: string;
}

export const PLOT_MARGINS = { top: 22, right: 12, bottom: 26, left: 52 } as const;

function usable(v: number): boolean {
  return Number.isFinite(v) && v > 0;
}

function windowValues(series: Series, lo: number, hi: number): number[] {
  const out: number[] = [];
  let below = -1;
  let above = -1;
  let reachesLo = false;
  let reachesHi = false;
  series.q.forEach((q, i) => {
    if (q >= lo) reachesLo = true;
    if (q <= hi) reachesHi = true;
    if (q >= lo && q <= hi) out.push(series.I[i]!);
    else if (q < lo && (below < 0 || q > series.q[below]!)) below = i;
    else if (q > hi && (above < 0 || q < series.q[above]!)) above = i;
  });
  if (below >= 0 && reachesLo) out.push(series.I[below]!);
  if (above >= 0 && reachesHi) out.push(series.I[above]!);
  return out.filter(usable);
}

export function fitYDomain(main: Series, extra: readonly Series[], xDomain: [number, number], type: ScaleType): [number, number] {
  const lo = Math.min(xDomain[0], xDomain[1]);
  const hi = Math.max(xDomain[0], xDomain[1]);
  const own = windowValues(main, lo, hi);
  const values = own.length < 2 ? main.I.filter(usable) : [...own, ...extra.flatMap((s) => windowValues(s, lo, hi))];
  if (values.length === 0) return type === "log" ? [1, 10] : [0, 1];
  const min = Math.min(...values);
  const max = Math.max(...values);
  if (type === "log") {
    const a = Math.log10(min);
    const b = Math.log10(max);
    const pad = (b - a) * 0.06 || 0.1;
    return [10 ** (a - pad), 10 ** (b + pad * 1.6)];
  }
  const pad = (max - min) * 0.06 || Math.abs(max) * 0.1 || 1;
  return [min - pad, max + pad * 1.6];
}

const SUPERSCRIPT: Record<string, string> = { "-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹" };

export function decadeLabel(exponent: number): string {
  return "10" + String(exponent).split("").map((c) => SUPERSCRIPT[c] ?? c).join("");
}

function labelled(values: number[]): AxisTick[] {
  const labels = formatTicks(values);
  return values.map((value, i) => ({ value, label: labels[i]! }));
}

export function linearTicks(domain: [number, number], count: number): AxisTick[] {
  return labelled(scaleLinear().domain(domain).ticks(count));
}

export function logTicks(domain: [number, number], count: number, decades: boolean): AxisTick[] {
  const lo = Math.min(domain[0], domain[1]);
  const hi = Math.max(domain[0], domain[1]);
  if (!(lo > 0) || !(hi > 0)) return [];
  if (decades) {
    const exps: number[] = [];
    for (let e = Math.ceil(Math.log10(lo) - 1e-9); e <= Math.floor(Math.log10(hi) + 1e-9); e++) exps.push(e);
    if (exps.length >= 2) {
      const stride = Math.ceil(exps.length / (count + 2));
      return exps.filter((_, i) => i % stride === 0).map((e) => ({ value: 10 ** e, label: decadeLabel(e) }));
    }
  }
  return labelled(axisTicks(makeAxis([lo, hi], [0, 1], "log"), count).filter((t) => t.kind !== "minor").map((t) => t.value));
}

export function ticksFor(domain: [number, number], type: ScaleType, count: number, decades: boolean): AxisTick[] {
  return type === "log" ? logTicks(domain, count, decades) : linearTicks(domain, count);
}
