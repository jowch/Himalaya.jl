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

function edgeValue(series: Series, at: number, log: boolean): number | null {
  for (let i = 1; i < series.q.length; i++) {
    const q0 = series.q[i - 1]!;
    const q1 = series.q[i]!;
    if ((q0 - at) * (q1 - at) > 0) continue;
    const a = series.I[i - 1]!;
    const b = series.I[i]!;
    if (!usable(a) || !usable(b) || q1 === q0) return null;
    const t = (at - q0) / (q1 - q0);
    return log ? 10 ** (Math.log10(a) + t * (Math.log10(b) - Math.log10(a))) : a + t * (b - a);
  }
  return null;
}

function windowValues(series: Series, lo: number, hi: number, log: boolean): number[] {
  const inside = series.I.filter((v, i) => series.q[i]! >= lo && series.q[i]! <= hi && usable(v));
  const edges = [edgeValue(series, lo, log), edgeValue(series, hi, log)].filter((v): v is number => v !== null);
  return [...inside, ...edges];
}

export function fitYDomain(main: Series, extra: readonly Series[], xDomain: [number, number], type: ScaleType): [number, number] {
  const lo = Math.min(xDomain[0], xDomain[1]);
  const hi = Math.max(xDomain[0], xDomain[1]);
  const log = type === "log";
  const own = windowValues(main, lo, hi, log);
  const values = own.length < 2 ? main.I.filter(usable) : [...own, ...extra.flatMap((s) => windowValues(s, lo, hi, log))];
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
    if (exps.length >= 3) {
      const stride = Math.ceil(exps.length / (count + 2));
      return exps.filter((_, i) => i % stride === 0).map((e) => ({ value: 10 ** e, label: decadeLabel(e) }));
    }
  }
  return labelled(axisTicks(makeAxis([lo, hi], [0, 1], "log"), count).filter((t) => t.kind !== "minor").map((t) => t.value));
}

export function ticksFor(domain: [number, number], type: ScaleType, count: number, decades: boolean): AxisTick[] {
  return type === "log" ? logTicks(domain, count, decades) : linearTicks(domain, count);
}
