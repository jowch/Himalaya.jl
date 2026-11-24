import { scaleLinear } from "d3-scale";
import type { ScaleType } from "../plot/projection";
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

function finite(v: number): boolean {
  return Number.isFinite(v);
}

function edgeValue(series: Series, at: number, log: boolean, ok: (v: number) => boolean): number | null {
  for (let i = 1; i < series.q.length; i++) {
    const q0 = series.q[i - 1]!;
    const q1 = series.q[i]!;
    if (!finite(q0) || !finite(q1) || (q0 - at) * (q1 - at) > 0) continue;
    const a = series.I[i - 1]!;
    const b = series.I[i]!;
    if (!ok(a) || !ok(b) || q1 === q0) return null;
    const t = (at - q0) / (q1 - q0);
    return log ? 10 ** (Math.log10(a) + t * (Math.log10(b) - Math.log10(a))) : a + t * (b - a);
  }
  return null;
}

function windowValues(series: Series, lo: number, hi: number, log: boolean, ok: (v: number) => boolean): number[] {
  const inside = series.I.filter((v, i) => series.q[i]! >= lo && series.q[i]! <= hi && ok(v));
  const edges = [edgeValue(series, lo, log, ok), edgeValue(series, hi, log, ok)].filter((v): v is number => v !== null);
  return [...inside, ...edges];
}

export function fitYDomain(main: Series, extra: readonly Series[], xDomain: [number, number], type: ScaleType, signed = false): [number, number] {
  const lo = Math.min(xDomain[0], xDomain[1]);
  const hi = Math.max(xDomain[0], xDomain[1]);
  const log = type === "log";
  const ok = signed && !log ? finite : usable;
  const own = windowValues(main, lo, hi, log, ok);
  const values = own.length < 2 ? main.I.filter(ok) : [...own, ...extra.flatMap((s) => windowValues(s, lo, hi, log, ok))];
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
  const values = scaleLinear().domain(domain).ticks(count);
  if (values.every((v) => Number.isInteger(v))) return values.map((value) => ({ value, label: value.toFixed(0) }));
  return labelled(values);
}

const MANTISSAS = [[1, 2, 5], [1, 2, 3, 4, 5, 6, 7, 8, 9]] as const;

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
  for (const set of MANTISSAS) {
    const values: number[] = [];
    for (let e = Math.floor(Math.log10(lo)); e <= Math.ceil(Math.log10(hi)); e++) {
      for (const m of set) {
        const v = Number((m * 10 ** e).toPrecision(12));
        if (v >= lo && v <= hi) values.push(v);
      }
    }
    if (values.length >= 3) return labelled(values.length > count + 3 ? values.filter((_, i) => i % Math.ceil(values.length / (count + 1)) === 0) : values);
  }
  return linearTicks([lo, hi], count);
}

export function ticksFor(domain: [number, number], type: ScaleType, count: number, decades: boolean): AxisTick[] {
  return type === "log" ? logTicks(domain, count, decades) : linearTicks(domain, count);
}
