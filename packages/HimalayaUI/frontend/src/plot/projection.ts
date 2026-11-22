import { scaleLinear, scaleLog } from "d3-scale";

export type ScaleType = "log" | "linear";

export interface Axis1D {
  to(value: number): number;
  invert(px: number): number;
  ticks(count?: number): number[];
  readonly domain: readonly [number, number];
  readonly range: readonly [number, number];
  readonly type: ScaleType;
}

export function makeAxis(
  domain: [number, number],
  range: [number, number],
  type: ScaleType,
): Axis1D {
  const scale =
    type === "log" ? scaleLog<number, number>() : scaleLinear<number, number>();
  scale.domain(domain);
  scale.range(range);
  return {
    to: (v: number) => scale(v),
    invert: (px: number) => scale.invert(px),
    ticks: (count = 6) => scale.ticks(count),
    domain,
    range,
    type,
  };
}

export function positiveExtent(
  values: number[],
  fallback: [number, number] = [1, 10],
): [number, number] {
  let min = Infinity;
  let max = -Infinity;
  for (const v of values) {
    if (!Number.isFinite(v) || v <= 0) continue;
    if (v < min) min = v;
    if (v > max) max = v;
  }
  if (min === Infinity) return fallback;
  if (min === max) return [min * 0.9, max * 1.1];
  return [min, max];
}