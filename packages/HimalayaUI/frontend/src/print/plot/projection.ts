import { scaleLinear, scaleLog } from "d3-scale";

export type ScaleType = "log" | "linear";

/** One axis: data↔pixel mapping plus "nice" tick values. Backed by d3-scale. */
export interface Axis1D {
  /** data value → pixel. */
  to(value: number): number;
  /** pixel → data value. */
  invert(px: number): number;
  /** "Nice" tick values across the domain (default ~6). */
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

export interface Projection {
  x: Axis1D;
  y: Axis1D;
}

/**
 * Build x (left→right) and y (top→bottom *inverted*: domain max at px 0) axes
 * for a plot area of `plotWidth × plotHeight` pixels.
 */
export function makeProjection(opts: {
  xDomain: [number, number];
  yDomain: [number, number];
  plotWidth: number;
  plotHeight: number;
  xType?: ScaleType;
  yType?: ScaleType;
}): Projection {
  const {
    xDomain,
    yDomain,
    plotWidth,
    plotHeight,
    xType = "log",
    yType = "log",
  } = opts;
  return {
    x: makeAxis(xDomain, [0, plotWidth], xType),
    y: makeAxis(yDomain, [plotHeight, 0], yType),
  };
}

/**
 * Tight [min, max] over the positive, finite values — the domain a log scale
 * needs. Falls back to [1, 10] when there is no positive data, and pads a
 * degenerate single-value extent so the scale is non-empty.
 */
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
