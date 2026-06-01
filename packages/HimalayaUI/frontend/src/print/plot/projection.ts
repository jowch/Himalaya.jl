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

export type TickKind = "major" | "mid" | "minor";
export interface Tick { value: number; kind: TickKind; }

/** Decade-anchored ticks for log axes (major = 10ⁿ, mid = ×2/×5, minor = rest);
 *  linear axes pass d3's nice ticks through as majors.
 *
 *  Density is adaptive so the axis is always legible AND always labelled:
 *   - a moderate range keeps the full detailed set (decades + ×2/×5 + minors);
 *   - a wide range (intensity over many decades) thins to majors, then strides
 *     the decades, so the label count stays near `count`;
 *   - a degenerate sub-decade window (e.g. after a tight wheel-zoom) where no
 *     decade multiple lands in range falls back to d3's nice ticks, and finally
 *     to the geometric endpoints — never an empty, unlabelled axis. */
export function axisTicks(axis: Axis1D, count = 6): Tick[] {
  if (axis.type === "linear") {
    return axis.ticks(count).map((value) => ({ value, kind: "major" as const }));
  }
  const lo = Math.min(axis.domain[0], axis.domain[1]);
  const hi = Math.max(axis.domain[0], axis.domain[1]);
  if (!(lo > 0) || !(hi > 0)) return [];

  const inRange = (v: number): boolean =>
    v >= lo * (1 - 1e-9) && v <= hi * (1 + 1e-9);
  const eLo = Math.floor(Math.log10(lo));
  const eHi = Math.floor(Math.log10(hi));

  const detailed: Tick[] = [];
  for (let e = eLo; e <= eHi; e++) {
    const decade = Math.pow(10, e);
    for (const m of [1, 2, 3, 4, 5, 6, 7, 8, 9]) {
      const value = m * decade;
      if (!inRange(value)) continue;
      const kind: TickKind = m === 1 ? "major" : m === 2 || m === 5 ? "mid" : "minor";
      detailed.push({ value, kind });
    }
  }

  // Sub-decade window with no decade multiple in range → never leave the axis
  // bare: use d3's nice log ticks, and the geometric endpoints as a last resort.
  if (detailed.length === 0) {
    const nice = axis.ticks(count);
    if (nice.length > 0) {
      return nice.map((value) => ({ value, kind: "major" as const }));
    }
    return [lo, Math.sqrt(lo * hi), hi].map((value) => ({
      value,
      kind: "major" as const,
    }));
  }

  const cap = count + 2;
  const labelled = detailed.filter((t) => t.kind !== "minor").length;
  if (labelled <= cap) return detailed;

  // Too dense → decades only, strided to keep ~`count` labels.
  const majors = detailed.filter((t) => t.kind === "major");
  if (majors.length === 0) {
    return detailed.filter((t) => t.kind === "mid"); // sub-decade-ish: keep ×2/×5
  }
  if (majors.length <= cap) return majors;
  const stride = Math.ceil(majors.length / count);
  return majors.filter((_, i) => i % stride === 0);
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
