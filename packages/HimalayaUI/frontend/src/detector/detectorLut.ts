export interface LutStop { t: number; L: number; C: number; h: number }

export type DetectorLutVariant = "neutral" | "warm" | "magma";

export const MAGMA_SRGB_STOPS: readonly (readonly [number, number, number])[] = [
  [0, 0, 4], [40, 11, 84], [101, 21, 110], [159, 42, 99],
  [212, 72, 66], [245, 125, 21], [250, 193, 39], [252, 255, 164],
];

function buildSrgbRamp(stops: readonly (readonly [number, number, number])[]): Uint8ClampedArray {
  const lut = new Uint8ClampedArray(256 * 3);
  const n = stops.length - 1;
  for (let i = 0; i < 256; i++) {
    const x = (i / 255) * n;
    const k = Math.min(n - 1, Math.floor(x));
    const f = x - k;
    const a = stops[k]!, b = stops[k + 1]!;
    for (let c = 0; c < 3; c++) lut[i * 3 + c] = Math.round(a[c]! + f * (b[c]! - a[c]!));
  }
  return lut;
}

export const PRINT_DETECTOR_STOPS_NEUTRAL: readonly LutStop[] = [
  { t: 0.00, L: 0.130, C: 0.005, h: 60 },
  { t: 0.30, L: 0.320, C: 0.008, h: 65 },
  { t: 0.60, L: 0.550, C: 0.009, h: 70 },
  { t: 0.82, L: 0.780, C: 0.007, h: 78 },
  { t: 1.00, L: 0.950, C: 0.005, h: 85 },
];

export const PRINT_DETECTOR_STOPS_WARM: readonly LutStop[] = [
  { t: 0.00, L: 0.150, C: 0.010, h: 55 },
  { t: 0.35, L: 0.320, C: 0.060, h: 50 },
  { t: 0.62, L: 0.550, C: 0.130, h: 45 },
  { t: 0.85, L: 0.780, C: 0.090, h: 70 },
  { t: 1.00, L: 0.930, C: 0.020, h: 85 },
];

export const PRINT_DETECTOR_STOPS = PRINT_DETECTOR_STOPS_NEUTRAL;

function stopsFor(variant: DetectorLutVariant): readonly LutStop[] {
  return variant === "warm" ? PRINT_DETECTOR_STOPS_WARM : PRINT_DETECTOR_STOPS_NEUTRAL;
}

function clamp01(v: number): number { return v < 0 ? 0 : v > 1 ? 1 : v; }

function oklchToSrgb(L: number, C: number, h: number): [number, number, number] {
  const hr = (h * Math.PI) / 180;
  const a = C * Math.cos(hr);
  const b = C * Math.sin(hr);
  const l_ = L + 0.3963377774 * a + 0.2158037573 * b;
  const m_ = L - 0.1055613458 * a - 0.0638541728 * b;
  const s_ = L - 0.0894841775 * a - 1.2914855480 * b;
  const l = l_ * l_ * l_, m = m_ * m_ * m_, s = s_ * s_ * s_;
  const lr =  4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s;
  const lg = -1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s;
  const lb = -0.0041960863 * l - 0.7034186147 * m + 1.7076147010 * s;
  const gamma = (v: number): number => {
    const x = v <= 0.0031308 ? 12.92 * v : 1.055 * Math.pow(v, 1 / 2.4) - 0.055;
    return clamp01(x);
  };
  return [
    Math.round(gamma(lr) * 255),
    Math.round(gamma(lg) * 255),
    Math.round(gamma(lb) * 255),
  ];
}

const cache: Partial<Record<DetectorLutVariant, Uint8ClampedArray>> = {};

export const DETECTOR_LUT_VARIANTS: readonly DetectorLutVariant[] = ["neutral", "warm", "magma"];

export function detectorLutFilterId(variant: DetectorLutVariant): string {
  return `detector-lut-${variant}`;
}

export function detectorLutFilter(variant: DetectorLutVariant): string {
  return `url(#${detectorLutFilterId(variant)})`;
}

export function buildPrintDetectorLut(variant: DetectorLutVariant = "neutral"): Uint8ClampedArray {
  const hit = cache[variant];
  if (hit) return hit;
  if (variant === "magma") return (cache.magma = buildSrgbRamp(MAGMA_SRGB_STOPS));
  const lut = new Uint8ClampedArray(256 * 3);
  const stops = stopsFor(variant);
  for (let i = 0; i < 256; i++) {
    const t = i / 255;
    let a = stops[0], b = stops[stops.length - 1];
    for (let k = 0; k < stops.length - 1; k++) {
      if (t >= stops[k].t && t <= stops[k + 1].t) { a = stops[k]; b = stops[k + 1]; break; }
    }
    const span = b.t - a.t || 1;
    const f = (t - a.t) / span;
    const [r, g, bl] = oklchToSrgb(
      a.L + f * (b.L - a.L),
      a.C + f * (b.C - a.C),
      a.h + f * (b.h - a.h),
    );
    lut[i * 3] = r; lut[i * 3 + 1] = g; lut[i * 3 + 2] = bl;
  }
  cache[variant] = lut;
  return lut;
}

export function detectorLutTableValues(
  variant: DetectorLutVariant = "neutral",
): { r: string; g: string; b: string } {
  const lut = buildPrintDetectorLut(variant);
  const r: string[] = [], g: string[] = [], b: string[] = [];
  for (let i = 0; i < 256; i++) {
    r.push((lut[i * 3] / 255).toFixed(4));
    g.push((lut[i * 3 + 1] / 255).toFixed(4));
    b.push((lut[i * 3 + 2] / 255).toFixed(4));
  }
  return { r: r.join(" "), g: g.join(" "), b: b.join(" ") };
}