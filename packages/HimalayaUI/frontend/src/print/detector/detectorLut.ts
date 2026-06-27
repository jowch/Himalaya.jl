// The Print detector colormap. SAXS images arrive log-normalized grayscale; this
// maps intensity → an OKLCH ramp, converted to sRGB in code (Ottosson) so the LUT
// is pure + JSDOM-testable.
//
// TWO variants:
//  - "neutral" (DEFAULT): a warm-GRAY ink ramp (near-zero chroma, black → warm
//    paper-cream). The detector image is raw MEASUREMENT, so it stays neutral and
//    the phase-coloured ring overlay owns all the hue — maximises ring legibility
//    for alignment / missing-peak judgement, is colour-vision-deficiency safe
//    (the image carries luminance only), and follows the Print's "colour = meaning"
//    rule. This is what the alignment panel uses.
//  - "warm": a saturated terracotta ramp (black → umber → terracotta → gold →
//    cream) — a film-exposure "beauty" look for an un-annotated detector view.
//    Kept available but NOT the default, because saturated warm data clashes with
//    the (often warm-hued) phase rings.

export interface LutStop { t: number; L: number; C: number; h: number } // OKLCH

export type DetectorLutVariant = "neutral" | "warm";

// Warm-GRAY ink ramp: lightness rises monotonically; chroma is tiny (~0.005–0.009)
// so it reads as a warm neutral, leaving the colour channel to the rings.
export const PRINT_DETECTOR_STOPS_NEUTRAL: readonly LutStop[] = [
  { t: 0.00, L: 0.130, C: 0.005, h: 60 }, // warm near-black (window backing)
  { t: 0.30, L: 0.320, C: 0.008, h: 65 },
  { t: 0.60, L: 0.550, C: 0.009, h: 70 },
  { t: 0.82, L: 0.780, C: 0.007, h: 78 },
  { t: 1.00, L: 0.950, C: 0.005, h: 85 }, // warm paper-cream
];

// Saturated terracotta ramp (the original "beauty" look). Hues ~45–85° so no
// hue-wrap handling is needed.
export const PRINT_DETECTOR_STOPS_WARM: readonly LutStop[] = [
  { t: 0.00, L: 0.150, C: 0.010, h: 55 }, // warm near-black (window backing)
  { t: 0.35, L: 0.320, C: 0.060, h: 50 }, // deep umber
  { t: 0.62, L: 0.550, C: 0.130, h: 45 }, // terracotta bloom (the Print accent heat)
  { t: 0.85, L: 0.780, C: 0.090, h: 70 }, // warm gold
  { t: 1.00, L: 0.930, C: 0.020, h: 85 }, // cream highlight
];

// Default stops (the neutral analysis ramp). Exposed for the Storybook swatch.
export const PRINT_DETECTOR_STOPS = PRINT_DETECTOR_STOPS_NEUTRAL;

function stopsFor(variant: DetectorLutVariant): readonly LutStop[] {
  return variant === "warm" ? PRINT_DETECTOR_STOPS_WARM : PRINT_DETECTOR_STOPS_NEUTRAL;
}

function clamp01(v: number): number { return v < 0 ? 0 : v > 1 ? 1 : v; }

// OKLCH -> sRGB (0..255). Bjorn Ottosson's oklab matrices + sRGB gamma.
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

/** Build (and memoize) the 256×RGB Print detector LUT. Defaults to the neutral
 *  warm-gray ramp; pass "warm" for the saturated beauty ramp. */
export function buildPrintDetectorLut(variant: DetectorLutVariant = "neutral"): Uint8ClampedArray {
  const hit = cache[variant];
  if (hit) return hit;
  const lut = new Uint8ClampedArray(256 * 3);
  const stops = stopsFor(variant);
  for (let i = 0; i < 256; i++) {
    const t = i / 255;
    // Find the bracketing stop pair.
    let a = stops[0], b = stops[stops.length - 1];
    for (let k = 0; k < stops.length - 1; k++) {
      if (t >= stops[k].t && t <= stops[k + 1].t) { a = stops[k]; b = stops[k + 1]; break; }
    }
    const span = b.t - a.t || 1;
    const f = (t - a.t) / span;
    // Linear in L/C/h - the stops sit in a narrow warm hue band so no hue
    // wrap-around handling is needed.
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

/**
 * The same 256-entry LUT expressed as SVG `<feFunc* type="table">` ramps — three
 * space-separated 0..1 strings (R, G, B). Lets the browser apply the colormap as
 * a GPU `<feComponentTransfer>` filter on the grayscale PNG instead of a
 * per-pixel JS loop + canvas readback. Derived from `buildPrintDetectorLut`, so
 * the filter and the (retired) loop share one colour source — no drift.
 *
 * Input is grayscale (R=G=B=t), so each output channel is just its table indexed
 * by t — exactly `lut[t]`.
 */
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
