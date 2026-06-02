// The Print detector colormap. SAXS images arrive log-normalized grayscale;
// this warm OKLCH ramp (black -> umber -> terracotta bloom -> gold -> cream) is the
// Print-native replacement for the old two-endpoint frame-edge->frame-signal lerp.
// OKLCH->sRGB is done in code (Ottosson) so the LUT is pure + JSDOM-testable.

export interface LutStop { t: number; L: number; C: number; h: number } // OKLCH

// Anchored to the Print warm axis (hues ~45-85 deg = terracotta->gold). Lightness
// rises monotonically so the ramp is non-inverting.
export const PRINT_DETECTOR_STOPS: readonly LutStop[] = [
  { t: 0.00, L: 0.150, C: 0.010, h: 55 }, // warm near-black (window backing)
  { t: 0.35, L: 0.320, C: 0.060, h: 50 }, // deep umber
  { t: 0.62, L: 0.550, C: 0.130, h: 45 }, // terracotta bloom (the Print accent heat)
  { t: 0.85, L: 0.780, C: 0.090, h: 70 }, // warm gold
  { t: 1.00, L: 0.930, C: 0.020, h: 85 }, // cream highlight
];

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

let cached: Uint8ClampedArray | null = null;

/** Build (and memoize) the 256xRGB Print detector LUT. */
export function buildPrintDetectorLut(): Uint8ClampedArray {
  if (cached) return cached;
  const lut = new Uint8ClampedArray(256 * 3);
  const stops = PRINT_DETECTOR_STOPS;
  for (let i = 0; i < 256; i++) {
    const t = i / 255;
    // Find the bracketing stop pair.
    let a = stops[0], b = stops[stops.length - 1];
    for (let k = 0; k < stops.length - 1; k++) {
      if (t >= stops[k].t && t <= stops[k + 1].t) { a = stops[k]; b = stops[k + 1]; break; }
    }
    const span = b.t - a.t || 1;
    const f = (t - a.t) / span;
    // Linear in L/C/h - the warm stops sit in a narrow hue band so no hue
    // wrap-around handling is needed.
    const [r, g, bl] = oklchToSrgb(
      a.L + f * (b.L - a.L),
      a.C + f * (b.C - a.C),
      a.h + f * (b.h - a.h),
    );
    lut[i * 3] = r; lut[i * 3 + 1] = g; lut[i * 3 + 2] = bl;
  }
  cached = lut;
  return lut;
}
