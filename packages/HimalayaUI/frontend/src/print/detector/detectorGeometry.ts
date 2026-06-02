import { RING_R_MIN, RING_R_SPAN, RING_VIEWBOX } from "../../lib/qRing";

/** Standard SAXS/WAXS calibration (PONI-style). All lengths share one unit
 *  system; beamCenter + imageSize are in the displayed image's pixel space. */
export interface DetectorCalibration {
  beamCenterPx: { x: number; y: number };
  imageSizePx: { w: number; h: number };
  sampleDistanceMm: number;
  pixelSizeMm: number;
  energyKeV: number;
}

/** A ring to draw, in NORMALIZED image coordinates: radius is a fraction of the
 *  image width, so the component is resolution-independent. */
export interface RingPlacement {
  q: number;
  r: number;
  color?: string;
  ghost?: boolean;
}

/** A q-value plus optional render hints. A bare number is sugar for { q }. */
export type RingInput = number | { q: number; color?: string; ghost?: boolean };

const PLANCK_KEV_ANGSTROM = 12.39842; // hc in keV*Angstrom

/** q (1/Angstrom) -> ring radius as a fraction of image width, via SAXS geometry. */
export function qToImageRadius(q: number, cal: DetectorCalibration): number {
  const lambda = PLANCK_KEV_ANGSTROM / cal.energyKeV;          // Angstrom
  const sinTheta = Math.min(1, (q * lambda) / (4 * Math.PI));
  const theta = Math.asin(sinTheta);
  const rMm = cal.sampleDistanceMm * Math.tan(2 * theta);
  const rPx = rMm / cal.pixelSizeMm;
  return rPx / cal.imageSizePx.w;
}

function norm(input: RingInput): { q: number; color?: string; ghost?: boolean } {
  return typeof input === "number" ? { q: input } : input;
}

/** Build ring placements (normalized image coords) + the beam center.
 *  `null` calibration -> presentational fallback (centered, radius prop to q-range),
 *  matching the legacy behaviour until real geometry data is wired. */
export function buildRingPlacements(
  inputs: RingInput[],
  cal: DetectorCalibration | null,
): { beamCenter: { x: number; y: number }; rings: RingPlacement[] } {
  const items = inputs.map(norm);
  if (cal) {
    return {
      beamCenter: {
        x: cal.beamCenterPx.x / cal.imageSizePx.w,
        y: cal.beamCenterPx.y / cal.imageSizePx.h,
      },
      rings: items.map((it) => ({
        q: it.q,
        r: qToImageRadius(it.q, cal),
        ...(it.color !== undefined ? { color: it.color } : {}),
        ...(it.ghost !== undefined ? { ghost: it.ghost } : {}),
      })),
    };
  }
  // Presentational fallback: normalize the old 100-unit viewBox radii to 0..1.
  const qs = items.map((it) => it.q);
  const qLo = qs.length ? Math.min(...qs) : 0;
  const qHi = qs.length ? Math.max(...qs) : 1;
  const span = qHi - qLo || 1;
  return {
    beamCenter: { x: 0.5, y: 0.5 },
    rings: items.map((it) => {
      const t = (it.q - qLo) / span;
      const rViewbox = RING_R_MIN + t * RING_R_SPAN; // 12..45 in the 100 box
      return {
        q: it.q,
        r: rViewbox / RING_VIEWBOX, // -> 0.12..0.45 normalized
        ...(it.color !== undefined ? { color: it.color } : {}),
        ...(it.ghost !== undefined ? { ghost: it.ghost } : {}),
      };
    }),
  };
}
