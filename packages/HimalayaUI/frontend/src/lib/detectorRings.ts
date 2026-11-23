import type { Setup } from "../api";

export interface DetectorGeometry {
  beamX: number;
  beamY: number;
  pixelM: number;
  wavelengthA: number;
  distanceM: number;
}

export interface ImageSize {
  width: number;
  height: number;
}

export interface RingCircles {
  cx: number;
  cy: number;
  radii: number[];
}

const HC_KEV_ANGSTROM = 12.39842;

const finite = (v: number | null | undefined): v is number => typeof v === "number" && Number.isFinite(v);
const positive = (v: number | null | undefined): v is number => finite(v) && v > 0;

export const perAngstrom = (qUnits: string | null | undefined): number => (/nm/i.test(qUnits ?? "") ? 0.1 : 1);

export function detectorGeometry(
  setup: Pick<Setup, "beam_center_x" | "beam_center_y" | "pixel_size_um" | "energy_kev" | "flight_path_m"> | undefined,
): DetectorGeometry | null {
  if (!setup) return null;
  const { beam_center_x: bx, beam_center_y: by, pixel_size_um: px, energy_kev: e, flight_path_m: l } = setup;
  if (!finite(bx) || !finite(by) || !positive(px) || !positive(e) || !positive(l)) return null;
  return {
    beamX: bx,
    beamY: by,
    pixelM: px * 1e-6,
    wavelengthA: HC_KEV_ANGSTROM / e,
    distanceM: l,
  };
}

// physics: Bragg's law, q = 4π sin θ / λ; a ring at scattering angle 2θ lands L tan 2θ from the beam on a flat detector.
function radiusPx(qPerAngstrom: number, g: DetectorGeometry): number | null {
  const s = (qPerAngstrom * g.wavelengthA) / (4 * Math.PI);
  if (!(s > 0) || s >= Math.SQRT1_2) return null;
  return (g.distanceM * Math.tan(2 * Math.asin(s))) / g.pixelM;
}

// physics: setup_info gives the beam centre in detector pixels from the bottom-left corner with y up; the image's rows run top-down.
export function ringCircles(
  qs: readonly number[], qUnits: string | null | undefined, g: DetectorGeometry, size: ImageSize | undefined,
): RingCircles | null {
  if (!size || !(size.width > 0) || !(size.height > 0)) return null;
  const toAngstrom = perAngstrom(qUnits);
  const radii = qs.map((q) => radiusPx(q * toAngstrom, g)).filter((r): r is number => r !== null);
  return { cx: g.beamX + 0.5, cy: size.height - g.beamY - 0.5, radii };
}
