import type { Setup } from "../api";

export interface DetectorGeometry {
  beamX: number;
  beamY: number;
  pixelM: number;
  wavelengthA: number;
  distanceM: number;
  perAngstrom: number;
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

const positive = (v: number | null | undefined): v is number => typeof v === "number" && Number.isFinite(v) && v > 0;

export function detectorGeometry(
  setup: Pick<Setup, "beam_center_x" | "beam_center_y" | "pixel_size_um" | "energy_kev" | "flight_path_m" | "q_units"> | undefined,
): DetectorGeometry | null {
  if (!setup) return null;
  const { beam_center_x: bx, beam_center_y: by, pixel_size_um: px, energy_kev: e, flight_path_m: l } = setup;
  if (!positive(bx) || !positive(by) || !positive(px) || !positive(e) || !positive(l)) return null;
  return {
    beamX: bx,
    beamY: by,
    pixelM: px * 1e-6,
    wavelengthA: HC_KEV_ANGSTROM / e,
    distanceM: l,
    perAngstrom: /nm/i.test(setup.q_units ?? "") ? 0.1 : 1,
  };
}

// physics: Bragg's law, q = 4π sin θ / λ; a ring at scattering angle 2θ lands L tan 2θ from the beam on a flat detector.
function radiusPx(q: number, g: DetectorGeometry): number | null {
  const s = (q * g.perAngstrom * g.wavelengthA) / (4 * Math.PI);
  if (!(s > 0) || s >= Math.SQRT1_2) return null;
  return (g.distanceM * Math.tan(2 * Math.asin(s))) / g.pixelM;
}

// physics: setup_info gives the beam centre in detector pixels from the bottom-left corner with y up; the image's rows run top-down.
export function ringCircles(qs: readonly number[], g: DetectorGeometry, size: ImageSize | undefined): RingCircles | null {
  if (!size || !(size.width > 0) || !(size.height > 0)) return null;
  const radii = qs.map((q) => radiusPx(q, g)).filter((r): r is number => r !== null);
  return { cx: g.beamX + 0.5, cy: size.height - g.beamY - 0.5, radii };
}
