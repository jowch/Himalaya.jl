// Normalized q -> ring-radius map for the focus-workspace detector overlay.
//
// There is NO detector geometry in the data model (no beam center, pixel size,
// sample distance, or wavelength reach the client — the detector image is a
// pre-rendered PNG). So the rings are PRESENTATIONAL: radius proportional to q
// over the peak q-range, exactly like the mockup `focus-workspace.html`'s
// `renderDetector`. They convey "this peak is this far out", not a calibrated
// projection. (#180 open decision 1.)

/** SVG viewBox edge — the overlay draws in a 100x100 space, centered (50,50). */
export const RING_VIEWBOX = 100;
/** Innermost ring radius (clears the beamstop), mockup parity. */
export const RING_R_MIN = 12;
/** Radial span from qLo..qHi, mockup parity (12..45). */
export const RING_R_SPAN = 33;

/** Map a q-value to a normalized ring radius in the 100x100 viewBox. */
export function qToRadius(q: number, qLo: number, qHi: number): number {
  if (qHi <= qLo) return RING_R_MIN;
  const t = (q - qLo) / (qHi - qLo);
  return RING_R_MIN + t * RING_R_SPAN;
}

/** The peak q nearest `hoveredQ` within `tol`, else undefined (match-highlight). */
export function nearestRingQ(
  hoveredQ: number | undefined,
  peakQs: number[],
  tol: number,
): number | undefined {
  if (hoveredQ === undefined) return undefined;
  let best: number | undefined;
  let bestD = Infinity;
  for (const q of peakQs) {
    const d = Math.abs(q - hoveredQ);
    if (d <= tol && d < bestD) { best = q; bestD = d; }
  }
  return best;
}
