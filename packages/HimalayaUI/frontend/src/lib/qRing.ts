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

