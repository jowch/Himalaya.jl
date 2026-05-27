// Shared detector orient-decision — used by DetectorImage (the image canvas)
// and DetectorRingOverlay (the q-link rings), so the ring overlay rotates in
// lockstep with the image instead of re-deriving the rule (#180).

/** Rotate when the container is meaningfully wider than the image's aspect. */
export const ROTATE_THRESHOLD = 1.25;
/** Auto-rotate is gated to viewports >= this width (the 1400px three-column
 *  breakpoint that the retired WorkspaceGrid used; kept as the shared cutoff). */
export const ROTATE_MIN_VIEWPORT = 1400;

export interface OrientInput {
  containerW: number;
  containerH: number;
  imageW: number;
  imageH: number;
  viewportW: number;
}
export interface OrientResult {
  orient: "portrait" | "landscape";
  /** Pre-rotation caps (px). Null in portrait. */
  caps: { maxW: number; maxH: number } | null;
}

export function decideOrient(i: OrientInput): OrientResult {
  const { containerW, containerH, imageW, imageH, viewportW } = i;
  if (!containerW || !containerH || !imageW || !imageH) {
    return { orient: "portrait", caps: null };
  }
  if (viewportW < ROTATE_MIN_VIEWPORT) {
    return { orient: "portrait", caps: null };
  }
  const containerAspect = containerW / containerH;
  const imageAspect = imageW / imageH;
  if (containerAspect > imageAspect * ROTATE_THRESHOLD) {
    // Pre-rotation max-width is capped by container HEIGHT (becomes visual
    // height after the 90° rotation), and max-height by container WIDTH.
    return { orient: "landscape", caps: { maxW: containerH, maxH: containerW } };
  }
  return { orient: "portrait", caps: null };
}
