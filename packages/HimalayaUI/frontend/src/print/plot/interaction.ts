/**
 * Pure interaction math ported from the legacy PlotSurface so the greenfield
 * engine reproduces its feel exactly. No DOM, no d3 — caller supplies the q→px
 * projection.
 */

/** Click within this many pixels of a peak to act on it. (PlotSurface #180.) */
export const PEAK_HIT_PX = 10;
/** Drags shorter than this are treated as clicks, not brushes. */
export const BRUSH_DEADZONE_PX = 4;
/** Refuse zooms narrower than this fraction of the full extent. */
export const MIN_SPAN_FRAC = 1e-4;

export interface HitPeak {
  id: number;
  q: number;
}

/**
 * Nearest clickable peak (id ≥ 0) to a pixel x, within `tolerancePx`, or null.
 * On an exact tie the later peak wins (`<=`), matching PlotSurface:104-123.
 */
export function hitTestPeaks<T extends HitPeak>(
  peaks: T[],
  clickX: number,
  toPx: (q: number) => number,
  tolerancePx: number,
): T | null {
  let best: T | null = null;
  let bestDist = tolerancePx;
  for (const p of peaks) {
    if (p.id < 0) continue;
    const px = toPx(p.q);
    if (!Number.isFinite(px)) continue;
    const d = Math.abs(px - clickX);
    if (d <= bestDist) {
      best = p;
      bestDist = d;
    }
  }
  return best;
}

/**
 * Wheel-zoom the x domain about the cursor (log-aware). Returns the new
 * [min, max] clamped to `extent`, or null when the zoom would collapse the
 * span below MIN_SPAN_FRAC of the extent. Ported from PlotSurface:258-284.
 */
export function zoomXDomain(opts: {
  cursorQ: number;
  deltaY: number;
  current: [number, number];
  extent: [number, number];
  type: "log" | "linear";
}): [number, number] | null {
  const { cursorQ, deltaY, current, extent, type } = opts;
  const [curMin, curMax] = current;
  const [q0, qN] = extent;
  const factor = Math.exp(deltaY * 0.001);
  let newMin: number;
  let newMax: number;
  if (type === "log") {
    const logMin = Math.log(curMin);
    const logMax = Math.log(curMax);
    const logCur = Math.log(Math.max(cursorQ, 1e-6));
    newMin = Math.max(q0, Math.exp(logCur - (logCur - logMin) * factor));
    newMax = Math.min(qN, Math.exp(logCur + (logMax - logCur) * factor));
  } else {
    newMin = Math.max(q0, cursorQ - (cursorQ - curMin) * factor);
    newMax = Math.min(qN, cursorQ + (curMax - cursorQ) * factor);
  }
  if (newMax - newMin < (qN - q0) * MIN_SPAN_FRAC) return null;
  return [newMin, newMax];
}
