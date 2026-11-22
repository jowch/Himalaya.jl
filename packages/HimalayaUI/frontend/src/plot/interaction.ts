export const PEAK_HIT_PX = 10;
export const MIN_SPAN_FRAC = 1e-4;
export interface HitPeak {
  id: number;
  q: number;
}
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