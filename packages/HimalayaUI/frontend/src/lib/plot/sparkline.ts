import type { Trace } from "../../api";

export const SPARK_W = 76;
export const SPARK_H = 28;
const PAD_X = 4;
const PAD_B = 4;
const AMP = 17;

/**
 * Hand-rolled mini-trace path for the scoping worksheet sparkline
 * (series-scoping.html `sparkline()`): log-x, linear-y, peak-relative scale.
 * Pure — returns just the SVG path `d` (the React wrapper supplies stroke/fill
 * and the phase colour). Returns null when there is nothing to draw.
 */
export function sparklinePath(trace: Trace): string | null {
  const pts: Array<[number, number]> = [];
  for (let i = 0; i < trace.q.length; i++) {
    const q = trace.q[i], I = trace.I[i];
    if (Number.isFinite(q) && q > 0 && Number.isFinite(I)) pts.push([q, Math.max(0, I)]);
  }
  if (pts.length < 2) return null;
  const qMin = pts[0][0], qMax = pts[pts.length - 1][0];
  if (qMax <= qMin) return null;
  const lnMin = Math.log(qMin), lnSpan = Math.log(qMax) - lnMin;
  const xOf = (q: number) => PAD_X + ((Math.log(q) - lnMin) / lnSpan) * (SPARK_W - 2 * PAD_X);
  const peak = Math.max(...pts.map((p) => p[1]), 1e-6);
  const hS = AMP / peak;
  const y0 = SPARK_H - PAD_B;
  let d = "";
  pts.forEach(([q, I], i) => {
    d += (i ? "L" : "M") + xOf(q).toFixed(1) + " " + (y0 - I * hS).toFixed(1) + " ";
  });
  return d.trim();
}
