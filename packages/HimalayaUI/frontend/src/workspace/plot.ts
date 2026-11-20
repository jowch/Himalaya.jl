import type { Peak, Trace } from "../api";

export interface CombLine {
  id: number;
  qs: readonly number[];
  color: string;
  emphasis: boolean;
}

export interface Geometry {
  plot: { left: number; right: number; top: number; bottom: number };
  path: [number, number][];
  peaks: { id: number; x: number; y: number; excluded: boolean }[];
  combs: { id: number; xs: number[]; color: string; emphasis: boolean }[];
  xTicks: { x: number; label: string }[];
}

const PAD = { left: 12, right: 12, top: 16, bottom: 40 };

function niceStep(span: number, target: number): number {
  const raw = span / target;
  const mag = 10 ** Math.floor(Math.log10(raw));
  const f = raw / mag;
  return (f < 1.5 ? 1 : f < 3.5 ? 2 : f < 7.5 ? 5 : 10) * mag;
}

export function traceGeometry(
  trace: Trace | undefined,
  peaks: readonly Peak[],
  combs: readonly CombLine[],
  width: number,
  height: number,
): Geometry | null {
  if (!trace) return null;
  const pts: [number, number][] = [];
  for (let i = 0; i < trace.q.length; i++) {
    const q = trace.q[i]!;
    const I = trace.I[i]!;
    if (Number.isFinite(q) && Number.isFinite(I) && I > 0) pts.push([q, Math.log10(I)]);
  }
  if (pts.length < 2) return null;
  const q0 = pts[0]![0];
  const q1 = pts[pts.length - 1]![0];
  let y0 = Infinity;
  let y1 = -Infinity;
  for (const [, y] of pts) {
    if (y < y0) y0 = y;
    if (y > y1) y1 = y;
  }
  if (y1 === y0) y1 = y0 + 1;
  const plot = { left: PAD.left, right: width - PAD.right, top: PAD.top, bottom: height - PAD.bottom };
  const sx = (q: number) => plot.left + ((q - q0) / (q1 - q0)) * (plot.right - plot.left);
  const sy = (ly: number) => plot.bottom - ((ly - y0) / (y1 - y0)) * (plot.bottom - plot.top);
  const path = pts.map(([q, ly]) => [sx(q), sy(ly)] as [number, number]);

  const yAt = (q: number): number => {
    let lo = 0;
    let hi = pts.length - 1;
    while (hi - lo > 1) {
      const mid = (lo + hi) >> 1;
      if (pts[mid]![0] < q) lo = mid;
      else hi = mid;
    }
    const a = pts[lo]!;
    const b = pts[hi]!;
    const t = b[0] === a[0] ? 0 : (q - a[0]) / (b[0] - a[0]);
    return sy(a[1] + t * (b[1] - a[1]));
  };

  const inRange = (q: number) => q >= q0 && q <= q1;
  const step = niceStep(q1 - q0, Math.max(2, Math.floor((plot.right - plot.left) / 90)));
  const xTicks: Geometry["xTicks"] = [];
  const digits = Math.max(0, -Math.floor(Math.log10(step)));
  for (let t = Math.ceil(q0 / step) * step; t <= q1 + 1e-12; t += step) {
    xTicks.push({ x: sx(t), label: t.toFixed(digits) });
  }

  return {
    plot,
    path,
    peaks: peaks.filter((p) => inRange(p.q)).map((p) => ({ id: p.id, x: sx(p.q), y: yAt(p.q), excluded: p.excluded })),
    combs: combs.map((c) => ({ id: c.id, color: c.color, emphasis: c.emphasis, xs: c.qs.filter(inRange).map(sx) })),
    xTicks,
  };
}
