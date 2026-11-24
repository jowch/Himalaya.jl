import type { Trace, TraceMonitors } from "../api";

export const MAX_COMPARES = 5;
export const COMPARES_FULL = `${MAX_COMPARES} comparisons at most. Remove one to compare another.`;

export function interpolateAt(trace: Pick<Trace, "q" | "I">, q: number): number {
  const xs = trace.q;
  const n = xs.length;
  if (n === 0 || !(q >= xs[0]!) || !(q <= xs[n - 1]!)) return Number.NaN;
  let lo = 0;
  let hi = n - 1;
  while (hi - lo > 1) {
    const mid = (lo + hi) >> 1;
    if (xs[mid]! <= q) lo = mid;
    else hi = mid;
  }
  const x0 = xs[lo]!;
  const x1 = xs[hi]!;
  const y0 = trace.I[lo]!;
  const y1 = trace.I[hi]!;
  return x1 === x0 ? y0 : y0 + ((y1 - y0) * (q - x0)) / (x1 - x0);
}

// physics: without beam monitors the comparison is scaled to the current trace by
// physics: the median intensity ratio over the top 15% of their shared q range, on
// physics: the assumption that by high q both curves are mostly solvent and capillary
// physics: background. Bragg reflections often still sit in that window (HII and
// physics: higher Lα orders), so q within PEAK_MASK of the current trace's own peaks
// physics: is left out. A sample's own smooth high-q scattering (bilayer form factor,
// physics: more lipid) still reads as extra background and over-scales the reference,
// physics: so this is only the fallback for exposures that lack monitors.
export const SCALE_WINDOW = 0.15;
export const PEAK_MASK = 0.008;

export function matchScale(
  main: Pick<Trace, "q" | "I">,
  ref: Pick<Trace, "q" | "I">,
  peakQs: readonly number[] = [],
): number {
  if (main.q.length < 2 || ref.q.length < 2) return 1;
  const lo = Math.max(main.q[0]!, ref.q[0]!);
  const hi = Math.min(main.q[main.q.length - 1]!, ref.q[ref.q.length - 1]!);
  if (!(hi > lo)) return 1;
  const from = hi - (hi - lo) * SCALE_WINDOW;
  const ratios: number[] = [];
  main.q.forEach((q, i) => {
    if (q < from || q > hi) return;
    if (peakQs.some((p) => Math.abs(q - p) < PEAK_MASK)) return;
    const r = interpolateAt(ref, q);
    const m = main.I[i]!;
    if (r > 0 && m > 0) ratios.push(m / r);
  });
  if (ratios.length === 0) return 1;
  ratios.sort((a, b) => a - b);
  return ratios[Math.floor(ratios.length / 2)]!;
}

// physics: SASTool usually divides each trace by its own I_2, the beamstop photodiode
// physics: behind the sample (the setup's integ.mpp has `-i yes i2`). I_2 is the beam
// physics: that got through the sample, so a normalised sample and its buffer are
// physics: already on one scale and the buffer scale is 1. A setup integrated without
// physics: `-i` writes raw counts, which scale with the flux that crossed the sample
// physics: while the detector counted: I_2 × exposure time / counting time, because
// physics: the counter runs longer than the exposure. There the ratio of the two
// physics: fluxes is the scale. Any other or unknown normalisation falls back to the
// physics: high-q match, as does a beam-off image, whose I_2 is only dark noise.
export type BeamMonitors = TraceMonitors;

const positive = (v: number | null | undefined): v is number => typeof v === "number" && Number.isFinite(v) && v > 0;

export function transmittedFlux(m: BeamMonitors | undefined): number | null {
  if (!m || !positive(m.transmitted_counts) || !positive(m.exposure_time) || !positive(m.counting_time)) return null;
  return (m.transmitted_counts * m.exposure_time) / m.counting_time;
}

export interface CompareScale {
  scale: number;
  source: "normalised" | "monitor" | "matched" | "other_setup" | "beam_off" | "unknown_normalisation";
}

export const isMatchedScale = (s: CompareScale): boolean => s.source !== "normalised" && s.source !== "monitor";

type ScaledTrace = Pick<Trace, "q" | "I" | "monitors">;

// physics: counts per pixel also depend on the sample-to-detector distance and on
// physics: where the beamstop diode sits, which change between setups. A flux ratio
// physics: puts two images on one scale only within a setup, so across setups the
// physics: high-q match is the better estimate.
export function compareScale(main: ScaledTrace, ref: ScaledTrace, peakQs: readonly number[], sameSetup: boolean): CompareScale {
  const matched = (source: CompareScale["source"]): CompareScale => ({ scale: matchScale(main, ref, peakQs), source });
  if (!sameSetup) return matched("other_setup");
  if (main.monitors?.beam_off || ref.monitors?.beam_off) return matched("beam_off");
  const by = main.monitors?.normalised_by;
  if (by === undefined || by !== ref.monitors?.normalised_by) return matched("unknown_normalisation");
  if (by === "i2") return { scale: 1, source: "normalised" };
  if (by !== "none") return matched("unknown_normalisation");
  const fm = transmittedFlux(main.monitors);
  const fr = transmittedFlux(ref.monitors);
  return fm !== null && fr !== null ? { scale: fm / fr, source: "monitor" } : matched("matched");
}

export function scaledTrace(ref: Pick<Trace, "q" | "I">, scale: number): { q: number[]; I: number[] } {
  return { q: [...ref.q], I: ref.I.map((v) => v * scale) };
}

export function subtractTrace(
  main: Pick<Trace, "q" | "I">,
  ref: Pick<Trace, "q" | "I">,
  scale: number,
): { q: number[]; d: number[] } {
  const q: number[] = [];
  const d: number[] = [];
  main.q.forEach((x, i) => {
    const r = interpolateAt(ref, x);
    if (Number.isFinite(r)) {
      q.push(x);
      d.push(main.I[i]! - scale * r);
    }
  });
  return { q, d };
}

export function robustRange(values: readonly number[], lowFrac = 0.01, highFrac = 0.99): [number, number] {
  const finite = values.filter(Number.isFinite).sort((a, b) => a - b);
  if (finite.length === 0) return [-1, 1];
  const at = (f: number): number => finite[Math.min(finite.length - 1, Math.max(0, Math.round(f * (finite.length - 1))))]!;
  const lo = Math.min(0, at(lowFrac));
  const hi = Math.max(0, at(highFrac));
  return hi > lo ? [lo, hi] : [lo - 1, hi + 1];
}

// physics: below the first Bragg peak the difference is dominated by the
// physics: beamstop and low-q upturn, which rarely subtracts cleanly and swings far
// physics: past the peaks. The strip's y range therefore starts just below the
// physics: first peak in view; with no peak in view it uses the upper half of the
// physics: visible q span, where that residual has died away.
export const LOW_Q_SKIP = 0.5;

export function stripRange(
  q: readonly number[],
  d: readonly number[],
  domain: [number, number],
  peakQs: readonly number[],
): [number, number] {
  const firstPeak = Math.min(...peakQs.filter((p) => p >= domain[0] && p <= domain[1]));
  const from = Number.isFinite(firstPeak)
    ? Math.max(domain[0], firstPeak - 2 * PEAK_MASK)
    : domain[0] + (domain[1] - domain[0]) * LOW_Q_SKIP;
  const inView = d.filter((_, i) => q[i]! >= from && q[i]! <= domain[1]);
  const [lo, hi] = robustRange(inView);
  const top = Math.max(hi, ...inView.filter(Number.isFinite));
  return [lo, top];
}

export function clampDomain(
  domain: [number, number] | null,
  extent: [number, number] | null,
): [number, number] | null {
  if (!domain || !extent) return null;
  const lo = Math.max(domain[0], extent[0]);
  const hi = Math.min(domain[1], extent[1]);
  return hi > lo ? [lo, hi] : null;
}
