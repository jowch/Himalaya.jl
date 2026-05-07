/**
 * Adaptive normalization library for the Compare page multi-trace plot
 * (Plan §Phase 6, Task 6.1; spec §Plot rendering / Working band vs total band).
 *
 * Produces the per-member reference value used to scale a trace into its
 * y-band, and the y-mapped points (with tail clipping at the total band
 * envelope so a tall low-q dropoff doesn't bleed into the next member's
 * working band).
 *
 * `yBand` is treated as a numeric *envelope* `[topPx, bottomPx]` in plot
 * coordinates where `topPx < bottomPx` (small-y is up — matches the SVG /
 * Observable Plot convention). The reference value maps to the TOP of the
 * inner working band (default 70%); zero intensity maps to the bottom of
 * the working band; over/underflow clips at the total envelope edges.
 *
 * Why peak-fit-with-signal-fallback for `qwindow`:
 * - Peak-bearing traces (typical exposures) → reference = max peak intensity
 *   in the window, so peak-position alignment reads cleanly across the stack.
 * - Peakless / all-manual-peak traces (buffers, controls) → fall back to
 *   the in-window signal max so the trace's overall envelope still fills
 *   the band rather than collapsing flat. Manual peaks have intensity:null
 *   (server contract: see api.ts Peak.intensity) and are excluded from the
 *   peak-fit max — if every peak in the window is manual, we go to signal.
 */

export type Normalization = "none" | "max" | "area" | "qwindow";
export type QWindow = [number, number] | null;

/**
 * Subset of `MemberSnapshotPeak` consumed by normalization.
 *
 * Kept structural (not importing the full snapshot type) so this library
 * stays self-contained and unit-testable without TanStack/router context.
 */
export interface Peak {
  q: number;
  intensity: number | null;
}

interface TraceLike {
  q: number[];
  I: number[];
}

/**
 * Default fraction of the total band height occupied by the working band.
 * Outer (1 - default) is split evenly as headroom / footroom; see spec
 * §Working band vs total band for the SAXS rationale.
 */
export const DEFAULT_WORKING_BAND_FRACTION = 0.7;

/** Minimum positive reference returned to callers — guards division-by-zero. */
const MIN_REFERENCE = 1e-12;

// ── computeReference ─────────────────────────────────────────────────────────

export function computeReference(
  trace: TraceLike,
  peaks: Peak[],
  qWindow: QWindow,
  normalization: Normalization,
): number {
  switch (normalization) {
    case "none":
      return 1;
    case "max":
      return positiveOr(signalMaxInWindow(trace, qWindow), MIN_REFERENCE);
    case "area":
      return positiveOr(trapezoidalArea(trace, qWindow), MIN_REFERENCE);
    case "qwindow":
      return adaptiveQWindowReference(trace, peaks, qWindow);
  }
}

/**
 * Adaptive peak-fit/signal-fit reference for the `qwindow` strategy.
 *
 * Inclusive boundaries: `q_window_min ≤ q ≤ q_window_max` (matches the spec
 * §Default q_window_min/max and the plot rendering rules).
 */
function adaptiveQWindowReference(
  trace: TraceLike,
  peaks: Peak[],
  qWindow: QWindow,
): number {
  const peaksInWindow = peaks.filter((p) => inWindow(p.q, qWindow));
  // Peak-fit only counts peaks with non-null intensity (manual peaks
  // contribute nothing — their intensity is null on the server contract).
  let peakMax = 0;
  let havePeakIntensity = false;
  for (const p of peaksInWindow) {
    if (typeof p.intensity === "number" && Number.isFinite(p.intensity)) {
      havePeakIntensity = true;
      if (p.intensity > peakMax) peakMax = p.intensity;
    }
  }
  if (havePeakIntensity && peakMax > 0) return peakMax;
  // Fallback: in-window signal max.
  return positiveOr(signalMaxInWindow(trace, qWindow), MIN_REFERENCE);
}

// ── applyNormalization ───────────────────────────────────────────────────────

export function applyNormalization(
  trace: TraceLike,
  reference: number,
  yBand: [number, number],
  workingBandFraction: number = DEFAULT_WORKING_BAND_FRACTION,
): Array<{ q: number; y: number }> {
  if (trace.q.length === 0) return [];

  // Total band envelope: [top, bottom] in plot coords (small y = up).
  const [bandTop, bandBottom] = yBand[0] <= yBand[1] ? yBand : [yBand[1], yBand[0]];
  const totalHeight = bandBottom - bandTop;

  // Working band: centered inside the total band, occupying
  // `workingBandFraction` of its height. Padding split evenly above/below.
  const padding = (totalHeight * (1 - workingBandFraction)) / 2;
  const workTop    = bandTop + padding;
  const workBottom = bandBottom - padding;
  const workHeight = workBottom - workTop;

  // Guard against division-by-zero from a degenerate reference. If the
  // reference is non-positive (e.g. blank trace), collapse the trace to
  // the working band's bottom — every point lands at the floor.
  const safeRef = reference > 0 ? reference : MIN_REFERENCE;

  const out: Array<{ q: number; y: number }> = new Array(trace.q.length);
  for (let i = 0; i < trace.q.length; i++) {
    const q = trace.q[i]!;
    const I = trace.I[i]!;
    // Map I to a fraction of the working band: 0 → bottom, reference → top.
    const frac = Number.isFinite(I) ? I / safeRef : 0;
    let y = workBottom - frac * workHeight;
    // Clip at total band envelope (so a low-q dropoff doesn't bleed into
    // the adjacent member's working band).
    if (y < bandTop)    y = bandTop;
    if (y > bandBottom) y = bandBottom;
    out[i] = { q, y };
  }
  return out;
}

// ── helpers ─────────────────────────────────────────────────────────────────

function inWindow(q: number, w: QWindow): boolean {
  if (w === null) return true;
  return q >= w[0] && q <= w[1];
}

function signalMaxInWindow(trace: TraceLike, qWindow: QWindow): number {
  let max = 0;
  for (let i = 0; i < trace.q.length; i++) {
    const q = trace.q[i]!;
    if (!inWindow(q, qWindow)) continue;
    const v = trace.I[i]!;
    if (Number.isFinite(v) && v > max) max = v;
  }
  return max;
}

function trapezoidalArea(trace: TraceLike, qWindow: QWindow): number {
  // Integrate I dq over the (clamped) q-window using the trapezoidal rule.
  // Falls back to 0 when fewer than two in-window points exist.
  const qs: number[] = [];
  const Is: number[] = [];
  for (let i = 0; i < trace.q.length; i++) {
    const q = trace.q[i]!;
    if (!inWindow(q, qWindow)) continue;
    qs.push(q);
    Is.push(trace.I[i]!);
  }
  if (qs.length < 2) return 0;
  let area = 0;
  for (let i = 1; i < qs.length; i++) {
    const dq = qs[i]! - qs[i - 1]!;
    const avg = (Is[i]! + Is[i - 1]!) / 2;
    if (Number.isFinite(dq) && Number.isFinite(avg)) area += dq * avg;
  }
  return area;
}

function positiveOr(v: number, fallback: number): number {
  return Number.isFinite(v) && v > 0 ? v : fallback;
}
