/**
 * folioFigure — pure geometry for the series folio card's read-only miniature
 * (R6, finding F-A / F-B; mockup `series-folio.html` `figSVG()` + `.ps`).
 *
 * Two builders, both pure (no DOM, no React) so they unit-test in isolation:
 *
 *  - `buildMiniWaterfall(members)` reproduces the mockup's `figSVG()` geometry
 *    — a small stacked log-x waterfall — but fed by LIVE data: each member's
 *    `snapshot.effective_peaks` (the curation-frozen peak set) supplies the
 *    Bragg comps, and `snapshot.confirmed_index.phase` supplies the colour.
 *    This is the deliberate "read-only variant" the issue calls for: it reuses
 *    the geometry, NOT the heavy interactive `MultiTracePlot` (one Observable-
 *    Plot instance per card across a masonry would be the wrong tool).
 *
 *  - `buildPhaseStrip(members)` produces one cell PER sample coloured by that
 *    member's confirmed phase, plus a transition caption ("Pn3m → Lamellar").
 *
 * The intensity model (decaying background + Gaussian Bragg peaks) is copied
 * byte-for-byte from the mockup so the live miniature reads identically.
 */
import type { SeriesMember } from "../../api";
import { phaseColor } from "../../phases";

// ── Phase strip ─────────────────────────────────────────────────────────────

export interface StripSegment {
  /** Dominant phase for this member, or null = unindexed. */
  phase: string | null;
  /** Second phase for a 2-stop coexistence gradient (null = single phase). */
  coexistWith: string | null;
}

export interface PhaseStripModel {
  segments: StripSegment[];
  /** "transition" (first ≠ last indexed), "throughout" (one phase), "none". */
  kind: "transition" | "throughout" | "none";
  first: string | null;
  last: string | null;
}

/** Dominant phase of a member = its confirmed index phase (null if unindexed). */
export function memberPhase(m: SeriesMember): string | null {
  return m.snapshot?.confirmed_index?.phase ?? null;
}

/**
 * Build the per-sample phase strip + caption.
 *
 * `coexistResolver` is an optional escape hatch for true two-phase coexistence:
 * a member snapshot carries a single `confirmed_index`, so production passes no
 * resolver (every segment is single-phase). The gradient rendering path stays
 * real and tested — and ready for a future multi-index snapshot — by letting a
 * caller (today, only the unit test) supply the second phase.
 */
export function buildPhaseStrip(
  members: SeriesMember[],
  coexistResolver?: (m: SeriesMember) => string | null,
): PhaseStripModel {
  const segments: StripSegment[] = members.map((m) => ({
    phase: memberPhase(m),
    coexistWith: coexistResolver ? coexistResolver(m) : null,
  }));

  const indexed = segments.map((s) => s.phase).filter((p): p is string => p !== null);
  const first = indexed.length > 0 ? indexed[0]! : null;
  const last = indexed.length > 0 ? indexed[indexed.length - 1]! : null;

  let kind: PhaseStripModel["kind"];
  if (first === null || last === null) {
    kind = "none";
  } else {
    // A transition exists when the strip spans more than one distinct phase.
    const distinct = new Set(indexed);
    kind = distinct.size > 1 ? "transition" : "throughout";
  }

  return { segments, kind, first, last };
}

// ── Mini-waterfall ───────────────────────────────────────────────────────────

export interface WaterfallTick {
  x: number;
  /** Tick height in px (taller for the principal peak). */
  h: number;
  color: string;
  opacity: number;
}

export interface WaterfallTrace {
  /** Baseline y for this member's row (top member smallest y). */
  baselineY: number;
  /** SVG path `d` for the integration line. */
  linePath: string;
  /** SVG path `d` for the filled area (line closed down to the baseline). */
  fillPath: string;
  color: string;
  ticks: WaterfallTick[];
}

export interface WaterfallModel {
  width: number;
  height: number;
  traces: WaterfallTrace[];
}

const UNINDEXED_COLOR = "var(--color-ink-faint)";

// Geometry constants (mockup `figSVG()`).
const W = 340;
const STEP_Y = 31;
const AMP = 24;
const PAD_T = 12;
const PAD_B = 12;
const PAD_L = 14;
const PAD_R = 14;
const Q_MIN = 0.02;
const Q_MAX = 0.42;
const N_SAMP = 150;

interface Comp {
  /** Principal peak q (the smallest effective-peak q for the member). */
  q0: number;
  /** Amplitude, normalised to the member's strongest peak. */
  amp: number;
}

/** Decaying background + Gaussian Bragg peaks (mockup `intensity()`). */
function intensity(q: number, comps: Comp[]): number {
  let v = 0.26 * Math.pow(0.02 / q, 0.9);
  for (const c of comps) {
    // The mockup expands a ratio series; live data already lists the peaks, so
    // each effective peak is its own Gaussian at its measured q.
    const qp = c.q0;
    const h = c.amp;
    const sg = 0.03 * qp + 0.0006;
    v += h * Math.exp(-0.5 * Math.pow((q - qp) / sg, 2));
  }
  return v;
}

/**
 * Build the stacked mini-waterfall from live member snapshots. Each member's
 * `effective_peaks` become Gaussian comps (amplitude normalised to the
 * member's strongest peak); the line colour is the confirmed phase colour.
 */
export function buildMiniWaterfall(members: SeriesMember[]): WaterfallModel {
  const n = members.length;
  const height = PAD_T + PAD_B + Math.max(n, 1) * STEP_Y;
  const plotW = W - PAD_L - PAD_R;

  const logMin = Math.log(Q_MIN);
  const logMax = Math.log(Q_MAX);
  const xOf = (q: number): number =>
    PAD_L + ((Math.log(q) - logMin) / (logMax - logMin)) * plotW;

  // Sample q-grid in log space (shared across traces).
  const qs: number[] = [];
  for (let i = 0; i <= N_SAMP; i++) {
    const f = i / N_SAMP;
    qs.push(Math.exp(logMin + f * (logMax - logMin)));
  }

  // Per-member comps, plus a global max for amplitude normalisation.
  const compsByMember: Comp[][] = members.map((m) => {
    const peaks = m.snapshot?.effective_peaks ?? [];
    if (peaks.length === 0) return [];
    const intensities = peaks.map((p) => (p.intensity ?? 0) || 0);
    const peakMax = Math.max(1e-9, ...intensities);
    return peaks.map((p) => ({
      q0: p.q,
      amp: ((p.intensity ?? 0) || 0.5) / peakMax,
    }));
  });

  let gmax = 1e-6;
  const curves = members.map((_, i) => {
    const comps = compsByMember[i]!;
    const ys = qs.map((q) => intensity(q, comps));
    gmax = Math.max(gmax, ...ys);
    return ys;
  });
  const hScale = AMP / gmax;
  // Stack top→bottom in sample order: member 0 sits at the top (smallest y),
  // each later member one step lower. This reads in the same direction as the
  // per-sample phase strip below it, so the card scans as one figure. (The
  // mockup stacked bottom-up; we flip it for strip/figure parity.) The curve
  // bulges UPWARD from its baseline, so the baseline leaves a STEP_Y headroom
  // above it for the peaks.
  const baseY = (i: number): number => PAD_T + STEP_Y + i * STEP_Y;

  const traces: WaterfallTrace[] = members.map((m, i) => {
    const y0 = baseY(i);
    const phase = memberPhase(m);
    const color = phase ? phaseColor(phase) : UNINDEXED_COLOR;

    let d = "";
    qs.forEach((q, k) => {
      const x = xOf(q);
      const y = y0 - curves[i]![k]! * hScale;
      d += (k === 0 ? "M" : "L") + x.toFixed(1) + " " + y.toFixed(1) + " ";
    });
    const linePath = d.trimEnd();
    const fillPath = `${d}L${(W - PAD_R).toFixed(1)} ${y0.toFixed(1)} L${PAD_L.toFixed(1)} ${y0.toFixed(1)} Z`;

    // Up to three short peak ticks hugging the baseline (principal peaks first,
    // by q ascending — the smallest few q's carry the indexing signal).
    const peaks = (m.snapshot?.effective_peaks ?? [])
      .slice()
      .sort((a, b) => a.q - b.q)
      .slice(0, 3);
    const ticks: WaterfallTick[] = [];
    peaks.forEach((p, k) => {
      if (p.q < Q_MIN || p.q > Q_MAX) return;
      ticks.push({
        x: xOf(p.q),
        h: 4 + 3 * Math.pow(0.6, k),
        color,
        opacity: Number((0.5 + 0.35 * Math.pow(0.7, k)).toFixed(2)),
      });
    });

    return { baselineY: y0, linePath, fillPath, color, ticks };
  });

  return { width: W, height, traces };
}
