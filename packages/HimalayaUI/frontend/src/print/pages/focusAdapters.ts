// focusAdapters.ts (Phase-4 Focus) — PURE (no React) functions that turn carried
// backend data (api.ts `Peak` / `IndexEntry` / `Trace`) into the prop shapes the
// greenfield "The Print" composites expect (TraceModel / RingInput / CombSeries).
//
// The custom-index physics is NOT re-implemented here — it is reused wholesale
// from `src/lib/customIndex.ts` (SYMS / customRefls / landsOn /
// latticeForFirstOrderOnPeak). The legacy React components this ports the
// presentation-free logic out of are cited inline.

import type { Peak, IndexEntry, Trace } from "../../api";
import { phaseColor } from "../../phases";
import type { TraceModel } from "../plot/TracePlot";
import type { PlotPeak } from "../plot/marks/PlotPeaks";
import type { RingInput } from "../detector/detectorGeometry";
import type { CombSeries, CombTooth } from "../comb/combModel";
import {
  SYMS,
  customRefls,
  landsOn,
  latticeForFirstOrderOnPeak,
} from "../../lib/customIndex";

// The q-link / claim tolerance the legacy surfaces use (CombPanel,
// FocusDetectorPanel): span-relative, floored at 1e-6.
function spanTol(qs: number[]): number {
  const lo = qs.length ? Math.min(...qs) : 0;
  const hi = qs.length ? Math.max(...qs) : 1;
  return Math.max((hi - lo) * 0.02, 1e-6);
}

// ── trace peaks ───────────────────────────────────────────────────────────────

/**
 * Observed `Peak[]` → trace `PlotPeak[]`; `phase` drives the trace colour
 * downstream (TracePlot resolves it via phaseColor).
 *
 * Minimal field map ({ id, q, intensity?, source, excluded }). The legacy
 * PlotCard → TraceViewer mapping does NOT decorate trace peaks with
 * `predictedAbsent`/`color` — those are index-overlay concerns (the comb /
 * detector rings carry the assignment colouring), so this adapter takes no
 * `activeIndices` param. `intensity` is conditionally spread (never set to
 * `undefined`) to satisfy `exactOptionalPropertyTypes`.
 */
export function toTraceModel(
  trace: Trace,
  peaks: Peak[],
  phase: string | null,
): TraceModel {
  return {
    trace,
    phase,
    peaks: peaks.map((p): PlotPeak => ({
      id: p.id,
      q: p.q,
      source: p.source,
      excluded: p.excluded,
      ...(p.intensity != null ? { intensity: p.intensity } : {}),
    })),
  };
}

// ── losing / complement peak sets (PlotCard.tsx 250-263, verbatim) ───────────

/**
 * The peaks an active phase would orphan if the hovered candidate were added.
 * Ported VERBATIM from PlotCard's `losingPeakIds` useMemo: for each active
 * index that shares ≥1 claimed peak with the hovered candidate, every claimed
 * peak it does NOT share is "losing". Empty set when there is no hover, the
 * hovered candidate is already active, or it overlaps no active phase
 * (independent phases coexist and lose nothing).
 */
export function losingPeakIds(
  hovered: IndexEntry | undefined,
  activeIndices: IndexEntry[],
): Set<number> {
  const losing = new Set<number>();
  if (!hovered) return losing;
  const alreadyActive = activeIndices.some((ix) => ix.id === hovered.id);
  if (alreadyActive) return losing;
  const claim = new Set(hovered.peaks.map((p) => p.peak_id));
  for (const active of activeIndices) {
    const activePeakIds = active.peaks.map((p) => p.peak_id);
    const overlaps = activePeakIds.some((id) => claim.has(id));
    if (!overlaps) continue; // independent phase → coexists, nothing lost
    for (const id of activePeakIds) if (!claim.has(id)) losing.add(id);
  }
  return losing;
}

/**
 * The complement of `losing` within `allPeakIds`. TracePlot's
 * `highlightPeakIds` KEEPS the peaks in its set (dimming the rest); to dim the
 * LOSING peaks we therefore pass everything EXCEPT the losing set.
 */
export function complementPeakIds(
  allPeakIds: Iterable<number>,
  losing: Set<number>,
): Set<number> {
  const out = new Set<number>();
  for (const id of allPeakIds) if (!losing.has(id)) out.add(id);
  return out;
}

// ── detector rings (FocusDetectorPanel.tsx 52-77) ────────────────────────────

/**
 * Active indices + observed peaks → greenfield `RingInput[]`:
 *  - each claimed peak q (ix.peaks[].q_observed) → a phase-coloured ring,
 *  - each predicted-but-absent order (predicted_q with no observed peak within
 *    tol) → a ghost ring,
 *  - each leftover observed peak (claimed by no active index) → a neutral ring
 *    (bare `{ q }`).
 * No active indices → no rings (the panel falls back to plain peak rings).
 */
export function toDetectorRings(
  activeIndices: IndexEntry[],
  peaks: Peak[],
): RingInput[] {
  if (activeIndices.length === 0) return [];
  const peakQs = peaks.map((p) => p.q);
  const tol = spanTol(peakQs);
  const rings: RingInput[] = [];
  const claimed: number[] = [];
  for (const ix of activeIndices) {
    const color = phaseColor(ix.phase);
    for (const p of ix.peaks) {
      rings.push({ q: p.q_observed, color });
      claimed.push(p.q_observed);
    }
    for (const pq of ix.predicted_q) {
      const matched = peakQs.some((q) => Math.abs(q - pq) <= tol);
      if (!matched) rings.push({ q: pq, color, ghost: true });
    }
  }
  for (const q of peakQs) {
    if (!claimed.some((cq) => Math.abs(cq - q) <= tol)) rings.push({ q });
  }
  return rings;
}

// ── comb series (CombPanel.tsx) ──────────────────────────────────────────────

// PHYSICS-REVIEWED √N labelling.
// The legacy CombPanel draws teeth from `ix.predicted_q` but labels them by
// ORDER INDEX (1,2,3…), not √N. The greenfield comb wants √N labels consistent
// with the custom-index physics, so we recover them: from the index's own
// lattice we regenerate the reflection list and label each predicted_q tooth by
// the √N of the NEAREST reflection (nearest by |q|).
//
// We reconstruct the lattice via `latticeForFirstOrderOnPeak(phase, basis)`:
// `basis` IS the first-order q the backend stored (the q₁ slope; see
// customIndex.basisFor), and customRefls' first reflection sits at that q, so
// passing `basis` as the "first-order peak q" reproduces exactly the index's
// own comb. For an unknown symmetry (phase ∉ SYMS) there is no reflection list,
// so we fall back to an order-index label ("√?·i") and never throw.
function labelTeeth(phase: string, predictedQ: number[]): string[] {
  const known = SYMS[phase] !== undefined;
  if (!known || predictedQ.length === 0) {
    // Fallback: no √N basis for an unknown symmetry — label by order, no crash.
    return predictedQ.map((_, i) => `√?·${i + 1}`);
  }
  const firstQ = predictedQ.find((q) => q > 0);
  if (firstQ === undefined) return predictedQ.map((_, i) => `√?·${i + 1}`);
  const val = latticeForFirstOrderOnPeak(phase, firstQ);
  const refls = customRefls(phase, val);
  return predictedQ.map((q) => {
    if (refls.length === 0) return "√?";
    let best = refls[0]!;
    let bestD = Math.abs(refls[0]!.q - q);
    for (const r of refls) {
      const d = Math.abs(r.q - q);
      if (d < bestD) { bestD = d; best = r; }
    }
    return `√${best.N}`;
  });
}

/**
 * Active indices + observed peaks → assigned comb rows + the leftover q set.
 *
 * Mirrors CombPanel: a predicted order is `observed` when a claimed peak
 * (ix.peaks[].q_observed) sits within tol of it; `leftover` is the observed
 * peaks no active index claims (by peak_id). Per claimed tooth the residual is
 * the index's fractional residual (q_obs − q_pred)/q_pred, derived from the
 * matched IndexPeakRef (`residual` is the absolute Δq; q_pred = q_obs −
 * residual). √N labels per `labelTeeth` above. `latticeLabel` mirrors the
 * legacy "a = … Å" (rounded) from `lattice_d`; `rSquared` from `r_squared`.
 */
export function toCombSeries(
  activeIndices: IndexEntry[],
  peaks: Peak[],
): { assigned: CombSeries[]; leftover: number[] } {
  const peakQs = peaks.map((p) => p.q);
  const tol = spanTol(peakQs);

  // Leftover: observed peaks claimed (by peak_id) by no active index.
  const claimedPeakIds = new Set<number>();
  for (const ix of activeIndices) for (const p of ix.peaks) claimedPeakIds.add(p.peak_id);
  const leftover = peaks.filter((p) => !claimedPeakIds.has(p.id)).map((p) => p.q);

  const assigned: CombSeries[] = activeIndices.map((ix) => {
    const labels = labelTeeth(ix.phase, ix.predicted_q);
    // claimed observed q-values + their fractional residual, keyed for lookup.
    const claimedQs = ix.peaks.map((p) => p.q_observed);
    const teeth: CombTooth[] = ix.predicted_q.map((q, i): CombTooth => {
      // Nearest claimed observed peak within tol → this predicted order is
      // "observed". The residual comes from the matched IndexPeakRef.
      let matched: IndexEntry["peaks"][number] | undefined;
      let bestD = tol;
      for (const p of ix.peaks) {
        const d = Math.abs(p.q_observed - q);
        if (d <= bestD) { bestD = d; matched = p; }
      }
      // Also require an actual observed peak near this q (mirror CombPanel,
      // which only marks a tooth observed when a member q_observed lands near a
      // detected peak). claimedQs already are member q_observed values.
      const observed = matched !== undefined
        && claimedQs.some((cq) => Math.abs(cq - q) <= tol);
      if (observed && matched !== undefined) {
        const predicted = matched.q_observed - matched.residual;
        const frac = Math.abs(predicted) > 1e-9
          ? matched.residual / predicted
          : 0;
        return { q, label: labels[i]!, observed: true, residual: frac };
      }
      return { q, label: labels[i]!, observed: false };
    });
    return {
      phase: ix.phase,
      color: phaseColor(ix.phase),
      teeth,
      ...(ix.lattice_d != null ? { latticeLabel: `a = ${ix.lattice_d.toFixed(0)} Å` } : {}),
      ...(ix.r_squared != null ? { rSquared: ix.r_squared } : {}),
    };
  });

  return { assigned, leftover };
}

// ── custom-index modal picker + live preview/fit ─────────────────────────────

/** Picker metadata for the custom-index modal, derived from `SYMS`. Order +
 *  membership mirror the modal's `SYM_OPTIONS` (cubic Pn3m/Im3m/Ia3d, then
 *  Lamellar, Hexagonal). The lattice slider step is 1 and the unit is Å (the
 *  modal hardcodes both); `paramName` is the SymmetrySpec `param` ("a"/"d"). */
export const CUSTOM_SYMS: ReadonlyArray<{
  name: string;
  paramName: string;
  unit: string;
  min: number;
  max: number;
  step?: number;
}> = (["Pn3m", "Im3m", "Ia3d", "Lamellar", "Hexagonal"] as const).map((name) => {
  const s = SYMS[name]!;
  return { name, paramName: s.param, unit: "Å", min: s.min, max: s.max, step: 1 };
});

/**
 * Live preview comb + "lands on N of M" fit for the custom-index modal, at
 * lattice `paramValue` against `observed` peak q-values. Reuses customRefls /
 * landsOn from customIndex.ts; the per-tooth `observed`/`residual` use the SAME
 * relTol (0.022) `landsOn` does, so the displayed fit count and the lit teeth
 * agree by construction.
 */
export function customIndexPreview(
  sym: string,
  paramValue: number,
  observed: number[],
): { previewSeries: CombSeries; fit: { landed: number; total: number; snapped?: boolean } } {
  const RELTOL = 0.022; // must equal customIndex.landsOn's default relTol.
  const refls = customRefls(sym, paramValue);
  const teeth: CombTooth[] = refls.map((r): CombTooth => {
    // Nearest observed peak within the same relTol landsOn uses.
    let nearest: number | undefined;
    let bestRel = RELTOL;
    for (const pq of observed) {
      const rel = Math.abs(pq - r.q) / r.q;
      if (rel < bestRel) { bestRel = rel; nearest = pq; }
    }
    if (nearest !== undefined) {
      return { q: r.q, label: `√${r.N}`, observed: true, residual: (nearest - r.q) / r.q };
    }
    return { q: r.q, label: `√${r.N}`, observed: false };
  });
  return {
    previewSeries: { phase: sym, color: phaseColor(sym), teeth },
    fit: { landed: landsOn(sym, paramValue, observed), total: refls.length },
  };
}
