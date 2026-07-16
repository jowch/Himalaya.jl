// focusAdapters.ts (Phase-4 Focus) — PURE (no React) functions that turn carried
// backend data (api.ts `Peak` / `IndexEntry` / `Trace`) into the prop shapes the
// greenfield "The Print" composites expect (TraceModel / RingInput / CombSeries).
//
// The custom-index physics is NOT re-implemented here — it is reused wholesale
// from `src/lib/customIndex.ts` (SYMS / customRefls / landsOn /
// latticeForFirstOrderOnPeak). The legacy React components this ports the
// presentation-free logic out of are cited inline.

import type { Peak, IndexEntry, Trace, Experiment } from "../../api";
import { phaseColor } from "../../phases";
import type { TraceModel } from "../plot/TracePlot";
import { traceIntensityAt } from "../../lib/plot/traceIntensity";
import type { PlotPeak } from "../plot/marks/PlotPeaks";
import type { RingInput, DetectorCalibration } from "../detector/detectorGeometry";
import type { CombSeries, CombTooth } from "../comb/combModel";
import {
  SYMS,
  customRefls,
  landsOn,
  latticeForFirstOrderOnPeak,
} from "../../lib/customIndex";

// The legacy span-relative q tolerance (floored at 1e-6). Only used where no
// backend ratio_position join exists: a claimless index (peaks: [] — a
// committed custom index) matching its predicted_q against observed peaks in
// toDetectorRings. Claimed reflections join by ratio_position, never by tol.
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
 * `activeIndices` param.
 *
 * A peak with no backend `intensity` (a manually-added curation peak) is
 * anchored to the trace curve at its q via {@link traceIntensityAt} so its
 * glyph sits on the data, not on the baseline. `intensity` is conditionally
 * spread (never set to `undefined`) to satisfy `exactOptionalPropertyTypes`;
 * it is omitted only when there is no curve to read (empty trace).
 */
export function toTraceModel(
  trace: Trace,
  peaks: Peak[],
  phase: string | null,
): TraceModel {
  return {
    trace,
    phase,
    peaks: peaks.map((p): PlotPeak => {
      const intensity =
        p.intensity != null ? p.intensity : traceIntensityAt(p.q, trace);
      return {
        id: p.id,
        q: p.q,
        source: p.source,
        excluded: p.excluded,
        ...(intensity != null ? { intensity } : {}),
      };
    }),
  };
}

// ── armed-edit peak click routing ────────────────────────────────────────────

/** What a click on a peak mark does, resolved from the peak's provenance. */
export type PeakClickAction =
  | { kind: "remove" }
  | { kind: "toggle-exclude"; excluded: boolean };

/**
 * Source-aware routing for a click on a peak in the armed trace editor.
 *
 * Auto peaks belong to the indexer — they are NOT deletable (the backend
 * rejects a remove on a `source:"auto"` peak), so a click toggles their
 * `excluded` flag: disable a spurious detection (it restyles struck-through in
 * place, still visible) or restore it. Manual (curation) peaks are user-authored
 * and a click removes them outright. This replaces the old plain-click=remove /
 * alt-click=exclude split, which left auto peaks un-disable-able by a plain
 * click. Returns null for an unknown id (a stale mark mid-reconcile).
 */
export function peakClickAction(
  peaks: ReadonlyArray<{ id: number; source: string; excluded: boolean }>,
  id: number,
): PeakClickAction | null {
  const p = peaks.find((pk) => pk.id === id);
  if (!p) return null;
  if (p.source === "auto") return { kind: "toggle-exclude", excluded: !p.excluded };
  return { kind: "remove" };
}

// ── detector rings (FocusDetectorPanel.tsx 52-77) ────────────────────────────

/**
 * Active indices + observed peaks → greenfield rings + their caption phases:
 *  - each claimed peak q (ix.peaks[].q_observed) → a phase-coloured ring,
 *  - each predicted-but-absent order → a ghost ring,
 *  - each leftover observed peak (claimed by no active index, by peak_id) → a
 *    neutral ring (bare `{ q }`).
 * No active indices → no rings (the panel falls back to plain peak rings).
 *
 * "Absent" is the backend's own join, not a q-tolerance rescan:
 * `IndexPeakRef.ratio_position` is the 1-based index into `predicted_q`
 * (pipeline.jl builds both from the same normalized ratio series, at most one
 * claimed peak per position), so a predicted order is absent iff its position
 * is claimed by NO ref of this index. Re-matching by q within spanTol could
 * disagree with the claim (a peak accepted at residual > tol would paint BOTH
 * an observed ring at q_observed AND a ghost at its predicted q for the same
 * reflection — the FO-RESCORE2-P2 F2 lie). The one place the q-tolerance scan
 * survives is a claimless index (peaks: [] — a committed custom index, which
 * insert_custom_index! stores with no index_peaks rows): there is no join to
 * read, so its absent orders are the predicted_q with no observed peak within
 * tol — a fully-landed custom index therefore still emits zero rings.
 *
 * `phases` is the ring-identity caption source (FO-RING): the distinct phases
 * that actually put a ring on the frame, appended in walk (= rail) order the
 * first time their index emits one. A ghost ring still carries the phase hue,
 * so it counts; an index that emits NOTHING (the fully-landed custom above)
 * contributes no phase. Deriving both from the same walk keeps the caption
 * honest by construction: it can only name hues that are on the frame.
 */
export function toDetectorRings(
  activeIndices: IndexEntry[],
  peaks: Peak[],
): { rings: RingInput[]; phases: string[] } {
  if (activeIndices.length === 0) return { rings: [], phases: [] };
  const peakQs = peaks.map((p) => p.q);
  const tol = spanTol(peakQs); // only for the claimless-index (no-join) fallback
  const rings: RingInput[] = [];
  const phases: string[] = [];
  const claimedPeakIds = new Set<number>();
  for (const ix of activeIndices) {
    const color = phaseColor(ix.phase);
    let emitted = false;
    const claimedPositions = new Set<number>(); // 0-based predicted_q indices
    for (const p of ix.peaks) {
      rings.push({ q: p.q_observed, color });
      claimedPeakIds.add(p.peak_id);
      claimedPositions.add(p.ratio_position - 1);
      emitted = true;
    }
    const hasJoin = ix.peaks.length > 0;
    for (let i = 0; i < ix.predicted_q.length; i++) {
      const pq = ix.predicted_q[i]!;
      const absent = hasJoin
        ? !claimedPositions.has(i)
        : !peakQs.some((q) => Math.abs(q - pq) <= tol);
      if (absent) {
        rings.push({ q: pq, color, ghost: true });
        emitted = true;
      }
    }
    if (emitted && !phases.includes(ix.phase)) phases.push(ix.phase);
  }
  for (const p of peaks) {
    if (!claimedPeakIds.has(p.id)) rings.push({ q: p.q });
  }
  return { rings, phases };
}

// ── detector calibration ─────────────────────────────────────────────────────

/**
 * Experiment beamline params + the RAW detector pixel size (from the image
 * route's X-Image-Width/Height headers) → a DetectorCalibration for the
 * geometry engine. Returns null unless ALL ingredients are present and finite
 * (a 0/NaN dim or null field would yield NaN/Infinity radii); null → the
 * DetectorPanel centered fallback. Pure + unit-tested; all the arithmetic and
 * guarding lives here, not in the component.
 */
export function buildDetectorCalibration(
  experiment: Experiment | undefined,
  rawSize: { w: number; h: number } | null,
): DetectorCalibration | null {
  if (!experiment || !rawSize) return null;
  const { beam_center_x, beam_center_y, pixel_size_um, energy_kev, flight_path_m } = experiment;
  const ok = (n: number | null): n is number => n !== null && Number.isFinite(n);
  if (!ok(beam_center_x) || !ok(beam_center_y) || !ok(pixel_size_um) ||
      !ok(energy_kev) || !ok(flight_path_m)) return null;
  if (!(rawSize.w > 0) || !(rawSize.h > 0)) return null;
  return {
    beamCenterPx: { x: beam_center_x, y: beam_center_y },
    imageSizePx: { w: rawSize.w, h: rawSize.h },
    sampleDistanceMm: flight_path_m * 1000,
    pixelSizeMm: pixel_size_um / 1000,
    energyKeV: energy_kev,
  };
}

// ── comb series (CombPanel.tsx) ──────────────────────────────────────────────

// PHYSICS-REVIEWED √N labelling.
// The legacy CombPanel draws teeth from `ix.predicted_q` but labels them by
// ORDER INDEX (1,2,3…), not √N. The greenfield comb wants √N labels consistent
// with the custom-index physics, so we recover N from the DIFFRACTION Q-LAW —
// NOT by matching against customRefls' truncated `SYMS.Ms` list. That list stops
// at a handful of orders, so nearest-matching clamps/mislabels every higher
// order (Pn3m √10..√16 → "√9", Im3m √14..√20 → "√12", Lamellar 6..11 → "√5"),
// and position-aligning fails too where SYMS skips an order the backend keeps
// (Hexagonal SYMS.Ms omits √11).
//
// q ∝ √N for cubic AND hexagonal (both q ∝ √(quadratic form)); q ∝ N for
// lamellar. We anchor on `predicted_q[0]` (== basis, the q₁ slope; see
// customIndex.basisFor) and the first allowed reflection's N (n1), then solve
// the q-law per order:
//   lamellar:       N = round(n1 · q/q0)     (q ∝ N)
//   cubic/hex/square: N = round(n1 · (q/q0)²) (q ∝ √N ⇒ N = N1·(q/q0)²)
// All eight canonical phases are now in SYMS (Square is kind "square", q ∝ √N
// like cubic), so each labels via the q-law branch. The order-index fallback
// ("√?·i") only fires for a symmetry genuinely absent from SYMS, and never throws.
function labelTeeth(phase: string, predictedQ: number[]): string[] {
  const spec = SYMS[phase];
  if (!spec || predictedQ.length === 0) {
    // Fallback: no q-law basis for an unknown symmetry — label by order, no crash.
    return predictedQ.map((_, i) => `√?·${i + 1}`);
  }
  const q0 = predictedQ[0]!;
  if (q0 <= 0) return predictedQ.map((_, i) => `√?·${i + 1}`);
  const val = latticeForFirstOrderOnPeak(phase, q0);
  const n1 = customRefls(phase, val)[0]?.N ?? 1; // first allowed N (e.g. Pn3m → 2)
  return predictedQ.map((q) => {
    const ratio = q / q0;
    const n = spec.kind === "lamellar"
      ? Math.round(n1 * ratio)          // q ∝ N
      : Math.round(n1 * ratio * ratio); // cubic + hex: q ∝ √N ⇒ N = N1·(q/q0)²
    return `√${n}`;
  });
}

/**
 * Active indices + observed peaks → assigned comb rows + the leftover q set.
 *
 * A predicted order `i` (0-based) is `observed` iff a claimed IndexPeakRef has
 * `ratio_position - 1 === i` — the backend's own join. `ratio_position` is the
 * 1-based index into the phase's normalized ratio series, and `predicted_q` is
 * basis × that same series in the SAME order (pipeline.jl), with at most one
 * claimed peak per position. Re-deriving the match by nearest-q within spanTol
 * was a SECOND tolerance regime that could reject claims the backend accepted
 * (a peak claimed at residual > tol rendered absent while the assignment cart
 * counted it — the FO-RESCORE2-P2 F2 lie). Refs whose position falls outside
 * `predicted_q` are ignored (defensive; should not happen).
 *
 * Per claimed tooth the residual is the SIGNED fraction
 * (q_obs − q_pred)/q_pred with q_pred = predicted_q[ratio_position − 1] — the
 * true ideal. (The backend `residual` field is abs(q_obs − ideal), so
 * reconstructing the ideal as q_obs − residual is wrong below the prediction;
 * ResidualChart plots signed fractions above/below its baseline.)
 *
 * `leftover` is the observed peaks no active index claims (by peak_id).
 * √N labels per `labelTeeth` above. `latticeLabel` mirrors the legacy
 * "a = … Å" (rounded) from `lattice_d`; `rSquared` from `r_squared`.
 */
export function toCombSeries(
  activeIndices: IndexEntry[],
  peaks: Peak[],
): { assigned: CombSeries[]; leftover: number[] } {
  // Leftover: observed peaks claimed (by peak_id) by no active index.
  const claimedPeakIds = new Set<number>();
  for (const ix of activeIndices) for (const p of ix.peaks) claimedPeakIds.add(p.peak_id);
  const leftover = peaks.filter((p) => !claimedPeakIds.has(p.id)).map((p) => p.q);

  const assigned: CombSeries[] = activeIndices.map((ix) => {
    const labels = labelTeeth(ix.phase, ix.predicted_q);
    // Join claimed refs to teeth by ratio_position (1-based into predicted_q).
    const claimedByTooth = new Map<number, IndexEntry["peaks"][number]>();
    for (const p of ix.peaks) {
      const i = p.ratio_position - 1;
      if (i >= 0 && i < ix.predicted_q.length) claimedByTooth.set(i, p);
    }
    const teeth: CombTooth[] = ix.predicted_q.map((q, i): CombTooth => {
      const claimed = claimedByTooth.get(i);
      if (claimed !== undefined) {
        const frac = Math.abs(q) > 1e-9 ? (claimed.q_observed - q) / q : 0;
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

/** Picker metadata for the custom-index modal, derived from `SYMS`. All eight
 *  canonical phases are offered (F13): the three cubics Pn3m/Im3m/Ia3d, then the
 *  inverse-micellar cubics Fm3m/Fd3m, then Lamellar, Hexagonal, and the 2D
 *  Square — cubics grouped first, low-dimensional phases last. The lattice
 *  slider step is 1 and the unit is Å (the modal hardcodes both); `paramName`
 *  is the SymmetrySpec `param` ("a"/"d"). */
export const CUSTOM_SYMS: ReadonlyArray<{
  name: string;
  paramName: string;
  unit: string;
  def: number;
  min: number;
  max: number;
  step?: number;
}> = (["Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m", "Lamellar", "Hexagonal", "Square"] as const).map((name) => {
  const s = SYMS[name]!;
  return { name, paramName: s.param, unit: "Å", def: s.def, min: s.min, max: s.max, step: 1 };
});

/**
 * Live preview comb + "N of M reflections land" fit for the custom-index modal, at
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
