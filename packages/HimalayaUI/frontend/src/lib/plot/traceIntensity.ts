import type { Trace } from "../../api";

/**
 * Trace intensity at `q` via nearest-sample lookup (the legacy TraceViewer
 * `interpolateI`).
 *
 * Manually-added ("curation") peaks come back from the backend with
 * `intensity: null` (routes_peaks.jl stores them with no metrics computed yet).
 * Anywhere a peak's glyph is placed by its y-intensity, a null would otherwise
 * drop the glyph to the axis baseline (PlotPeaks `intensity ?? baselineI`),
 * rendering it flat on the floor instead of on the curve. Both the Focus
 * trace-model adapter AND the series/folio waterfall model run null intensities
 * through this so every peak bead sits on the data.
 *
 * Returns `null` for an empty trace (no curve to read) — callers then omit the
 * intensity (conditional spread for `exactOptionalPropertyTypes`) and the glyph
 * falls back to the baseline, which is the only honest place with no curve.
 */
export function traceIntensityAt(q: number, trace: Trace): number | null {
  const qs = trace.q;
  if (qs.length === 0) return null;
  let nearest = 0;
  for (let i = 1; i < qs.length; i++) {
    if (Math.abs(qs[i]! - q) < Math.abs(qs[nearest]! - q)) nearest = i;
  }
  const iv = trace.I[nearest];
  return iv != null && Number.isFinite(iv) ? iv : null;
}
