/**
 * Per-q tolerance for matching peaks across optimistic / remote / server views.
 * Mirrors the Julia formula `max(1e-6, abs(q) * 0.001)` used by
 * `hash_peak_set_from_db` and curation lookups in pipeline.jl. Diverges from a
 * pure 1e-6 absolute tolerance for q values far from 1.0.
 */
export function peakQTol(q: number): number {
  return Math.max(1e-6, Math.abs(q) * 0.001);
}
