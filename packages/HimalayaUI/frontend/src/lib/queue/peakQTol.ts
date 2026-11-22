/**
 * Tolerance for optimistic display only. The server's matching is grid-aware
 * and nearest-first (peak_matching.jl); the refetch in refreshPeaks.ts settles
 * the peak list, so do not rely on this to decide what is excluded.
 */
export function peakQTol(q: number): number {
  return Math.max(1e-6, Math.abs(q) * 0.001);
}
