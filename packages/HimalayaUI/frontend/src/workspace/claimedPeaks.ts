import type { IndexEntry, Peak } from "../api";
import { peakKey } from "../lib/peakRef";
import { phaseColor } from "../phases";

const sameQ = (a: number, b: number): boolean => Math.abs(a - b) <= 1e-9 + 1e-6 * Math.abs(b);

export function claimedPeaks(hovered: readonly IndexEntry[], peaks: readonly Peak[]): Map<string, string> | null {
  if (hovered.length === 0) return null;
  const out = new Map<string, string>();
  for (const c of hovered) {
    for (const ref of c.peaks) {
      const p = peaks.find((pk) => sameQ(pk.q, ref.q_observed));
      if (p && !out.has(peakKey(p))) out.set(peakKey(p), phaseColor(c.phase));
    }
  }
  return out;
}
