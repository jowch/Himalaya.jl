import type { IndexEntry } from "../../api";

export interface PhaseRead {
  dominant: string | null; // highest-scored phase
  coexist: string | null; // second distinct phase, if any (two-phase coexistence)
}

/**
 * Read a sample's dominant phase from its candidate indices: highest `score`
 * wins; the next distinct phase (if any) is the coexistence partner that
 * drives the preview strip's two-phase gradient. Pure; null score sorts last.
 */
export function dominantPhase(indices: IndexEntry[]): PhaseRead {
  const ranked = [...indices].sort((a, b) => (b.score ?? -Infinity) - (a.score ?? -Infinity));
  const dominant = ranked[0]?.phase ?? null;
  const coexist = ranked.find((i) => i.phase !== dominant)?.phase ?? null;
  return { dominant, coexist };
}
