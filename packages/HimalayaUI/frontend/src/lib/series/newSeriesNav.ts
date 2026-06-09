// src/lib/series/newSeriesNav.ts
import type { NavigateFunction, Location } from "react-router-dom";

export interface NewSeriesLocationState {
  seedSampleIds: number[];
}

/**
 * Navigate to /series/new seeded with the checked sample ids. The seed is
 * passed as React Router location state — ephemeral navigational intent,
 * scoped to the history entry, not stored in Zustand or localStorage.
 */
export function navigateToNewSeries(
  sampleIds: Set<number>,
  navigate: NavigateFunction,
): void {
  navigate("/series/new", {
    state: { seedSampleIds: [...sampleIds] } satisfies NewSeriesLocationState,
  });
}

/**
 * Read the seed from the current location's state (produced by
 * `navigateToNewSeries`). Returns the array of sample ids if present and
 * well-formed, otherwise null (direct /series/new visit → full-corpus path).
 */
export function readNewSeriesSeed(location: Pick<Location, "state">): number[] | null {
  const s = location.state;
  if (s == null || typeof s !== "object") return null;
  const seed = (s as Record<string, unknown>)["seedSampleIds"];
  if (!Array.isArray(seed)) return null;
  return seed as number[];
}
