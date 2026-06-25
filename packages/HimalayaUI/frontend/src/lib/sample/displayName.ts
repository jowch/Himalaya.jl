// displayName.ts — single source of truth for "what string do we render for a sample".
import type { Sample } from "../../api";

// Uses `||` not `??` so an empty-string name falls through to the id-based
// fallback rather than rendering as a blank tile or a leading separator.
export const sampleDisplayName = (
  s: Pick<Sample, "id" | "name">,
): string => s.name || `Sample #${s.id}`;
