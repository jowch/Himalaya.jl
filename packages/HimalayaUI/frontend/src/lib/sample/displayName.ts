// displayName.ts — single source of truth for "what string do we render for a sample".
import type { Sample } from "../../api";

// Uses `||` not `??` so an empty-string display_name falls through rather
// than rendering as a blank tile or a leading separator (e.g. " · run.dat"
// in comparison labels). Matches the existing logic in lib/comparison/labels.ts.
export const sampleDisplayName = (
  s: Pick<Sample, "id" | "name" | "display_name">,
): string => s.display_name || s.name || `Sample #${s.id}`;
