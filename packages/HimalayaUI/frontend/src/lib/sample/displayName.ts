import type { Sample } from "../../api";

export const sampleDisplayName = (
  s: Pick<Sample, "id" | "name">,
): string => s.name || `Sample ${s.id}`;
