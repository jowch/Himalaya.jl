import type { CorpusSample } from "../api";
import { naturalCompare } from "../lib/naturalSort";

export type SampleSort = "measured" | "name";

export interface SampleSortOption {
  sort: SampleSort;
  name: string;
  hint: string;
}

export const SAMPLE_SORTS: readonly SampleSortOption[] = [
  { sort: "measured", name: "Measured", hint: "Order the images were taken" },
  { sort: "name", name: "Name", hint: "A to Z, with 2 before 10" },
];

export const SORT_WORDS = ["order", "alphabetical", "name", "time", "measured"];

type Sortable = Pick<CorpusSample, "id" | "name" | "setup_id" | "measured_at">;

function byMeasured(a: Sortable, b: Sortable): number {
  const at = a.measured_at ?? null;
  const bt = b.measured_at ?? null;
  if (at === bt) return a.id - b.id;
  if (at === null) return 1;
  if (bt === null) return -1;
  return at < bt ? -1 : 1;
}

function byName(a: Sortable, b: Sortable): number {
  return naturalCompare(a.name, b.name) || byMeasured(a, b);
}

export function sortSamples<T extends Sortable>(
  samples: readonly T[],
  sort: SampleSort,
  setupRank: (s: T) => number,
): T[] {
  const within = sort === "name" ? byName : byMeasured;
  return [...samples].sort((a, b) => setupRank(a) - setupRank(b) || within(a, b));
}
