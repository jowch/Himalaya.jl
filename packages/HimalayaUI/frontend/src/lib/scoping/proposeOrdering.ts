import type { SampleTagPair, PickerSampleRow } from "../../api";

export interface OrderingRow {
  sampleId: number;
  sampleName: string;
  value: string;
  flagged: boolean; // true when the sample has no value for the ordering key
  include: boolean; // D5: per-row member include toggle, default all-in
}
export interface OrderingProposal {
  orderingKey: string | undefined; // undefined on a cold corpus
  rows: OrderingRow[];
}

/**
 * Machine proposal for the series ordering variable. Picks the most frequent
 * tag `key` across the corpus tag pairs as the ordering variable, then reads
 * each member sample's value for that key. A sample missing the key is
 * `flagged` (the human must resolve it before the confirm-and-build gate
 * opens). Cold corpus (no tags) → no key, every row flagged (#174 accepts
 * this — scoping starts cold).
 */
export function proposeOrdering(
  corpusTags: SampleTagPair[],
  samples: PickerSampleRow[],
): OrderingProposal {
  const freq = new Map<string, number>();
  for (const t of corpusTags) freq.set(t.key, (freq.get(t.key) ?? 0) + 1);
  // Deterministic: highest count, ties broken by lexicographic key. Take the
  // first entry of the sorted array — no loop. `sorted[0]` is `undefined` on a
  // cold corpus (no tags), so `orderingKey` is `undefined` there.
  const sorted = [...freq.entries()].sort(
    (a, b) => b[1] - a[1] || (a[0] < b[0] ? -1 : 1),
  );
  const [orderingKey] = sorted[0] ?? [];

  const rows: OrderingRow[] = samples.map((s) => {
    const tag = orderingKey === undefined
      ? undefined
      : s.sample.tags.find((t) => t.key === orderingKey);
    const value = tag?.value ?? "";
    return {
      sampleId: s.sample.id,
      sampleName: s.sample.display_name ?? s.sample.name ?? "",
      value,
      flagged: value === "",
      include: true, // D5: every member candidate starts included
    };
  });
  return { orderingKey, rows };
}
