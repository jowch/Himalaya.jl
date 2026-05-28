import type { OrderingProposal, OrderingRow } from "./proposeOrdering";

export interface SplitProposal {
  members: OrderingRow[]; // rows that have a parsed value for the ordering key
  looseMatches: OrderingRow[]; // rows that lack the key entirely — surfaced, not assumed
}

/**
 * Partition a proposal into the members the machine WOULD assume into the
 * series (a parsed value for the ordering key) and the loose matches it
 * NOTICED but will not assume (no value for the key). Loose matches start
 * excluded — the human adds them deliberately. Cold corpus (no orderingKey)
 * → everything is a loose match (the page shows its cold fallback). Pure:
 * surfaces proposeOrdering's output, does not change its key selection.
 */
export function splitProposal(p: OrderingProposal): SplitProposal {
  if (p.orderingKey === undefined) {
    return { members: [], looseMatches: p.rows.map((r) => ({ ...r, include: false })) };
  }
  const members: OrderingRow[] = [];
  const looseMatches: OrderingRow[] = [];
  for (const r of p.rows) {
    if (r.value !== "") members.push(r);
    else looseMatches.push({ ...r, include: false });
  }
  return { members, looseMatches };
}

/** Present a tag key as a human label (underscores/hyphens → spaces). */
export function humanizeKey(key: string | undefined): string {
  if (key === undefined || key === "") return "—";
  return key.replace(/[_-]+/g, " ");
}
