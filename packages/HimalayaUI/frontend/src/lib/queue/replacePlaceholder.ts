/**
 * Replace a single negative-id optimistic placeholder with the server row,
 * deduping against any concurrent SSE that already inserted the server id.
 *
 * Algorithm:
 *   - Walk the list once.
 *   - The FIRST item with id < 0 whose `matches(item)` predicate returns true
 *     is replaced in place with `serverItem`. Order is preserved.
 *   - Any item where `isDuplicate(item)` returns true is dropped. By default
 *     this matches `item.id === serverItem.id` (the common case: an SSE-wins
 *     race already inserted the server row). Override for cases like
 *     `peakAdd` where the id namespace is shared across kinds and dedup
 *     must be scoped — e.g. only drop the duplicate if it's a manual peak.
 *   - If no placeholder matched, `serverItem` is appended.
 *
 * Invariants:
 *   - At most one row matching `isDuplicate` survives in the output (it
 *     becomes the serverItem if it also matched the placeholder predicate).
 *   - At most one placeholder is consumed per call. Stale placeholders left
 *     over (rare) survive — a later mutator call or invalidate sweeps them.
 *   - Rows that do not match `isDuplicate` are preserved verbatim.
 *
 * Behavior-equivalence note:
 *   This drops EVERY `isDuplicate` row regardless of position, whereas the
 *   hand-rolled loops it replaced seeded their `seen` set as they iterated —
 *   so a duplicate appearing BEFORE the placeholder was kept and the server
 *   row dropped. That divergence is unreachable: the optimistic placeholder
 *   is written synchronously in `onMutate` before any SSE-inserted server row
 *   can be appended, so the placeholder always precedes its duplicate. The
 *   equivalence holds only because of that ordering guarantee.
 */
export function replacePlaceholder<T extends { id: number }>(
  list: T[],
  serverItem: T,
  matches: (item: T) => boolean,
  options?: { isDuplicate?: (item: T) => boolean },
): T[] {
  const isDup = options?.isDuplicate ?? ((item: T) => item.id === serverItem.id);
  let replaced = false;
  const out: T[] = [];
  for (const item of list) {
    if (isDup(item)) continue;
    if (!replaced && item.id < 0 && matches(item)) {
      out.push(serverItem);
      replaced = true;
      continue;
    }
    out.push(item);
  }
  if (!replaced) out.push(serverItem);
  return out;
}
