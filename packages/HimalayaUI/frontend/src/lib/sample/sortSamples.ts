/**
 * Sortable-column model for the contact-sheet samples table (SA-SORT).
 *
 * The page owns sort state and applies a STABLE sort to the already-filtered
 * sample list before mapping it to rows. Ingest order (the array index) is the
 * stable tiebreaker — we decorate-sort-undecorate with the original index so
 * equal keys never reorder regardless of the engine's sort stability.
 */

/** Which column the table is sorted by. `null` = default ingest order. */
export type SortKey = "sample" | "exposures" | "kept" | "status" | null;
export type SortDir = "asc" | "desc";

export interface SortState {
  key: SortKey;
  dir: SortDir;
}

/** The four sortable column keys, in header order. Tags is intentionally absent
 *  (multi-valued → no total order). The checkbox column is not a data column. */
export const SORTABLE_KEYS = ["sample", "exposures", "kept", "status"] as const;
export type SortableKey = (typeof SORTABLE_KEYS)[number];

export function isSortableKey(k: string): k is SortableKey {
  return (SORTABLE_KEYS as readonly string[]).includes(k);
}

/**
 * Status rank for the "Status" column (ascending order).
 *
 * Rationale: ascending floats the most ACTIONABLE rows to the top — a triage
 * tool should surface the work that remains, not bury it. "Not indexed" (no
 * phase call yet) is the open work, so it ranks FIRST (0); any indexed sample
 * ranks AFTER it (1). Among indexed rows the secondary order is the row's name
 * (handled by the stable name tiebreaker below), so all the not-indexed samples
 * cluster at the top in ascending and at the bottom in descending.
 *
 * `phase` here is the same nullable phase string the StatusCell reads: a
 * non-empty string ⇒ indexed; null / "" ⇒ not indexed.
 */
export function statusRank(phase: string | null | undefined): number {
  const indexed = typeof phase === "string" && phase.length > 0;
  return indexed ? 1 : 0;
}

/** Accessors that pull each sortable column's comparable value off a row datum. */
export interface SortAccessors<T> {
  name: (item: T) => string;
  exposures: (item: T) => number;
  kept: (item: T) => number;
  phase: (item: T) => string | null | undefined;
}

/** Numeric-aware, locale-aware name compare so "Sample 2" < "Sample 10". */
function compareName(a: string, b: string): number {
  return a.localeCompare(b, undefined, { numeric: true, sensitivity: "base" });
}

/**
 * Stable sort of `items` by `sort`. Returns a NEW array; the input is untouched.
 * A `null` key (or an unknown key) returns ingest order (a shallow copy).
 *
 * The primary comparator decides the column; ties fall back to the original
 * ingest index so the sort is total and stable. `dir: "desc"` negates ONLY the
 * primary comparator — the ingest tiebreaker always stays ascending so equal
 * rows keep their corpus order in both directions.
 */
export function sortSampleRows<T>(
  items: readonly T[],
  sort: SortState,
  acc: SortAccessors<T>,
): T[] {
  const decorated = items.map((item, index) => ({ item, index }));
  const { key, dir } = sort;
  if (key === null || !isSortableKey(key)) {
    return decorated.map((d) => d.item);
  }
  const sign = dir === "desc" ? -1 : 1;

  const primary = (a: T, b: T): number => {
    switch (key) {
      case "sample":
        return compareName(acc.name(a), acc.name(b));
      case "exposures":
        return acc.exposures(a) - acc.exposures(b);
      case "kept":
        return acc.kept(a) - acc.kept(b);
      case "status":
        return statusRank(acc.phase(a)) - statusRank(acc.phase(b));
    }
  };

  decorated.sort((da, db) => {
    const p = sign * primary(da.item, db.item);
    if (p !== 0) return p;
    return da.index - db.index; // ingest-order tiebreaker (always ascending)
  });
  return decorated.map((d) => d.item);
}

/**
 * The header-click toggle cycle for a single column key:
 *   inactive → ascending → descending → cleared (back to ingest/null).
 * Activating a DIFFERENT column moves the sort there at ascending.
 * Only one column is ever sorted at a time.
 */
export function nextSortState(prev: SortState, clickedKey: SortableKey): SortState {
  if (prev.key !== clickedKey) return { key: clickedKey, dir: "asc" };
  if (prev.dir === "asc") return { key: clickedKey, dir: "desc" };
  return { key: null, dir: "asc" }; // desc → cleared
}
