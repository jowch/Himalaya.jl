/**
 * folioFilter — pure search / filter / sort for the series folio controls
 * (R6, finding F-D; mockup `series-folio.html` controls row).
 *
 * Kept out of the page component so the (data-shaped) decisions are unit-
 * tested without rendering. Both the cross-experiment filter and the variable
 * sort are now backed by real listing fields (M2 — `series_listing` projects
 * `spans_experiments` via a member→exposure→sample experiment-id join, and
 * `ordering_variable` straight off the series row):
 *
 *  - "cross" (cross-experiment) keeps only series whose members span more than
 *    one distinct experiment (valid because q is absolute — see redesign-notes
 *    architecture decision 1).
 *  - "variable" sort orders by the recipe ordering variable, tiebroken by title
 *    so series with no variable yet (null) stay in a stable order.
 */
import type { SeriesSummary } from "../../api";

export type FolioFilter = "all" | "transition" | "cross";
export type FolioSort = "recent" | "variable" | "size";

export interface FolioControls {
  search: string;
  filter: FolioFilter;
  sort: FolioSort;
}

export function filterSort(
  series: SeriesSummary[],
  c: FolioControls,
): SeriesSummary[] {
  let rows = series.slice();

  const needle = c.search.trim().toLowerCase();
  if (needle !== "") {
    rows = rows.filter((s) => s.title.toLowerCase().includes(needle));
  }

  if (c.filter === "transition") {
    // "Has transition" ≈ spans more than one distinct phase.
    rows = rows.filter((s) => s.member_phase_count > 1);
  } else if (c.filter === "cross") {
    // Members resolve to more than one distinct experiment.
    rows = rows.filter((s) => s.spans_experiments);
  }

  if (c.sort === "size") {
    rows.sort((a, b) => b.member_count - a.member_count);
  } else if (c.sort === "variable") {
    // Order by the recipe ordering variable; series with no variable yet (null)
    // sort last, and ties break by title so the order is stable.
    rows.sort((a, b) => {
      const av = a.ordering_variable;
      const bv = b.ordering_variable;
      if (av !== bv) {
        if (av === null) return 1;
        if (bv === null) return -1;
        const byVar = av.localeCompare(bv);
        if (byVar !== 0) return byVar;
      }
      return a.title.localeCompare(b.title);
    });
  }
  // "recent": backend already returns last_event_at DESC; preserve order.

  return rows;
}
