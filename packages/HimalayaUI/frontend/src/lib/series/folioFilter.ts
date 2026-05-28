/**
 * folioFilter — pure search / filter / sort for the series folio controls
 * (R6, finding F-D; mockup `series-folio.html` controls row).
 *
 * Kept out of the page component so the (data-shaped) decisions are unit-
 * tested without rendering. Two facets are present-but-data-starved on the
 * current corpus and documented as such:
 *
 *  - "cross" (cross-experiment) is not derivable from the `SeriesSummary`
 *    listing (no sample→experiment join on this path), so it yields nothing.
 *    The chip is wired per F-D and will light up once that join exists.
 *  - "variable" sort would key on the ordering variable, which the listing
 *    does not carry — so it sorts by title (the only stable text key on the
 *    summary). Labelled "Variable" per the mockup.
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
    // No cross-experiment signal on the listing (see module note).
    rows = [];
  }

  if (c.sort === "size") {
    rows.sort((a, b) => b.member_count - a.member_count);
  } else if (c.sort === "variable") {
    rows.sort((a, b) => a.title.localeCompare(b.title));
  }
  // "recent": backend already returns last_event_at DESC; preserve order.

  return rows;
}
