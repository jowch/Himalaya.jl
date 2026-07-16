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
 *    one distinct experiment (valid because q is in absolute Å⁻¹, so traces
 *    from different experiments share one x-axis and compare directly — see
 *    docs/himalayaui-design.md §2.4).
 *  - "variable" sort orders by the recipe ordering variable, tiebroken by title
 *    so series with no variable yet (null) stay in a stable order.
 */
import type { SeriesSummary } from "../../api";

export type FolioFilter = "all" | "transition" | "cross";
export type FolioSort = "recent" | "variable" | "size";
export type FolioDir = "asc" | "desc";

export interface FolioControls {
  search: string;
  filter: FolioFilter;
  sort: FolioSort;
  /** Sort direction. Optional on input (callers may omit it — the sort's
   *  default direction is then used); `folioControlsFromParams` always fills it
   *  so the parsed shape is fully determined. */
  dir?: FolioDir;
}

/**
 * The natural direction for each sort: most-recent-first, largest-first, and
 * variable A→Z. The direction control flips this; the URL only carries `dir`
 * when it differs from this default (the absent-defaults convention).
 */
export function defaultDirForSort(sort: FolioSort): FolioDir {
  return sort === "variable" ? "asc" : "desc";
}

/** Resolve the effective direction: the explicit `dir` when present, else the
 *  sort's natural default. */
function effectiveDir(c: FolioControls): FolioDir {
  return c.dir ?? defaultDirForSort(c.sort);
}

/**
 * One lowercase haystack per series: ALL the human-meaningful text on a card,
 * using only fields that exist on `SeriesSummary` — the title, the description,
 * the author, the member phases (e.g. "Im3m"), and the recipe ordering variable
 * (e.g. "LL37 : lipid ratio"). Beamtime / experiment NAME is deliberately NOT
 * searched: SeriesSummary carries `experiment_name` only when a series does not
 * span experiments (null otherwise), so matching it would be inconsistent — the
 * folio's experiment axis is the `spans_experiments` "Cross-experiment" filter,
 * not free-text search.
 */
function searchHaystack(s: SeriesSummary): string {
  return [
    s.title,
    s.description ?? "",
    s.author_username ?? "",
    s.member_phases.join(" "),
    s.ordering_variable ?? "",
  ]
    .join(" ")
    .toLowerCase();
}

/**
 * URL ⇄ controls codec (FOL P2-3 — the house permalink convention).
 *
 * Schema: `?q=<search>&filter=transition|cross&sort=variable|size`.
 * Defaults are ABSENT from the URL: an empty search never leaves `?q=`,
 * `filter=all` / `sort=recent` are the no-param state. Unknown values fall
 * back to the defaults (a shared link with a stale token still renders).
 */
export function folioControlsFromParams(params: URLSearchParams): FolioControls {
  const filter = params.get("filter");
  const sort = params.get("sort");
  const rawDir = params.get("dir");
  const resolvedSort: FolioSort =
    sort === "variable" || sort === "size" ? sort : "recent";
  // An explicit asc/desc wins; anything else (absent or unknown token) falls
  // back to the resolved sort's natural direction so a stale link still renders.
  const dir: FolioDir =
    rawDir === "asc" || rawDir === "desc"
      ? rawDir
      : defaultDirForSort(resolvedSort);
  return {
    search: params.get("q") ?? "",
    filter: filter === "transition" || filter === "cross" ? filter : "all",
    sort: resolvedSort,
    dir,
  };
}

/** Writes `controls` onto `params` in place (set non-defaults, delete
 *  defaults) and returns it. Non-folio params are left untouched. */
export function writeFolioControlsToParams(
  params: URLSearchParams,
  c: FolioControls,
): URLSearchParams {
  if (c.search === "") params.delete("q");
  else params.set("q", c.search);
  if (c.filter === "all") params.delete("filter");
  else params.set("filter", c.filter);
  if (c.sort === "recent") params.delete("sort");
  else params.set("sort", c.sort);
  // `dir` is absent when it equals the sort's natural default (matching the
  // absent-defaults convention); only a flipped direction leaves a param.
  if (effectiveDir(c) === defaultDirForSort(c.sort)) params.delete("dir");
  else params.set("dir", effectiveDir(c));
  return params;
}

export function filterSort(
  series: SeriesSummary[],
  c: FolioControls,
): SeriesSummary[] {
  let rows = series.slice();

  const needle = c.search.trim().toLowerCase();
  if (needle !== "") {
    rows = rows.filter((s) => searchHaystack(s).includes(needle));
  }

  if (c.filter === "transition") {
    // "Has transition" ≈ spans more than one distinct phase.
    rows = rows.filter((s) => s.member_phase_count > 1);
  } else if (c.filter === "cross") {
    // Members resolve to more than one distinct experiment.
    rows = rows.filter((s) => s.spans_experiments);
  }

  // Direction: each sort has a natural default (defaultDirForSort). The
  // comparator below is written for that default; for the flipped direction we
  // reverse it. "recent" has no comparator (it rides the backend's
  // last_event_at DESC order), so we reverse the array in place instead.
  const dir = effectiveDir(c);
  const flipped = dir !== defaultDirForSort(c.sort);

  if (c.sort === "size") {
    // Default desc: largest first. localeCompare-free numeric comparator;
    // flip negates it for smallest-first.
    rows.sort((a, b) =>
      flipped ? a.member_count - b.member_count : b.member_count - a.member_count,
    );
  } else if (c.sort === "variable") {
    // Order by the recipe ordering variable; series with no variable yet (null)
    // sort last (asc) / first (desc), and ties break by title so the order is
    // stable. The whole comparator is negated for the flipped direction.
    rows.sort((a, b) => {
      const base = compareByVariable(a, b);
      return flipped ? -base : base;
    });
  } else {
    // "recent": backend already returns last_event_at DESC (the default).
    // The flipped (asc) direction is the reverse of the input order.
    if (flipped) rows.reverse();
  }

  return rows;
}

/** Stable variable comparator in the natural (asc) direction: by ordering
 *  variable A→Z with nulls last, tiebroken by title. */
function compareByVariable(a: SeriesSummary, b: SeriesSummary): number {
  const av = a.ordering_variable;
  const bv = b.ordering_variable;
  if (av !== bv) {
    if (av === null) return 1;
    if (bv === null) return -1;
    const byVar = av.localeCompare(bv);
    if (byVar !== 0) return byVar;
  }
  return a.title.localeCompare(b.title);
}
