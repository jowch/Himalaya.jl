import { useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { useSeriesList } from "../queries";
import { SeriesFolioCard } from "../components/SeriesFolioCard";
import {
  filterSort,
  type FolioFilter,
  type FolioSort,
} from "../lib/series/folioFilter";

// Static skeleton shape for boneyard's headless capture (docs/boneyard.md
// Rule 2): a few placeholder cards in the same CSS-columns masonry the live
// page uses, so the captured bone matches the loaded layout's geometry.
const FOLIO_FIXTURE = (
  <div className="[column-count:3] [column-gap:1.25rem]">
    {[0, 1, 2, 3, 4].map((i) => (
      <div key={i} className="mb-5 break-inside-avoid rounded border border-hair bg-plate">
        <div className="h-24 bg-paper-sunk" />
        <div className="space-y-2 p-3">
          <div className="h-4 w-2/3 rounded bg-paper-sunk" />
          <div className="h-3 w-1/2 rounded bg-paper-sunk" />
        </div>
      </div>
    ))}
  </div>
);

const FILTER_CHIPS: { value: FolioFilter; label: string }[] = [
  { value: "all", label: "All" },
  { value: "transition", label: "Has transition" },
  { value: "cross", label: "Cross-experiment" },
];

const SORT_OPTIONS: { value: FolioSort; label: string }[] = [
  { value: "recent", label: "Recent" },
  { value: "variable", label: "Variable" },
  { value: "size", label: "Largest" },
];

/**
 * SeriesFolioPage — the series folio at /series (#173 / I3.3; R6 #229). A
 * corpus-wide masonry of saved series, the landing surface of the Series
 * stage. Mounted under the CorpusShell layout route.
 *
 * Read-only: owns the corpus useSeriesList() query plus client-side
 * search / filter chips / 3-way sort (the listing is small and fully
 * client-held, so no refetch). Recipe-match and cross-experiment facets are
 * present but data-starved on the current corpus — see `folioFilter.ts`.
 */
export function SeriesFolioPage(): JSX.Element {
  const navigate = useNavigate();
  const query = useSeriesList();
  const [search, setSearch] = useState("");
  const [filter, setFilter] = useState<FolioFilter>("all");
  const [sort, setSort] = useState<FolioSort>("recent");

  const series = query.data ?? [];

  const visible = useMemo(
    () => filterSort(series, { search, filter, sort }),
    [series, search, filter, sort],
  );

  const isFiltered = search.trim() !== "" || filter !== "all";

  return (
    <div data-testid="series-folio-page" className="mx-auto flex max-w-[1380px] flex-col gap-4 p-6">
      <header className="flex items-end justify-between gap-8">
        <div className="flex flex-col gap-1">
          <div className="text-[11px] font-bold uppercase tracking-[0.14em] text-print-accent">
            Folio
          </div>
          <h1 data-testid="series-folio-heading" className="text-headline text-ink">
            Saved series
          </h1>
          <p className="mt-1 max-w-[60ch] text-sm text-ink-soft">
            Every comparison you've built, across all beamtimes. Pick one up where you
            left off — or start a new one from the contact sheet.
          </p>
        </div>
        <div className="flex shrink-0 items-end gap-5">
          <div className="text-right">
            <div data-testid="series-folio-count" className="text-display leading-none text-ink">
              {visible.length}
            </div>
            <div className="mt-1 text-[10.5px] font-bold uppercase tracking-[0.08em] text-ink-faint">
              {isFiltered ? `of ${series.length} shown` : "series in the folio"}
            </div>
          </div>
          <button
            type="button"
            data-testid="series-folio-new"
            onClick={() => navigate("/series/new")}
            className="rounded-md border border-print-accent bg-print-accent px-3 py-1.5 text-xs font-semibold text-paper hover:opacity-90"
          >
            + New series
          </button>
        </div>
      </header>

      <div className="flex flex-wrap items-center gap-3.5 border-b border-hair pb-4">
        <input
          data-testid="series-folio-search"
          type="search"
          placeholder="Search series…"
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          className="min-w-[230px] rounded-md border border-hair-strong bg-plate px-3 py-1.5 text-sm text-ink"
        />
        <div className="flex gap-1.5">
          {FILTER_CHIPS.map((c) => (
            <button
              key={c.value}
              type="button"
              data-testid={`series-folio-chip-${c.value}`}
              aria-pressed={filter === c.value}
              onClick={() => setFilter(c.value)}
              className={[
                "rounded-full border px-3 py-1 text-[11.5px] font-semibold",
                filter === c.value
                  ? "border-ink bg-ink text-paper"
                  : "border-hair-strong bg-plate text-ink-soft hover:border-ink-faint",
              ].join(" ")}
            >
              {c.label}
            </button>
          ))}
        </div>
        <span className="ml-auto text-[10.5px] font-bold uppercase tracking-[0.07em] text-ink-faint">
          Sort
        </span>
        <div className="flex overflow-hidden rounded-md border border-hair-strong">
          {SORT_OPTIONS.map((s) => (
            <button
              key={s.value}
              type="button"
              data-testid={`series-folio-sort-${s.value}`}
              aria-pressed={sort === s.value}
              onClick={() => setSort(s.value)}
              className={`px-3 py-1.5 text-[11.5px] font-semibold ${
                sort === s.value ? "bg-ink text-paper" : "text-ink-faint"
              }`}
            >
              {s.label}
            </button>
          ))}
        </div>
      </div>

      {query.isError ? (
        <div data-testid="series-folio-error" className="px-4 py-8 text-sm text-ink-soft">
          Could not load series. Try reloading the page.
        </div>
      ) : (
        <Skeleton
          name="series-folio"
          className="w-full"
          loading={query.isLoading}
          stagger={50}
          transition={200}
          fixture={FOLIO_FIXTURE}
        >
          {series.length === 0 ? (
            <div data-testid="series-folio-empty" className="px-4 py-12 text-center text-sm text-ink-faint">
              No series yet. Select samples on the contact sheet to start one.
            </div>
          ) : visible.length === 0 ? (
            <div data-testid="series-folio-no-match" className="px-4 py-12 text-center text-sm text-ink-faint">
              No series match. Clear the search or filter to see the whole folio.
            </div>
          ) : (
            <div className="[column-count:1] sm:[column-count:2] lg:[column-count:3] [column-gap:1.25rem]">
              {visible.map((s, i) => (
                <SeriesFolioCard
                  key={s.id}
                  series={s}
                  figNumber={i + 1}
                  onOpen={(id) => navigate(`/series/${id}`)}
                />
              ))}
            </div>
          )}
        </Skeleton>
      )}
    </div>
  );
}
