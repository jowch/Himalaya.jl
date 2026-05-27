import { useMemo, useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { useSeriesList } from "../queries";
import { SeriesFolioCard } from "../components/SeriesFolioCard";

type SortMode = "recency" | "title";

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

/**
 * SeriesFolioPage — the series folio at /series (#173 / I3.3). A corpus-wide
 * masonry of saved series, the landing surface of the Series stage. Mounted
 * under the CorpusShell layout route, like SamplesPage/LoupePage.
 *
 * Read-only: owns the corpus useSeriesList() query, a client-side title
 * search, and a recency/title sort toggle. The listing is small and fully
 * client-held, so search/sort are local (no refetch). Recency is the default
 * (the backend already returns last_event_at DESC). Beamtime/phase filter
 * chips from the mockup are deferred (#173 scope note).
 */
export function SeriesFolioPage(): JSX.Element {
  const navigate = useNavigate();
  const query = useSeriesList();
  const [search, setSearch] = useState("");
  const [sort, setSort] = useState<SortMode>("recency");

  const series = query.data ?? [];

  const visible = useMemo(() => {
    const needle = search.trim().toLowerCase();
    const filtered =
      needle === ""
        ? series
        : series.filter((s) => s.title.toLowerCase().includes(needle));
    if (sort === "title") {
      return [...filtered].sort((a, b) => a.title.localeCompare(b.title));
    }
    // recency: backend order is already last_event_at DESC; preserve it.
    return filtered;
  }, [series, search, sort]);

  return (
    <div data-testid="series-folio-page" className="flex flex-col gap-4 p-6">
      <header className="flex items-end justify-between gap-6">
        <div className="flex flex-col gap-1">
          <div className="text-xs font-semibold uppercase tracking-wide text-print-accent">
            Series
          </div>
          <p className="text-sm text-ink-faint">
            Every comparison you've built, across all beamtimes.
          </p>
        </div>
        <div className="text-right">
          <div data-testid="series-folio-count" className="text-2xl font-medium text-ink">
            {series.length}
          </div>
          <div className="text-[10px] uppercase tracking-wide text-ink-faint">
            in the folio
          </div>
        </div>
      </header>

      <div className="flex items-center gap-3 border-b border-hair pb-3">
        <input
          data-testid="series-folio-search"
          type="search"
          placeholder="Search series…"
          value={search}
          onChange={(e) => setSearch(e.target.value)}
          className="rounded border border-hair-strong bg-plate px-3 py-1.5 text-sm text-ink"
        />
        <div className="ml-auto flex overflow-hidden rounded border border-hair-strong">
          <button
            type="button"
            data-testid="series-folio-sort-recency"
            aria-pressed={sort === "recency"}
            onClick={() => setSort("recency")}
            className={`px-3 py-1.5 text-xs ${sort === "recency" ? "bg-ink text-paper" : "text-ink-faint"}`}
          >
            Recent
          </button>
          <button
            type="button"
            data-testid="series-folio-sort-title"
            aria-pressed={sort === "title"}
            onClick={() => setSort("title")}
            className={`px-3 py-1.5 text-xs ${sort === "title" ? "bg-ink text-paper" : "text-ink-faint"}`}
          >
            Title
          </button>
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
          ) : (
            <div className="[column-count:1] sm:[column-count:2] lg:[column-count:3] [column-gap:1.25rem]">
              {visible.map((s) => (
                <SeriesFolioCard key={s.id} series={s} onOpen={(id) => navigate(`/series/${id}`)} />
              ))}
            </div>
          )}
        </Skeleton>
      )}
    </div>
  );
}
