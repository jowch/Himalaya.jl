import { useState } from "react";
import { useNavigate } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { FolioHeader } from "../components/FolioHeader";
import { Gallery } from "../components/Gallery";
import { SeriesCard } from "../components/SeriesCard";
import {
  SearchInput,
  FilterChip,
  SegmentedControl,
  EmptyState,
} from "../ui";
import type { SegmentOption } from "../ui";
import { useSeriesList, useSeries, useSeriesTraces } from "../../queries";
import type { SeriesSummary } from "../../api";
import { filterSort } from "../../lib/series/folioFilter";
import type { FolioControls, FolioFilter, FolioSort } from "../../lib/series/folioFilter";
import { membersToSegments, toCardChrome } from "./folioAdapters";
import { toWaterfallRows } from "../waterfall/waterfallModel";

// ── Boneyard fixture — inline placeholder so the Skeleton has a defined
// fallback while no folio.bones.json capture has been made yet.
// Token-only classes; no inline appearance literals (design-guard clean). ────
const FOLIO_FIXTURE = (
  <div className="columns-1 sm:columns-2 lg:columns-3 gap-5">
    {[0, 1, 2].map((i) => (
      <div key={i} className="break-inside-avoid mb-5">
        <div className="bg-plate border border-hair rounded overflow-hidden">
          <div className="px-2 pt-2 h-28 bg-paper-sunk" />
          <div className="px-4 py-3">
            <div className="h-3 w-20 rounded-md bg-paper-sunk mb-3" />
            <div className="h-5 w-full rounded-md bg-paper-sunk mb-2" />
            <div className="h-3 w-32 rounded-md bg-paper-sunk" />
          </div>
        </div>
      </div>
    ))}
  </div>
);

// ── Sort options ──────────────────────────────────────────────────────────────
const SORT_OPTIONS: ReadonlyArray<SegmentOption<FolioSort>> = [
  { value: "recent", label: "Recent" },
  { value: "variable", label: "Variable" },
  { value: "size", label: "Largest" },
];

// ── Per-card data container (page glue; NOT a print/components composite) ─────
function FolioCard({
  summary,
  position,
  onOpen,
  now,
}: {
  summary: SeriesSummary;
  position: number;
  onOpen: () => void;
  now: Date;
}): JSX.Element {
  const detail = useSeries(summary.id);
  const tracesQ = useSeriesTraces(summary.id);
  const members = detail.data?.members ?? [];
  const rows = toWaterfallRows(members, tracesQ.data ?? {});
  const segments = membersToSegments(members);
  const chrome = toCardChrome(summary, position, now);
  return (
    <SeriesCard
      rows={rows}
      segments={segments}
      {...chrome}
      onClick={onOpen}
    />
  );
}

/**
 * SeriesFolioPage (greenfield) — the series folio wall at /series.
 *
 * Assembled from src/print composites + the series-folio mockup: a
 * `PageFrame width="folio"` body with the folio head (kicker + serif
 * "Saved series" title + series count), a controls row (search + filter chips
 * + sort segmented control), and a `Gallery` of `FolioCard` containers.
 *
 * Carried logic only (useSeriesList / filterSort); no legacy presentation.
 */
export function SeriesFolioPage(): JSX.Element {
  const navigate = useNavigate();
  const listQ = useSeriesList();
  const summaries = listQ.data ?? [];
  const now = new Date();

  const [controls, setControls] = useState<FolioControls>({
    search: "",
    filter: "all",
    sort: "recent",
  });

  const visible = filterSort(summaries, controls);

  function setFilter(filter: FolioFilter): void {
    setControls((c) => ({ ...c, filter }));
  }

  // Honest error surface: a failed list fetch settles with data=undefined, which
  // would otherwise read as "No series match" (telling the user to clear a filter
  // when the request actually failed). Handle it before the normal body so a
  // failure never masquerades as a zero-results state.
  if (listQ.isError) {
    return (
      <PageFrame width="folio" className="px-10 py-8">
        <FolioHeader
          kicker="Folio"
          title="Saved series"
          count={0}
          countLabel="series in the folio"
          className="mb-5"
        />
        <EmptyState
          title="Couldn't load the folio"
          body="The series list failed to load. Try reloading the page."
        />
      </PageFrame>
    );
  }

  return (
    <PageFrame width="folio" className="px-10 py-8">
      {/* ── Header ────────────────────────────────────────────────────── */}
      <FolioHeader
        kicker="Folio"
        title="Saved series"
        subtitle="Every comparison you've built, across all beamtimes. Pick one up where you left off, or select samples on the contact sheet to start a new one."
        count={summaries.length}
        countLabel="series in the folio"
        className="mb-5"
      />

      {/* ── Controls ──────────────────────────────────────────────────── */}
      <div className="flex items-center gap-3.5 flex-wrap pb-4 mb-5 border-b border-hair">
        <SearchInput
          value={controls.search}
          onChange={(v) => setControls((c) => ({ ...c, search: v }))}
          placeholder="Search series…"
        />
        <div className="flex gap-1.5">
          <FilterChip
            label="All"
            active={controls.filter === "all"}
            onClick={() => setFilter("all")}
          />
          <FilterChip
            label="Has transition"
            active={controls.filter === "transition"}
            onClick={() => setFilter("transition")}
          />
          <FilterChip
            label="Cross-experiment"
            active={controls.filter === "cross"}
            onClick={() => setFilter("cross")}
          />
        </div>
        <div className="flex-1" />
        <SegmentedControl
          aria-label="Sort series"
          options={SORT_OPTIONS}
          value={controls.sort}
          onChange={(v) => setControls((c) => ({ ...c, sort: v }))}
        />
      </div>

      {/* ── Gallery ───────────────────────────────────────────────────── */}
      <Skeleton
        name="folio"
        className="block"
        loading={listQ.isLoading}
        stagger={50}
        transition={200}
        fixture={FOLIO_FIXTURE}
        fallback={
          <div className="p-8 text-sm text-ink-soft">Loading series…</div>
        }
      >
        <Gallery
          empty={
            <EmptyState
              title="No series match"
              body="Clear the search or filter to see the whole folio."
            />
          }
        >
          {visible.map((s, i) => (
            <FolioCard
              key={s.id}
              summary={s}
              position={i + 1}
              now={now}
              onOpen={() => navigate(`/series/${s.id}`)}
            />
          ))}
        </Gallery>
      </Skeleton>
    </PageFrame>
  );
}
