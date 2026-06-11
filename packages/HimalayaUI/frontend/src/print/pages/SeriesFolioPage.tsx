import { useEffect, useMemo, useRef, useState } from "react";
import type { RefObject } from "react";
import { useNavigate, useSearchParams } from "react-router-dom";
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
  Button,
} from "../ui";
import type { SegmentOption } from "../ui";
import { useSeriesList, useSeries, useSeriesTraces } from "../../queries";
import type { SeriesSummary } from "../../api";
import {
  filterSort,
  folioControlsFromParams,
  writeFolioControlsToParams,
} from "../../lib/series/folioFilter";
import type { FolioControls, FolioFilter, FolioSort } from "../../lib/series/folioFilter";
import { membersToSegments, stableFigNumbers, toCardChrome } from "./folioAdapters";
import { toWaterfallRows } from "../waterfall/waterfallModel";

// ── Folio loading skeleton (FOL-BONES) — a hand-rolled card-shaped placeholder
// grid doing double duty: it is the Skeleton `fixture` (what the boneyard
// capture CLI measures when a folio.bones.json capture is eventually made) AND
// the runtime `fallback` (what renders today, since src/bones/registry.ts has
// no "folio" entry yet — boneyard's documented no-bones path, Rule 3). Until a
// real capture lands, loading shows these card bones instead of bare text.
// Token-only classes; no inline appearance literals (design-guard clean). ────
const FOLIO_SKELETON = (
  <div
    className="columns-1 sm:columns-2 lg:columns-3 gap-5"
    data-testid="folio-bones-fallback"
    role="status"
    aria-label="Loading series"
  >
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

// ── Viewport-lazy gate (FOL-N+1, page-local by design) ───────────────────────
// One-shot latch: `near` flips true the first time the element comes within
// PREFETCH_MARGIN of the viewport, then stays true (queries stay mounted; no
// refetch churn on scroll-away). Fails open — no IntersectionObserver (JSDOM,
// very old browsers) means every card loads immediately, which is exactly the
// pre-fix behaviour and keeps a tiny corpus rendering at once.
const PREFETCH_MARGIN = "400px 0px";

function useNearViewport(): [RefObject<HTMLDivElement>, boolean] {
  const ref = useRef<HTMLDivElement>(null);
  const [near, setNear] = useState(false);
  useEffect(() => {
    if (near) return;
    const el = ref.current;
    if (el === null || typeof IntersectionObserver === "undefined") {
      setNear(true);
      return;
    }
    const io = new IntersectionObserver(
      (entries) => {
        if (entries.some((e) => e.isIntersecting)) setNear(true);
      },
      { rootMargin: PREFETCH_MARGIN },
    );
    io.observe(el);
    return () => io.disconnect();
  }, [near]);
  return [ref, near];
}

// ── Per-card data container (page glue; NOT a print/components composite) ─────
// The card CHROME (fig label, title, counts, provenance) is fully derivable
// from the LIST summary; only the figure rows + phase strip need the detail
// and trace fetches. Those queries stay disabled (id=undefined → enabled:false)
// until the card nears the viewport, so a corpus-scale folio fires 2×(visible)
// requests, not 2N.
function FolioCard({
  summary,
  figNumber,
  onOpen,
  now,
}: {
  summary: SeriesSummary;
  /** Stable corpus number (stableFigNumbers); 0 for drafts (ignored). */
  figNumber: number;
  onOpen: () => void;
  now: Date;
}): JSX.Element {
  const [ref, near] = useNearViewport();
  const detail = useSeries(near ? summary.id : undefined);
  const tracesQ = useSeriesTraces(near ? summary.id : undefined);
  const members = detail.data?.members ?? [];
  const rows = toWaterfallRows(members, tracesQ.data ?? {});
  const segments = membersToSegments(members);
  const chrome = toCardChrome(summary, figNumber, now);

  // FOL-HONEST-DERIVED: compute figure readiness so the card never renders
  // "No clear phase" (an affirmative false scientific claim) while detail/traces
  // are loading or have errored.
  const figureState: "ready" | "pending" | "error" =
    detail.isError || tracesQ.isError
      ? "error"
      : detail.data !== undefined && tracesQ.data !== undefined
        ? "ready"
        : "pending";

  return (
    <div ref={ref}>
      <SeriesCard
        rows={rows}
        segments={segments}
        figureState={figureState}
        {...chrome}
        onClick={onOpen}
      />
    </div>
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

  // FOL P2-3: the controls live in the URL (the house permalink convention —
  // ?q=&filter=&sort=, defaults absent), so reload/share keeps the view.
  // Writes use replace, not push: a keystroke must not become a history entry
  // (same reason SamplesPage's beamtime select doesn't spam Back).
  const [searchParams, setSearchParams] = useSearchParams();
  const controls = useMemo(
    () => folioControlsFromParams(searchParams),
    [searchParams],
  );

  function updateControls(patch: Partial<FolioControls>): void {
    setSearchParams(
      (prev) =>
        writeFolioControlsToParams(new URLSearchParams(prev), {
          ...folioControlsFromParams(prev),
          ...patch,
        }),
      { replace: true },
    );
  }

  const visible = filterSort(summaries, controls);

  // Stable fig numbers (FOL-FIGNUM): computed over the WHOLE committed corpus,
  // not the filtered view — a card keeps its number under every sort/filter.
  const figNumbers = useMemo(() => stableFigNumbers(summaries), [summaries]);

  function setFilter(filter: FolioFilter): void {
    updateControls({ filter });
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
          count={null}
          countLabel="series in the folio"
          className="mb-5"
        />
        <EmptyState
          title="Couldn't load the folio"
          body="The series list failed to load."
          action={
            <Button
              variant="outline"
              disabled={listQ.isFetching}
              onClick={() => void listQ.refetch()}
            >
              Try again
            </Button>
          }
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
        count={listQ.isLoading ? null : summaries.length}
        countLabel="series in the folio"
        className="mb-5"
      />

      {/* ── Controls ──────────────────────────────────────────────────── */}
      <div className="flex items-center gap-3.5 flex-wrap pb-4 mb-5 border-b border-hair">
        <SearchInput
          value={controls.search}
          onChange={(v) => updateControls({ search: v })}
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
          onChange={(v) => updateControls({ sort: v })}
        />
      </div>

      {/* ── Gallery ───────────────────────────────────────────────────── */}
      <Skeleton
        name="folio"
        className="block"
        loading={listQ.isLoading}
        stagger={50}
        transition={200}
        fixture={FOLIO_SKELETON}
        fallback={FOLIO_SKELETON}
      >
        <Gallery
          empty={
            summaries.length === 0 ? (
              // Genuinely empty folio (nothing saved yet, regardless of what's
              // typed in the controls): honest first-run state with a door to
              // the creation path, never the filtered no-match masquerade.
              <EmptyState
                title="No series yet"
                body="The folio holds every comparison you save. New series start from samples selected on the contact sheet."
                action={
                  <Button variant="outline" onClick={() => navigate("/samples")}>
                    Open the contact sheet
                  </Button>
                }
              />
            ) : (
              <EmptyState
                title="No series match"
                body="Clear the search or filter to see the whole folio."
                action={
                  <Button
                    variant="outline"
                    // Resets search + filter (and drops their URL params);
                    // the sort — and its param — survive. Pinned semantic.
                    onClick={() => updateControls({ search: "", filter: "all" })}
                  >
                    Show the whole folio
                  </Button>
                }
              />
            )
          }
        >
          {visible.map((s) => (
            <FolioCard
              key={s.id}
              summary={s}
              figNumber={figNumbers.get(s.id) ?? 0}
              now={now}
              onOpen={() => navigate(`/series/${s.id}`)}
            />
          ))}
        </Gallery>
      </Skeleton>
    </PageFrame>
  );
}
