import { useEffect, useMemo, useRef, useState } from "react";
import { useLocation, useNavigate, useSearchParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { SheetTable } from "../components/SheetTable";
import { SampleTableRow } from "../components/SampleTableRow";
import { CullBar } from "../components/CullBar";
import { Kicker, KbLegend, ProgressBar, ComposeBar, EmptyState, Button } from "../ui";
import {
  useCorpusSamples,
  useCorpusExposures,
  useScreenedProgress,
  useExperiments,
  useSetExposureStatusBatch,
} from "../../queries";
import { navigateToNewSeries } from "../../lib/series/newSeriesNav";
import { useFloatingDock } from "../shell/floatingDock";
import {
  resolveExperimentFilter,
  UNKNOWN_BEAMTIME_LABEL,
} from "../../lib/experimentFilter";
import { suppressGlobalKeys } from "../../lib/keys";
import { showToast } from "../../lib/toast";
import { toSampleRowModel } from "./samplesAdapters";
import {
  sortSampleRows,
  nextSortState,
  isSortableKey,
  type SortKey,
  type SortDir,
  type SortableKey,
} from "../../lib/sample/sortSamples";
import { sampleDisplayName } from "../../lib/sample/displayName";

// Boneyard fixture — a static skeleton shaped to the SheetTable rows region so
// the headless capture CLI measures the contact-sheet body. Four placeholder
// rows of token-only bars (no inline appearance literals: design-guard clean).
const CONTACT_SHEET_FIXTURE = (
  <div className="bg-paper-sunk">
    {[0, 1, 2, 3].map((i) => (
      <div key={i} className="border-b border-hair px-4 py-5">
        <div className="h-3 w-40 rounded-md bg-paper-sunk" />
        <div className="mt-3 h-12 w-full rounded-md bg-paper-sunk" />
      </div>
    ))}
  </div>
);

/**
 * Single shared predicate (SA-F3): how many distinct samples a selection
 * spans, counted EXACTLY the way the Drop/Keep/Restore handler routes the
 * batch — distinct owners among the ids `ownerOf` can map (unmappable ids
 * are skipped by the handler, so they don't count here either). Both the
 * CullBar disclosure and the toast receipt go through this one function so
 * the promise and the action can never count differently.
 */
function selectionSpread(
  selected: ReadonlySet<number>,
  ownerOf: ReadonlyMap<number, number>,
): number {
  const owners = new Set<number>();
  for (const id of selected) {
    const sampleId = ownerOf.get(id);
    if (sampleId !== undefined) owners.add(sampleId);
  }
  return owners.size;
}

/**
 * SamplesPage (greenfield) — the contact-sheet table at /samples.
 *
 * Reimplemented from src/print composites + the sample-table mockup: a
 * `PageFrame width="sheet"` body with the contact-sheet head (kicker + serif
 * beamtime title + screened progress), a slotted `SheetTable` of
 * `SampleTableRow`s, a footer keyboard legend, and a floating `CullBar`.
 *
 * Carried logic only (queries / mutators / the row adapter); no legacy
 * presentation. Selection is page-owned (a Set of exposure ids); a Shift-anchor
 * ref drives the contiguous-range cull within a single sample's exposure order.
 */
export function SamplesPage(): JSX.Element {
  const navigate = useNavigate();
  const [searchParams, setSearchParams] = useSearchParams();

  // ── filter + title (carried logic) ──────────────────────────────────────────
  const corpusQuery = useCorpusSamples();
  const experimentsQuery = useExperiments();
  // SA-F5: the SAME resolver the topbar select uses — the page body and the
  // select must agree on whether the URL names a real beamtime.
  const filter = resolveExperimentFilter(
    searchParams.get("beamtime"),
    experimentsQuery.data,
    corpusQuery.data,
  );
  const beamtime = filter.id;

  const samples = useMemo(
    () => corpusQuery.data ?? [],
    [corpusQuery.data],
  );
  const filtered = useMemo(
    () =>
      beamtime === undefined
        ? samples
        : samples.filter((s) => s.experiment_id === beamtime),
    [samples, beamtime],
  );

  // The unknown title never repeats the raw id from the address; the loading
  // fallback (`experiment N`) only shows before the lists settle the verdict.
  const title =
    beamtime === undefined
      ? "The corpus"
      : filter.unknown
        ? UNKNOWN_BEAMTIME_LABEL
        : (filter.name ?? `experiment ${beamtime}`);

  function clearBeamtimeFilter(): void {
    setSearchParams((prev) => {
      const next = new URLSearchParams(prev);
      next.delete("beamtime");
      return next;
    });
  }

  const corpusExposures = useCorpusExposures(filtered);
  const { screened, total } = useScreenedProgress(filtered);
  const batch = useSetExposureStatusBatch();

  // ── sort (SA-SORT, page-owned, URL-persisted) ───────────────────────────────
  // The sort key+dir live in the address the SAME way the beamtime filter does
  // (useSearchParams), so a sorted view is a shareable permalink. An unknown /
  // absent ?sort falls back to null = ingest order; ?dir defaults to asc.
  const sortKeyParam = searchParams.get("sort");
  const sortKey: SortKey =
    sortKeyParam && isSortableKey(sortKeyParam) ? sortKeyParam : null;
  const sortDir: SortDir = searchParams.get("dir") === "desc" ? "desc" : "asc";
  const sort = { key: sortKey, dir: sortDir };

  function applySort(key: SortableKey): void {
    const next = nextSortState(sort, key);
    setSearchParams((prev) => {
      const params = new URLSearchParams(prev);
      if (next.key === null) {
        params.delete("sort");
        params.delete("dir");
      } else {
        params.set("sort", next.key);
        params.set("dir", next.dir);
      }
      return params;
    });
  }

  // Apply the stable sort to the FILTERED list before it becomes rows. The
  // exposures/kept accessors read the same per-sample exposure rows the row
  // adapter does, so the column you sort by is exactly the column you see.
  const sortedSamples = useMemo(() => {
    if (sortKey === null) return filtered;
    return sortSampleRows(filtered, { key: sortKey, dir: sortDir }, {
      name: (s) => sampleDisplayName(s),
      exposures: (s) => (corpusExposures.byId.get(s.id) ?? []).length,
      kept: (s) =>
        (corpusExposures.byId.get(s.id) ?? []).filter((e) => e.status === "accepted")
          .length,
      phase: (s) => s.phase,
    });
  }, [filtered, sortKey, sortDir, corpusExposures.byId]);

  // LO-NEXT: hand the loupe the exact visible order (sorted + beamtime-filtered)
  // through router state, so its prev/next walk matches what's on screen. (The
  // "exposures"/"kept" sorts can't be reconstructed from the loupe's data, so
  // passing the resolved id list is the honest source.)
  const sampleOrder = useMemo(
    () => sortedSamples.map((s) => s.id),
    [sortedSamples],
  );

  // ── return-focus from the loupe (LO-FOCUSRET, WCAG 2.4.3) ────────────────────
  // The loupe Escape/back navigation carries the originating sample id in router
  // state. On this remount we resolve it to its (1-based) row in the CURRENT
  // sort and hand SheetTable a one-shot focus target so focus lands back on that
  // row's Sample cell instead of <body>. Captured once on mount: a later sort or
  // SSE re-render must not re-yank focus (SheetTable's seed is itself one-shot,
  // but pinning the row index here keeps it stable against a re-sort too).
  const location = useLocation();
  const focusOnMountRowRef = useRef<number | undefined>(undefined);
  const focusSeedConsumed = useRef(false);
  if (!focusSeedConsumed.current) {
    focusSeedConsumed.current = true;
    const st = location.state as { focusSampleId?: number } | null;
    if (st?.focusSampleId != null) {
      const idx = sortedSamples.findIndex((s) => s.id === st.focusSampleId);
      if (idx >= 0) focusOnMountRowRef.current = idx + 1;
    }
  }
  const focusOnMountRow = focusOnMountRowRef.current;

  // ── selection state (page-owned) ────────────────────────────────────────────
  // Exposure-grain cull selection (drives CullBar → Drop/Restore).
  const [selected, setSelected] = useState<Set<number>>(() => new Set());
  const anchorRef = useRef<{ sampleId: number; exposureId: number } | null>(null);
  const shiftRef = useRef(false);

  // Sample-grain pick selection — DISTINCT from the exposure-grain cull above.
  // Drives the ComposeBar → "+ New series" carry; never merged with `selected`.
  const [checkedSamples, setCheckedSamples] = useState<Set<number>>(() => new Set());

  function toggleSampleCheck(sampleId: number): void {
    setCheckedSamples((prev) => {
      const next = new Set(prev);
      if (next.has(sampleId)) next.delete(sampleId);
      else next.add(sampleId);
      return next;
    });
  }

  // SA-STALESELECT: a beamtime-filter change re-scopes the visible corpus, so a
  // working selection from the prior scope is no longer on screen to act on or
  // uncheck. Carrying it forward made the CullBar/ComposeBar counts, the cull
  // toast receipt, and the "+ New series" carry lie about what's actionable
  // (e.g. "1 sample" with no visible row, or a drop toast counting frames whose
  // sample left the filter). Clear both grains whenever the scope changes; on
  // mount this is a no-op (both sets start empty).
  useEffect(() => {
    setSelected(new Set());
    setCheckedSamples(new Set());
    anchorRef.current = null;
  }, [beamtime]);

  // LA-COLLIDE: the CullBar / ComposeBar are opaque, fixed, bottom-centre and
  // sit above the global InfrastructureBanner's z-index — so while either shows
  // it would paint over a "Saving…" banner. Publish the dock-lane occupancy so
  // the banner steps aside to a corner; release it on unmount so leaving the
  // sheet recentres the banner.
  const setCenterLaneOccupied = useFloatingDock((s) => s.setCenterLaneOccupied);
  const dockLaneOccupied = selected.size > 0 || checkedSamples.size > 0;
  useEffect(() => {
    setCenterLaneOccupied(dockLaneOccupied);
    return () => setCenterLaneOccupied(false);
  }, [dockLaneOccupied, setCenterLaneOccupied]);

  useEffect(() => {
    function onShiftDown(e: KeyboardEvent): void {
      if (e.key === "Shift") shiftRef.current = true;
    }
    function onShiftUp(e: KeyboardEvent): void {
      if (e.key === "Shift") shiftRef.current = false;
    }
    window.addEventListener("keydown", onShiftDown);
    window.addEventListener("keyup", onShiftUp);
    return () => {
      window.removeEventListener("keydown", onShiftDown);
      window.removeEventListener("keyup", onShiftUp);
    };
  }, []);

  // exposureId → owning sampleId (for routing a batch mutate back to its sample).
  const ownerOf = useMemo(() => {
    const m = new Map<number, number>();
    for (const s of filtered)
      for (const e of corpusExposures.byId.get(s.id) ?? []) m.set(e.id, s.id);
    return m;
  }, [filtered, corpusExposures.byId]);

  function toggleSelect(sampleId: number, exposureId: number): void {
    setSelected((prev) => {
      const next = new Set(prev);
      const anchor = anchorRef.current;
      if (shiftRef.current && anchor && anchor.sampleId === sampleId) {
        // Contiguous range within THIS sample's exposure order.
        const ids = (corpusExposures.byId.get(sampleId) ?? []).map((e) => e.id);
        const a = ids.indexOf(anchor.exposureId);
        const b = ids.indexOf(exposureId);
        if (a >= 0 && b >= 0) {
          const lo = Math.min(a, b);
          const hi = Math.max(a, b);
          for (let i = lo; i <= hi; i++) next.add(ids[i]!);
          return next;
        }
      }
      // Single toggle + (re)anchor.
      if (next.has(exposureId)) next.delete(exposureId);
      else next.add(exposureId);
      anchorRef.current = { sampleId, exposureId };
      return next;
    });
  }

  function batchSet(status: "accepted" | "rejected" | null): void {
    // SA-STALESELECT: only act on exposures still in the visible scope (an owner
    // resolves in `ownerOf`, which is built from `filtered`). The toast then
    // counts what was actually mutated, never the raw selection size — a
    // destructive action must never claim a larger drop/keep count than it made.
    const targets = [...selected].filter((id) => ownerOf.has(id));
    for (const id of targets) {
      batch.mutate({ sampleId: ownerOf.get(id)!, exposureId: id, status });
    }
    const n = targets.length;
    // Consequential, batch-scale change → visible toast. Drop, Keep and the
    // symmetric restore all announce so each action and its undo are equally
    // legible. Restore sends null unconditionally: it clears BOTH verdicts.
    if (n > 0) {
      const verb =
        status === "rejected" ? "dropped" : status === "accepted" ? "kept" : "restored";
      // SA-F3: the receipt matches the bar's promise — a cross-sample batch
      // discloses its spread with the same suffix, under the same predicate.
      const spread = selectionSpread(new Set(targets), ownerOf);
      const suffix = spread > 1 ? ` across ${spread} samples` : "";
      showToast(`${n} frame${n === 1 ? "" : "s"} ${verb}${suffix}`, "success");
    }
    setSelected(new Set());
    anchorRef.current = null;
  }

  // ── keyboard: X drops, K keeps, Esc clears (only with a selection) ──────────
  useEffect(() => {
    function onKeyDown(e: KeyboardEvent): void {
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      if (suppressGlobalKeys(e)) return;
      if (selected.size === 0) return;
      if (e.key === "x" || e.key === "X") batchSet("rejected");
      else if (e.key === "k" || e.key === "K") batchSet("accepted");
      else if (e.key === "Escape") {
        setSelected(new Set());
        anchorRef.current = null;
      }
    }
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [selected, ownerOf]);

  // SA-F2: a double-clicked frame carries its exposure id into the loupe via
  // ?exposure= so the loupe opens AT that frame. A search param (not router
  // state) keeps the frame on refresh/permalink, the house slug-permalink
  // direction. Frameless openings (name click) omit it → loupe default frame.
  function loupeHref(id: number, exposureId?: number): string {
    const params = new URLSearchParams();
    if (beamtime !== undefined) params.set("beamtime", String(beamtime));
    if (exposureId !== undefined) params.set("exposure", String(exposureId));
    const qs = params.toString();
    return `/samples/loupe/${id}${qs ? `?${qs}` : ""}`;
  }

  return (
    <PageFrame width="sheet" className="px-10 py-8">
      <div data-testid="samples-page">
        {/* ── Head ──────────────────────────────────────────────────────────── */}
        <div className="flex items-end justify-between gap-8 mb-5">
          <div>
            <Kicker tone="accent">Contact sheet</Kicker>
            <h1 className="text-display text-ink">{title}</h1>
            <p className="text-body text-ink-soft mt-2 max-w-[62ch]">
              Flip the frames and drop the ones with flares or artifacts. Tags
              are a light, optional note on what each sample is. The ordering
              variable is named later, when you scope a series.
            </p>
          </div>
          <div className="shrink-0 text-right">
            <div className="text-headline-lg text-ink leading-none">
              {screened}
              <b className="font-medium text-ink-soft"> / {total}</b>
            </div>
            <Kicker tone="faint" className="mt-0.5">samples screened</Kicker>
            <ProgressBar
              value={screened}
              total={total}
              label="samples screened"
              className="w-40 mt-1"
            />
          </div>
        </div>

        {/* ── Table ─────────────────────────────────────────────────────────── */}
        {corpusQuery.isError ? (
          <div data-testid="samples-error">
            <EmptyState
              title="Could not load the corpus"
              body="The sample list request failed."
              action={
                <Button
                  variant="outline"
                  disabled={corpusQuery.isFetching}
                  onClick={() => void corpusQuery.refetch()}
                >
                  Reload the corpus
                </Button>
              }
            />
          </div>
        ) : filter.unknown ? (
          // SA-F5: an address filtering by a beamtime that names nothing gets
          // an honest empty state, never the bare "No samples in this beamtime"
          // placeholder (which would imply the beamtime exists). The single
          // action embodies the way forward: clear the filter.
          <div data-testid="samples-unknown-beamtime">
            <EmptyState
              title={UNKNOWN_BEAMTIME_LABEL}
              body="The address filters by a beamtime that is not in this corpus."
              action={
                <Button variant="outline" onClick={clearBeamtimeFilter}>
                  Show all experiments
                </Button>
              }
            />
          </div>
        ) : (
          <Skeleton
            name="contact-sheet"
            className="block"
            loading={corpusQuery.isLoading}
            stagger={50}
            transition={200}
            fixture={CONTACT_SHEET_FIXTURE}
            fallback={
              <div className="p-8 text-sm text-ink-soft">Loading samples…</div>
            }
          >
            <SheetTable
              checkboxColumn
              roving
              dataRowCount={sortedSamples.length}
              {...(focusOnMountRow !== undefined ? { focusOnMountRow } : {})}
              sort={sort}
              onSort={applySort}
              empty={
                <div className="p-10 text-center text-ink-soft">
                  No samples{" "}
                  {beamtime === undefined ? "in the corpus" : "in this beamtime"}{" "}
                  yet.
                </div>
              }
            >
              {sortedSamples.map((s, i) => {
                const loadedExposures = corpusExposures.byId.get(s.id);
                const m = toSampleRowModel(s, loadedExposures);
                // SA-ZEROEXP: distinguish a CONFIRMED-empty sample (query resolved
                // with []) from one still loading (undefined). Only the former
                // gets the terminal "No exposures" status + suppressed door.
                const noExposures = loadedExposures !== undefined && loadedExposures.length === 0;
                return (
                  <SampleTableRow
                    key={s.id}
                    rowIndex={i + 1}
                    name={m.name}
                    sampleId={m.sampleId}
                    screened={m.screened}
                    exposures={m.exposures}
                    kept={m.kept}
                    total={m.total}
                    dropped={m.dropped}
                    noExposures={noExposures}
                    tags={m.tags}
                    {...(m.phase !== undefined ? { phase: m.phase } : {})}
                    checked={checkedSamples.has(s.id)}
                    onCheck={() => toggleSampleCheck(s.id)}
                    selectedExposureIds={selected}
                    onSelectExposure={(id) => toggleSelect(s.id, id)}
                    onActivateExposure={(id) =>
                      navigate(loupeHref(s.id, id), { state: { sampleOrder } })
                    }
                    onOpenLoupe={() =>
                      navigate(loupeHref(s.id), { state: { sampleOrder } })
                    }
                    onOpenFocus={() => navigate(`/sample/${s.id}`)}
                  />
                );
              })}
            </SheetTable>
          </Skeleton>
        )}

        {/* ── Footer keyboard legend ────────────────────────────────────────── */}
        <KbLegend
          className="mt-4"
          shortcuts={[
            { keyLabel: "↑ ↓ ← →", description: "move between cells" },
            { keyLabel: "Enter", description: "open, sort a header, or enter a cell" },
            { keyLabel: "click", description: "select a frame" },
            { keyLabel: "⇧ click", description: "extend the range" },
            { keyLabel: "X", description: "drop the selected frames" },
            { keyLabel: "K", description: "keep the selected frames" },
            { keyLabel: "double-click", description: "open the loupe" },
            { keyLabel: "Esc", description: "exit a cell or clear" },
            { keyLabel: "⌘K", description: "find a sample" },
          ]}
        />
      </div>

      {/* ── Floating cull bar (page root) ─────────────────────────────────────── */}
      <CullBar
        count={selected.size}
        sampleCount={selectionSpread(selected, ownerOf)}
        show={selected.size > 0}
        onReject={() => batchSet("rejected")}
        onKeep={() => batchSet("accepted")}
        onRestore={() => batchSet(null)}
        onClear={() => {
          setSelected(new Set());
          anchorRef.current = null;
        }}
      />

      {/* ── Floating sample-grain compose bar (sits above CullBar) ─────────────── */}
      <ComposeBar
        count={checkedSamples.size}
        show={checkedSamples.size > 0}
        onNewSeries={() => navigateToNewSeries(checkedSamples, navigate)}
        onClear={() => setCheckedSamples(new Set())}
      />
    </PageFrame>
  );
}
