import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { SeriesPlate } from "../components/SeriesPlate";
import { BuilderRail } from "../components/BuilderRail";
import { MemberList } from "../components/MemberList";
import { IconButton, Button, EmptyState, Input, GripHandle } from "../ui";
import { useDragReorder } from "../components/useDragReorder";
import {
  useSeries,
  useSeriesTraces,
  useCorpusSamples,
  useCorpusPickerSamples,
  useSaveSeries,
  useCommitSeriesPlate,
} from "../../queries";
import { useAppState } from "../../state";
import * as api from "../../api";
import type { Series, SeriesMember } from "../../api";
import { toWaterfallRows, waterfallQDomain } from "../waterfall/waterfallModel";
import { memberRowLabel } from "../../lib/series/memberRead";
import {
  membersToMemberData,
  groupingSummary,
  legendPhasesOf,
  addableSamples,
  recipeRowView,
} from "./builderAdapters";
import { buildSeriesSaveBody } from "../../lib/series/buildSeriesSaveBody";
import { buildPlateFromRecipe } from "../../lib/series/buildPlateFromRecipe";
import { buildMultiTraceExportSpec } from "../../lib/figure-export/adapters/multiTraceAdapter";
import { ExportButton } from "../components/ExportButton";
import { useFigureExport } from "../components/useFigureExport";
import { showToast } from "../../lib/toast";

type Scale = "log" | "lin";

// Token-only skeleton fixture (no inline appearance literals — design-guard clean).
const BUILDER_FIXTURE = (
  <div className="grid grid-cols-[1fr_336px] gap-0">
    <div className="bg-plate border border-hair rounded m-4 h-72" />
    <div className="border-l border-hair p-4">
      <div className="h-4 w-1/2 rounded bg-paper-sunk" />
    </div>
  </div>
);

/**
 * SeriesBuilderPage (greenfield) — the always-"Compose" series builder at
 * /series/:id. A LAZY-DRAFT surface: the committed series renders read-only;
 * the first recipe edit (title / description / membership / ordering) silently
 * starts a draft via `startSeriesDraftFromSeries`; "Confirm series" runs a
 * Save→Commit chain then discards the draft (stay, now read); Cancel discards.
 * View controls (offset / scale / grouping / annotations) are local and never
 * start a draft.
 *
 * Assembled from src/print composites (SeriesPlate + BuilderRail + MemberList)
 * + the carried logic plane (queries / state machine / body builders /
 * AnnotationToggles / figure-export). No legacy presentation is imported.
 */
export function SeriesBuilderPage(): JSX.Element {
  const navigate = useNavigate();
  const id = Number(useParams<{ id: string }>().id);
  const seriesQ = useSeries(Number.isFinite(id) ? id : undefined);
  const tracesQ = useSeriesTraces(Number.isFinite(id) ? id : undefined);
  const corpusQ = useCorpusSamples();
  // Sample → indexing-exposure resolution source for the Confirm chain's
  // recipe→plate step (BU-RECIPENOOP). Loaded WITH the page (not lazily at
  // Confirm) so the resolution map is ready by Confirm time — lazy loading
  // would stall the chain on a fetch the user never sees. Cost: one corpus
  // GET /api/picker-samples per builder visit, shared (same query key) with
  // the scoping worksheet's picker.
  const pickerQ = useCorpusPickerSamples();
  // sample.id → indexing_exposure_id (the scoping page's resolution
  // precedent). The picker's indexing_exposure_id encodes the SAME
  // representative-exposure rule the backend's create-path resolver uses
  // (highest-id selected, else highest-id, else null).
  const exposureBySample = useMemo(() => {
    const m = new Map<number, number | null>();
    for (const r of pickerQ.data ?? []) m.set(r.sample.id, r.indexing_exposure_id);
    return m;
  }, [pickerQ.data]);
  // Confirm is gated on the picker having LOADED: without it every sample is
  // "unresolvable" and Confirm would honestly have to publish an empty plate.
  const resolverReady = pickerQ.data !== undefined;

  // ── Draft (lazy) + view-pref state ──────────────────────────────────────
  const draft = useAppState((s) => s.seriesDraft);
  const startDraft = useAppState((s) => s.startSeriesDraftFromSeries);
  const discardDraft = useAppState((s) => s.discardSeriesDraft);
  const setTitle = useAppState((s) => s.setSeriesDraftTitle);
  const addSample = useAppState((s) => s.addSeriesSample);
  const removeSample = useAppState((s) => s.removeSeriesSample);
  const reorderSample = useAppState((s) => s.reorderSeriesSample);
  const showPeakTicks = useAppState((s) => s.showPeakTicks);
  const showPeakLabels = useAppState((s) => s.showPeakLabels);
  const setShowPeakTicks = useAppState((s) => s.setShowPeakTicks);
  const setShowPeakLabels = useAppState((s) => s.setShowPeakLabels);

  // ── Local view state (never starts a draft) ─────────────────────────────
  const [offset, setOffset] = useState(1.2);
  const [scale, setScale] = useState<Scale>("log");
  const [hoveredKey, setHoveredKey] = useState<string | undefined>(undefined);

  // ── Confirm chain (Save → await fresh cache → Commit → discard) ─────────
  const save = useSaveSeries();
  const commit = useCommitSeriesPlate();
  const stage = useRef<"idle" | "saving" | "awaiting-fresh" | "committing">("idle");
  // Confirm-time snapshot of the series query's dataUpdatedAt. The commit only
  // fires once the cache holds a series NEWER than this watermark. That is a
  // recency gate, not a provenance gate: normally the newer data is our own
  // save's reconciliation, but under concurrent writes it is whatever the
  // server's latest projection holds (LWW — the commit route itself is
  // last-write-wins; conflict UI cancelled by decision).
  const watermark = useRef(0);
  // Confirm-time snapshot of the sample→exposure resolution map. The chain
  // reads THIS ref (not pickerQ directly) so a picker refetch/error mid-chain
  // cannot stall or skew the commit step — the recipe being saved was drafted
  // against exactly this resolution state.
  const resolveRef = useRef<Map<number, number | null>>(new Map());
  // How many recipe samples the LAST commit left off the plate (no resolvable
  // exposure). Read by the commit-success effect so the terminal toast tells
  // the truth instead of announcing a clean "Series confirmed".
  const leftOutRef = useRef(0);

  const series = seriesQ.data;
  // The active draft only counts when it targets THIS series.
  const liveDraft = draft !== null && series !== undefined && draft.id === series.id ? draft : null;

  // LOAD-BEARING: the page never calls save.reset()/commit.reset(). A lingering
  // save.isSuccess/commit.isSuccess === true between Confirm runs is INERT
  // because every chain effect below is `stage`-ref gated — they fire only while
  // stage.current is the matching phase, and onConfirm re-arms the chain by
  // flipping stage to "saving" itself. Do NOT add a reset() or remove the stage
  // gating: a reset would race the SSE-vs-HTTP deferred resolution, and dropping
  // the gate would let a stale isSuccess re-fire the commit on the next
  // unrelated re-render.
  //
  // Save landed → resolve the plate FROM THE FRESH RECIPE in the SERIES QUERY
  // CACHE, never from the stale draft and never from `save.data`. Under the
  // queue's SSE-wins race the own-op SSE confirmation resolves `save.data` as
  // the synthesized PARTIAL ({id}, no members/state —
  // saveSeriesMutator.synthesizeFromSse), so gating the commit on
  // isFullSeries(save.data) stalls the chain on live backends
  // (P0 BU-CONFIRMSTALL). The cache is the one source that is correct on BOTH
  // race outcomes: the mutator's onSuccess either setQueryData's the full HTTP
  // response or invalidates the key so TanStack refetches the canonical
  // projection — either way the series key gets a full Series and its
  // dataUpdatedAt advances past the Confirm-time watermark.
  //
  // COMMIT BODY (P0 BU-RECIPENOOP): the PATCH persists the recipe but does NOT
  // rebuild `series_members` from it, and POST /commit takes members verbatim
  // from the body — so committing `fresh.members` re-publishes the OLD plate
  // byte-for-byte (reorders never surface; an added sample never becomes a
  // member). buildPlateFromRecipe resolves `fresh.samples` (the just-saved,
  // position-ordered recipe) into the plate: each sample → its indexing
  // exposure (Confirm-time snapshot in resolveRef), display props carried over
  // from the old member with the same exposure_id, new members server-default.
  //
  // Single transition effect: "saving" + save.isSuccess → "awaiting-fresh";
  // "awaiting-fresh" + fresh full series in cache → "committing" + mutate.
  // One effect (not two) so the same pass that observes save.isSuccess also
  // re-checks the cache — on the HTTP-wins path setQueryData lands BEFORE the
  // mutation flips to success, so a split second effect keyed only on
  // dataUpdatedAt could miss its trigger and stall.
  useEffect(() => {
    if (stage.current === "saving" && save.isSuccess) stage.current = "awaiting-fresh";
    if (stage.current !== "awaiting-fresh") return;
    const fresh = seriesQ.data;
    if (fresh === undefined || !api.isFullSeries(fresh)) return;
    if (seriesQ.dataUpdatedAt <= watermark.current) return;
    stage.current = "committing";
    const { members, unresolvedSampleIds } = buildPlateFromRecipe(
      fresh.samples,
      fresh.members,
      (sampleId) => resolveRef.current.get(sampleId),
    );
    leftOutRef.current = unresolvedSampleIds.length;
    commit.mutate({ id, members });
  }, [save.isSuccess, seriesQ.data, seriesQ.dataUpdatedAt, id, commit]);

  // Stall exit: the refetch we are awaiting ERRORED (useSeries has retry:
  // false, so this surfaces fast). Reset so Confirm re-arms, and tell the
  // truth: the PATCH landed server-side; only the confirm step is missing.
  // Stage-ref guarded → single-fire (the first run flips stage to "idle").
  // save.isSuccess is in the deps so an error that PRE-DATES the stage flip
  // (isError already true when we enter awaiting-fresh) still gets a pass.
  useEffect(() => {
    if (stage.current !== "awaiting-fresh" || !seriesQ.isError) return;
    stage.current = "idle";
    showToast("Order saved, but confirming failed. Try Confirm again.", "error");
  }, [seriesQ.isError, save.isSuccess]);

  // Commit landed → drop the draft; stay (now read state).
  useEffect(() => {
    if (stage.current !== "committing" || !commit.isSuccess) return;
    stage.current = "idle";
    discardDraft();
    // Consequential terminal success of the Save→Commit chain → visible toast.
    // Honest variant when the plate resolution had to leave samples out
    // (BU-RECIPENOOP policy: commit the resolvable members, SAY what was
    // dropped — a sample with no exposure has nothing renderable to plate,
    // matching the backend create-path resolver, so blocking the whole commit
    // would dead-end the series instead).
    const leftOut = leftOutRef.current;
    showToast(
      leftOut === 0
        ? "Series confirmed"
        : leftOut === 1
          ? "Confirmed. 1 sample has no usable exposure and was left out."
          : `Confirmed. ${leftOut} samples have no usable exposure and were left out.`,
      "success",
    );
  }, [commit.isSuccess, discardDraft]);

  // Either error → reset so the user can retry. Stage-ref guarded like its
  // success sibling: the toast fires once per in-flight attempt (stage !=
  // "idle") and can't double-fire if both save.error and commit.error resolve
  // in the same cycle — the first run flips stage to "idle", so the re-run
  // early-returns.
  useEffect(() => {
    if (stage.current === "idle") return;
    if (save.error || commit.error) {
      stage.current = "idle";
      // The rail already shows a role=alert div; also surface it assertively as
      // a toast so the failure reaches a user whose focus is elsewhere.
      showToast("Couldn't confirm the series. Try again.", "error");
    }
  }, [save.error, commit.error]);

  // ── Honest states ───────────────────────────────────────────────────────
  // BU-BADID: a non-numeric :id (NaN) is knowable synchronously from the route
  // param — the gated query never enables, so without this branch the page
  // would sit on the loading skeleton forever. It exits through the SAME
  // not-found surface as a missing numeric id (the FocusPage routeStatus
  // "unknown" precedent), one shared branch for both.
  if (!Number.isFinite(id) || seriesQ.isError) {
    return (
      <PageFrame width="builder" className="px-8 py-8">
        <div data-testid="builder-not-found">
          <EmptyState
            title="Couldn't load this series"
            body="It may have been deleted, the address may be wrong, or the request failed."
            action={
              <Button variant="outline" onClick={() => navigate("/series")}>
                Back to the series folio
              </Button>
            }
          />
        </div>
      </PageFrame>
    );
  }

  // ── Recipe-edit handlers (lazy-draft: ensure then mutate) ───────────────
  const ensureDraft = (): void => {
    if (!liveDraft && series) startDraft(series);
  };
  const onEditTitle = (value: string): void => {
    ensureDraft();
    setTitle(value);
  };
  const onAddSample = (sampleId: number): void => {
    ensureDraft();
    addSample(sampleId);
  };
  // Guard-enforce the draft invariant rather than relying on the recipe editor
  // only rendering under a live draft (layout-enforced). `ensureDraft()` is a
  // no-op when a draft already exists.
  const onRemoveSample = (rowId: number): void => {
    ensureDraft();
    removeSample(rowId);
  };
  const onReorderSample = (from: number, to: number): void => {
    ensureDraft();
    reorderSample(from, to);
  };

  const onConfirm = (): void => {
    // `resolverReady` mirrors the rail's disabled gate: without the picker
    // projection every recipe sample is "unresolvable" and the chain would
    // publish an empty plate.
    if (!liveDraft || stage.current !== "idle" || !resolverReady) return;
    // Snapshot the resolution map for the whole chain (see resolveRef).
    resolveRef.current = exposureBySample;
    // Watermark BEFORE the save fires: any cache update at or before this
    // instant predates the save and must not be committed.
    watermark.current = seriesQ.dataUpdatedAt;
    stage.current = "saving";
    save.mutate({ id, ...buildSeriesSaveBody(liveDraft) });
  };
  const onCancel = (): void => {
    stage.current = "idle";
    discardDraft();
  };

  // Busy = the chain is in flight, derived from the SAME stage/mutation
  // sources the chain effects read (never a parallel flag). Each stage clause
  // carries its own truthful-exit term (the error that ends that stage) so the
  // busy register reverts IN THE SAME RENDER the failure surfaces — the exits
  // flip the stage REF in an effect, which triggers no re-render of its own,
  // so a plain `stage.current !== "idle"` would leave a lying "Confirming…"
  // stuck next to the error notice until some unrelated render.
  const confirmBusy =
    save.isPending ||
    commit.isPending ||
    // The saving clause also exits when the series query's error PRE-DATES
    // the save settling: in that render the ref still reads "saving" (the
    // effects flip it to awaiting-fresh -> idle without re-rendering), so
    // the awaiting-fresh clause below never gets to apply its exit.
    (stage.current === "saving" && save.error == null && !(save.isSuccess && seriesQ.isError)) ||
    (stage.current === "awaiting-fresh" && !seriesQ.isError) ||
    (stage.current === "committing" && commit.error == null);
  const chainError = save.error || commit.error;

  return (
    <Skeleton
      name="series-builder"
      className="block"
      loading={seriesQ.isLoading || series === undefined}
      fixture={BUILDER_FIXTURE}
      fallback={<div className="p-8 text-sm text-ink-soft">Loading series…</div>}
    >
      {series && (
        <BuilderBody
          series={series}
          tracesById={tracesQ.data ?? {}}
          corpus={corpusQ.data ?? []}
          liveDraft={liveDraft}
          offset={offset}
          onOffsetChange={setOffset}
          scale={scale}
          onScaleChange={setScale}
          hoveredKey={hoveredKey}
          onHoverRow={setHoveredKey}
          showPeakTicks={showPeakTicks}
          showPeakLabels={showPeakLabels}
          onToggleTicks={() => setShowPeakTicks(!showPeakTicks)}
          onToggleLabels={() => setShowPeakLabels(!showPeakLabels)}
          onEditTitle={onEditTitle}
          onAddSample={onAddSample}
          onRemoveSample={onRemoveSample}
          onReorderSample={onReorderSample}
          ensureDraft={ensureDraft}
          onConfirm={onConfirm}
          onCancel={onCancel}
          confirmBusy={confirmBusy}
          confirmReady={resolverReady}
          resolverError={pickerQ.isError}
          resolverLoading={!resolverReady && !pickerQ.isError}
          chainError={chainError}
        />
      )}
    </Skeleton>
  );
}

interface BuilderBodyProps {
  series: Series;
  tracesById: Record<number, api.Trace>;
  corpus: api.CorpusSample[];
  liveDraft: ReturnType<typeof useAppState.getState>["seriesDraft"];
  offset: number;
  onOffsetChange: (v: number) => void;
  scale: Scale;
  onScaleChange: (s: Scale) => void;
  hoveredKey: string | undefined;
  onHoverRow: (key?: string) => void;
  showPeakTicks: boolean;
  showPeakLabels: boolean;
  onToggleTicks: () => void;
  onToggleLabels: () => void;
  onEditTitle: (value: string) => void;
  onAddSample: (sampleId: number) => void;
  onRemoveSample: (rowId: number) => void;
  onReorderSample: (from: number, to: number) => void;
  ensureDraft: () => void;
  onConfirm: () => void;
  onCancel: () => void;
  confirmBusy: boolean;
  /** Picker projection loaded — the recipe→plate resolution source is ready.
   *  Confirm stays disabled until then (an unresolvable recipe would commit
   *  an empty plate). */
  confirmReady: boolean;
  /** Picker projection FAILED — Confirm can never become ready this session. */
  resolverError: boolean;
  resolverLoading: boolean;
  chainError: unknown;
}

function BuilderBody({
  series,
  tracesById,
  corpus,
  liveDraft,
  offset,
  onOffsetChange,
  scale,
  onScaleChange,
  hoveredKey,
  onHoverRow,
  showPeakTicks,
  showPeakLabels,
  onToggleTicks,
  onToggleLabels,
  onEditTitle,
  onAddSample,
  onRemoveSample,
  onReorderSample,
  ensureDraft,
  onConfirm,
  onCancel,
  confirmBusy,
  confirmReady,
  resolverError,
  resolverLoading,
  chainError,
}: BuilderBodyProps): JSX.Element {
  // Render model: the figure plate ALWAYS shows the committed plate (members);
  // a draft edits the recipe (membership/title), which only re-resolves the
  // plate after a Save→Commit round-trip. So the waterfall/MemberList read the
  // committed members; the draft drives only the title + the editable rail.
  const members: SeriesMember[] = series.members;
  const rows = toWaterfallRows(members, tracesById);
  // MemberList contract: rows in DISPLAY order top-down, but the waterfall
  // paints display order bottom-up — reverse so the rail's top row is the
  // plate's top trace (membersToMemberData returns a fresh array; in-place
  // reverse is safe). Keys are member ids, so hover-sync is order-agnostic.
  const memberData = membersToMemberData(members).reverse();
  const legendPhases = legendPhasesOf(members);

  // Peak-annotation gating: ticks/labels annotate indexed-peak anchors. A
  // form-factor / unindexed series has no anchors on any row → nothing to
  // annotate, so the annotation toggles are inert (not armed, disabled).
  // controls-don't-lie: don't arm a toggle that has nothing to show.
  const seriesHasPeaks = rows.some((r) => r.anchors.length > 0);

  const effectiveTitle = liveDraft ? liveDraft.title : series.title;

  // ── Figure export (plate-head ExportButton, wired via useFigureExport) ───
  // WYSIWYG contract (BU-EXPORTDIVERGE): the plate footnote promises "what
  // you compose is what you publish", so every axis the plate renders flows
  // into the spec from the SAME sources the plate reads — phase identity
  // (byPhase via the shared memberRead.dominantPhase predicate), the plate's
  // padded q-domain, the log/linear toggle, the offset slider, the annotation
  // flags, and the plate's row-label register. Publication furniture (title
  // block, white background, margins) is the only sanctioned delta.
  const exportSpec = useCallback(
    (): import("../../lib/figure-export/types").ExportSpec =>
      buildMultiTraceExportSpec({
        members,
        traces: new Map(Object.entries(tracesById).map(([k, v]) => [Number(k), v])),
        comparisonTitle: effectiveTitle,
        xDomain: waterfallQDomain(rows),
        showPeakTicks,
        showPeakLabels,
        groupingMode: "byPhase",
        sampleIdFor: () => null,
        displayLabelByMemberId: new Map(members.map((m) => [m.id, memberRowLabel(m)])),
        xType: scale === "log" ? "log" : "linear",
        offsetScale: offset,
        preset: "clean",
      }),
    [members, rows, tracesById, effectiveTitle, showPeakTicks, showPeakLabels, scale, offset],
  );
  const fx = useFigureExport(exportSpec, effectiveTitle, "series figure");

  // ── Traces slot — read MemberList, or the editable recipe in draft ──────
  const tracesSlot = liveDraft ? (
    <RecipeEditor
      recipe={liveDraft.recipe}
      sampleNameById={sampleNameMap(corpus)}
      onRemove={onRemoveSample}
      onReorder={onReorderSample}
    />
  ) : (
    <MemberList
      members={memberData}
      {...(hoveredKey !== undefined ? { hoveredKey } : {})}
      onHoverMember={onHoverRow}
    />
  );

  // Add-sample: a native select of corpus samples not already in the recipe.
  // `addableSamples` only reads `sample_id`; project the draft recipe rows up to
  // the `SeriesSample` shape it expects (the local recipe row is a subset).
  const recipeAsSamples: api.SeriesSample[] = (liveDraft ? liveDraft.recipe : series.samples).map(
    (r) => ({
      id: r.id,
      series_id: series.id,
      sample_id: r.sample_id,
      position: r.position,
      pinned: r.pinned,
      excluded: r.excluded,
    }),
  );
  const addable = addableSamples(corpus, recipeAsSamples);

  // Full-bleed layout MIRRORING FocusPage: a [work 1fr · rail 336px] grid with
  // NO outer max-width and NO `items-start`. The grid fills the scroll
  // container's height (min-h-full) so the rail — a direct grid child — stretches
  // flush from under the header to the bottom edge (default align-items:stretch),
  // matching the Focus AssignmentRail. The work column centres its plate at the
  // builder width (1180px); the rail (bg-paper-sunk + left hairline) pins right.
  return (
    <div
      data-testid="builder-workspace"
      className="grid min-h-full grid-cols-1 lg:grid-cols-[minmax(0,1fr)_336px]"
    >
      {/* work column — full-bleed; inner content capped at the builder width */}
      <div className="min-w-0 px-6 py-6">
        <div className="mx-auto w-full max-w-[1180px]">
          <SeriesPlate
            kicker="Series"
            title={
              <PlateTitle value={effectiveTitle} onChange={onEditTitle} />
            }
            // Honest state: the plate keeps rendering the COMMITTED members
            // mid-draft (lazy-draft render model) while the title above it
            // updates live — say so, precisely: only membership/order lag.
            {...(liveDraft
              ? {
                  subtitle:
                    "Editing the recipe. Membership and order changes appear here after you confirm.",
                }
              : {})}
            rows={rows}
            offsetScale={offset}
            scale={scale}
            onScaleChange={onScaleChange}
            // The annotation toggles drive the PLATE through the same two
            // Zustand flags that feed the export spec below — one source, so
            // the on-screen figure and the export cannot diverge on these axes.
            showPeakTicks={showPeakTicks}
            showPeakLabels={showPeakLabels}
            {...(hoveredKey !== undefined ? { hoveredKey } : {})}
            onHoverRow={onHoverRow}
            legendPhases={legendPhases}
            footHint={
              <AnnotationToggleRow
                enabled={seriesHasPeaks}
                showPeakTicks={showPeakTicks}
                showPeakLabels={showPeakLabels}
                onToggleTicks={onToggleTicks}
                onToggleLabels={onToggleLabels}
              />
            }
            footNote={`offset ×${offset.toFixed(2)} · ${scale === "log" ? "log" : "linear"} q`}
            // WYSIWYG-honest gate: disabled when the exported figure would
            // contain no data. toWaterfallRows pads missing traces with
            // EMPTY_TRACE (rows.length === members.length always), so the
            // predicate must look INSIDE the rows — every-trace-empty covers
            // no members, traces still loading, a failed trace fetch, and a
            // legitimately all-empty series alike.
            actions={
              <ExportButton
                {...fx}
                ariaContext="series figure"
                disabled={!rows.some((r) => r.trace.q.length > 0)}
              />
            }
          />
        </div>
      </div>
      <BuilderRail
        grouping={groupingSummary(series)}
        {...(liveDraft && !confirmBusy && confirmReady ? { onConfirm } : {})}
        confirmBusy={confirmBusy}
        {...(liveDraft ? {} : { onAdjust: ensureDraft })}
        // copy-doesn't-lie: the rail's default WYSIWYG caption is false
        // mid-draft, so a live draft swaps in the honest variant. Precision:
        // only membership/order on the plate are stale — the draft title and
        // the offset slider still render live.
        {...(liveDraft
          ? {
              caption:
                "Membership and order on the plate are from the last confirmed figure. Confirm the series to publish your edits.",
            }
          : {})}
        orderedBy={series.ordering_variable ?? "—"}
        offset={offset}
        onOffsetChange={onOffsetChange}
        traces={
            <div className="flex flex-col gap-2">
              {chainError != null && (
                <div role="alert" className="text-caption text-error">
                  Couldn't confirm the series. Try again.
                </div>
              )}
              {/* Truthful disabled-Confirm explanation: the resolution source
                  failed to load, so the recipe cannot be turned into a plate
                  this session. Only shown mid-draft, where Confirm matters. */}
              {resolverError && liveDraft != null && (
                <div role="alert" className="text-caption text-error">
                  Couldn't load exposure data, so the series can't be confirmed. Reload the page to retry.
                </div>
              )}
              {/* While the resolution source is still loading, the gated
                  Confirm states its reason instead of sitting mute. */}
              {resolverLoading && liveDraft != null && (
                <div className="text-caption text-ink-soft">
                  Loading exposure data…
                </div>
              )}
              {tracesSlot}
              {addable.length > 0 && (
                <AddSampleSelect options={addable} onAdd={onAddSample} />
              )}
              {liveDraft && (
                <div className="flex items-center gap-2 pt-1">
                  <Button variant="ghost" onClick={onCancel}>
                    Cancel
                  </Button>
                </div>
              )}
            </div>
          }
      />
    </div>
  );
}

// ── Annotation toggles — local view-prefs (Zustand-backed, never a draft) ──
// Inlined from the carried AnnotationToggles using ui/Button (the design guard
// forbids src/print/** importing src/components/**); the Zustand flags + the
// export-spec read are identical.
function AnnotationToggleRow({
  enabled,
  showPeakTicks,
  showPeakLabels,
  onToggleTicks,
  onToggleLabels,
}: {
  /** Whether the series has indexed peaks to annotate. When false the toggles
   *  are disabled and never read as armed (nothing to show). */
  enabled: boolean;
  showPeakTicks: boolean;
  showPeakLabels: boolean;
  onToggleTicks: () => void;
  onToggleLabels: () => void;
}): JSX.Element {
  return (
    <span data-testid="annotation-toggles" role="group" aria-label="Annotation toggles" className="inline-flex items-center gap-1">
      <Button variant="ghost" armed={enabled && showPeakTicks} disabled={!enabled} onClick={onToggleTicks}>
        Peak ticks
      </Button>
      <Button variant="ghost" armed={enabled && showPeakLabels} disabled={!enabled} onClick={onToggleLabels}>
        Peak labels
      </Button>
    </span>
  );
}

// ── Plate title — inline-editable (starts a draft on first edit) ───────────
function PlateTitle({ value, onChange }: { value: string; onChange: (v: string) => void }): JSX.Element {
  return (
    <Input
      testId="builder-title-input"
      aria-label="Series title"
      value={value}
      onValueChange={onChange}
      placeholder="Untitled series"
    />
  );
}

// ── Editable recipe rows (grip + name + reorder + remove) ──────────────────
function RecipeEditor({
  recipe,
  sampleNameById,
  onRemove,
  onReorder,
}: {
  recipe: import("../../lib/series/seriesDraft").SeriesRecipeRow[];
  sampleNameById: Record<number, string>;
  onRemove: (rowId: number) => void;
  onReorder: (from: number, to: number) => void;
}): JSX.Element {
  // Rail-mirrors-plate (BU-INVERT): the waterfall paints recipe position 0 at
  // the BOTTOM of the stack, so the editor renders the recipe REVERSED (top
  // row = plate top). ONE visual-to-recipe index mapping feeds every reorder
  // path (drag drops AND the arrow buttons) so "Move down" always moves the
  // trace down IN THE FIGURE. Mapping both endpoints is sufficient because a
  // splice-move commutes with reversal (remove-at + insert-at land on the
  // mirrored positions). Removal stays id-based and needs no mapping.
  const visual = [...recipe].reverse();
  const toRecipe = (v: number): number => recipe.length - 1 - v;
  const reorderVisual = (fromV: number, toV: number): void =>
    onReorder(toRecipe(fromV), toRecipe(toV));
  const { dragItemProps, dropEdge } = useDragReorder(reorderVisual);
  return (
    <div className="flex flex-col gap-0.5" data-testid="builder-recipe">
      {visual.map((row, i) => {
        const view = recipeRowView(
          { id: row.id, series_id: 0, sample_id: row.sample_id, position: row.position, pinned: row.pinned, excluded: row.excluded },
          sampleNameById,
        );
        const dprops = dragItemProps(i);
        const edge = dropEdge(i);
        return (
          <div
            key={row.id}
            data-testid="builder-recipe-row"
            {...dprops}
            className={`group relative cursor-grab flex items-center gap-2 px-2 py-1.5 rounded border border-transparent hover:bg-plate${dprops["data-dragging"] ? " opacity-50" : ""}`}
          >
            {edge && (
              <span
                aria-hidden="true"
                className={`pointer-events-none absolute left-0 right-0 z-10 h-0.5 rounded-full bg-accent ${edge === "top" ? "-top-px" : "-bottom-px"}`}
              />
            )}
            <GripHandle />
            <span className="flex-1 min-w-0 truncate text-meta text-ink">{view.name}</span>
            <IconButton
              label="Move up"
              tone="ghost"
              data-testid="builder-recipe-up"
              disabled={i === 0}
              onClick={() => reorderVisual(i, i - 1)}
            >
              &#9650;
            </IconButton>
            <IconButton
              label="Move down"
              tone="ghost"
              data-testid="builder-recipe-down"
              disabled={i === visual.length - 1}
              onClick={() => reorderVisual(i, i + 1)}
            >
              &#9660;
            </IconButton>
            <IconButton
              label="Remove sample"
              tone="ghost"
              dismiss
              data-testid="builder-recipe-remove"
              onClick={() => onRemove(row.id)}
            />
          </div>
        );
      })}
    </div>
  );
}

// ── Add-sample native select (addable corpus samples) ──────────────────────
function AddSampleSelect({
  options,
  onAdd,
}: {
  options: api.CorpusSample[];
  onAdd: (sampleId: number) => void;
}): JSX.Element {
  return (
    <select
      data-testid="builder-add-sample-select"
      aria-label="Add a sample to the series"
      className="w-full border border-hair-strong bg-plate rounded px-3 py-2 text-meta text-ink"
      value=""
      onChange={(e) => {
        const v = Number(e.target.value);
        if (Number.isFinite(v) && v > 0) onAdd(v);
      }}
    >
      <option value="">+ Add sample…</option>
      {options.map((s) => (
        <option key={s.id} value={s.id}>
          {s.display_name ?? s.name}
        </option>
      ))}
    </select>
  );
}

function sampleNameMap(corpus: api.CorpusSample[]): Record<number, string> {
  const out: Record<number, string> = {};
  for (const s of corpus) out[s.id] = s.display_name ?? s.name ?? `Sample ${s.id}`;
  return out;
}
