import { useCallback, useEffect, useRef, useState } from "react";
import { useParams } from "react-router-dom";
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
  useSaveSeries,
  useCommitSeriesPlate,
} from "../../queries";
import { useAppState } from "../../state";
import * as api from "../../api";
import type { Series, SeriesMember } from "../../api";
import { toWaterfallRows } from "../waterfall/waterfallModel";
import {
  membersToMemberData,
  groupingSummary,
  legendPhasesOf,
  addableSamples,
  recipeRowView,
} from "./builderAdapters";
import { buildSeriesSaveBody } from "../../lib/series/buildSeriesSaveBody";
import { buildSeriesCommitBody } from "../../lib/series/buildSeriesCommitBody";
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
  const id = Number(useParams<{ id: string }>().id);
  const seriesQ = useSeries(Number.isFinite(id) ? id : undefined);
  const tracesQ = useSeriesTraces(Number.isFinite(id) ? id : undefined);
  const corpusQ = useCorpusSamples();

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
  // Save landed → commit the FRESH plate from the SERIES QUERY CACHE, never the
  // stale draft and never `save.data`. Under the queue's SSE-wins race the
  // own-op SSE confirmation resolves `save.data` as the synthesized PARTIAL
  // ({id}, no members/state — saveSeriesMutator.synthesizeFromSse), so gating
  // the commit on isFullSeries(save.data) stalls the chain on live backends
  // (P0 BU-CONFIRMSTALL). The cache is the one source that is correct on BOTH
  // race outcomes: the mutator's onSuccess either setQueryData's the full HTTP
  // response or invalidates the key so TanStack refetches the canonical
  // projection — either way the series key gets a full Series and its
  // dataUpdatedAt advances past the Confirm-time watermark. The cache is also
  // the only holder of server-resolved SeriesMember[] (exposure_ids included)
  // that buildSeriesCommitBody builds from — the local draft is a sample-level
  // recipe, not resolved members (the body itself is positional; member ids
  // are stripped on the wire).
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
    commit.mutate({ id, ...buildSeriesCommitBody(fresh.members) });
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
    showToast("Series confirmed", "success");
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
  if (seriesQ.isError) {
    return (
      <PageFrame width="builder" className="px-8 py-8">
        <EmptyState
          title="Couldn't load this series"
          body="It may have been deleted, or the request failed. Try reloading the page."
        />
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
    if (!liveDraft || stage.current !== "idle") return;
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

  const confirmBusy = save.isPending || commit.isPending || stage.current !== "idle";
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
  chainError,
}: BuilderBodyProps): JSX.Element {
  // Render model: the figure plate ALWAYS shows the committed plate (members);
  // a draft edits the recipe (membership/title), which only re-resolves the
  // plate after a Save→Commit round-trip. So the waterfall/MemberList read the
  // committed members; the draft drives only the title + the editable rail.
  const members: SeriesMember[] = series.members;
  const rows = toWaterfallRows(members, tracesById);
  const memberData = membersToMemberData(members);
  const legendPhases = legendPhasesOf(members);

  // Peak-annotation gating: ticks/labels annotate indexed-peak anchors. A
  // form-factor / unindexed series has no anchors on any row → nothing to
  // annotate, so the annotation toggles are inert (not armed, disabled).
  // controls-don't-lie: don't arm a toggle that has nothing to show.
  const seriesHasPeaks = rows.some((r) => r.anchors.length > 0);

  const effectiveTitle = liveDraft ? liveDraft.title : series.title;

  // ── Figure export (plate-head ExportButton, wired via useFigureExport) ───
  const exportSpec = useCallback(
    (): import("../../lib/figure-export/types").ExportSpec =>
      buildMultiTraceExportSpec({
        members,
        traces: new Map(Object.entries(tracesById).map(([k, v]) => [Number(k), v])),
        comparisonTitle: effectiveTitle,
        xDomain: null,
        showPeakTicks,
        showPeakLabels,
        groupingMode: "bySample",
        sampleIdFor: () => null,
        preset: "clean",
      }),
    [members, tracesById, effectiveTitle, showPeakTicks, showPeakLabels],
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
        {...(liveDraft && !confirmBusy ? { onConfirm } : {})}
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
  const { dragItemProps, dropEdge } = useDragReorder(onReorder);
  return (
    <div className="flex flex-col gap-0.5" data-testid="builder-recipe">
      {recipe.map((row, i) => {
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
              onClick={() => onReorder(i, i - 1)}
            >
              &#9650;
            </IconButton>
            <IconButton
              label="Move down"
              tone="ghost"
              data-testid="builder-recipe-down"
              disabled={i === recipe.length - 1}
              onClick={() => onReorder(i, i + 1)}
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
