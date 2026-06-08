import { useEffect, useRef, useState } from "react";
import { useParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { SeriesPlate } from "../components/SeriesPlate";
import { BuilderRail } from "../components/BuilderRail";
import { MemberList } from "../components/MemberList";
import { IconButton, Button, EmptyState, Input } from "../ui";
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
import { buildExportPng, canExportPng } from "../../lib/figure-export/renderer";
import { canCopyPngToClipboard, copyPngToClipboard } from "../../lib/figure-export/clipboard";
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

  // ── Confirm chain (Save → Commit → discard) ─────────────────────────────
  const save = useSaveSeries();
  const commit = useCommitSeriesPlate();
  const stage = useRef<"idle" | "saving" | "committing">("idle");

  const series = seriesQ.data;
  // The active draft only counts when it targets THIS series.
  const liveDraft = draft !== null && series !== undefined && draft.id === series.id ? draft : null;

  // LOAD-BEARING: the page never calls save.reset()/commit.reset(). A lingering
  // save.isSuccess/commit.isSuccess === true between Confirm runs is INERT
  // because every chain effect below is `stage`-ref gated — they fire only while
  // stage.current is the matching phase ("saving"/"committing"), and onConfirm
  // re-arms the chain by flipping stage to "saving" itself. Do NOT add a reset()
  // or remove the stage gating: a reset would race the SSE-vs-HTTP deferred
  // resolution, and dropping the gate would let a stale isSuccess re-fire the
  // commit on the next unrelated re-render.
  //
  // Save landed → commit the FRESH plate from the save response, never the
  // stale draft. `save.data` is the updated full Series on the HTTP-wins path
  // (useSaveSeries surfaces mutation.data); guard with isFullSeries so the
  // SSE-wins partial shape (no members) can't be committed.
  useEffect(() => {
    if (stage.current !== "saving" || !save.isSuccess || !save.data) return;
    if (!api.isFullSeries(save.data)) {
      // SSE-wins partial: cannot commit from a memberless shape. Reset so the
      // user can retry; the save itself already landed.
      stage.current = "idle";
      return;
    }
    stage.current = "committing";
    commit.mutate({ id, ...buildSeriesCommitBody(save.data.members) });
  }, [save.isSuccess, save.data, id, commit]);

  // Commit landed → drop the draft; stay (now read state).
  useEffect(() => {
    if (stage.current !== "committing" || !commit.isSuccess) return;
    stage.current = "idle";
    discardDraft();
  }, [commit.isSuccess, discardDraft]);

  // Either error → reset so the user can retry.
  useEffect(() => {
    if (save.error || commit.error) stage.current = "idle";
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
      fallback={<div className="p-8 text-sm text-ink-faint">Loading series…</div>}
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

  const effectiveTitle = liveDraft ? liveDraft.title : series.title;

  // ── Figure export (Copy as PNG) ─────────────────────────────────────────
  const exportSpec = (): import("../../lib/figure-export/types").ExportSpec =>
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
    });
  const onCopyPng = (): void => {
    if (!canCopyPngToClipboard() || !canExportPng()) {
      showToast("This browser can't copy figures. Try Save instead.", "error");
      return;
    }
    void (async (): Promise<void> => {
      try {
        const blob = await buildExportPng(exportSpec());
        await copyPngToClipboard(blob);
        showToast("Copied figure to clipboard", "success");
      } catch {
        showToast("Couldn't copy figure.", "error");
      }
    })();
  };

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

  return (
    <PageFrame width="builder" className="px-6 py-6">
      <div className="grid grid-cols-[1fr_336px] gap-0 items-start">
        <SeriesPlate
          kicker="Series"
          title={
            <PlateTitle value={effectiveTitle} onChange={onEditTitle} />
          }
          rows={rows}
          offsetScale={offset}
          scale={scale}
          onScaleChange={onScaleChange}
          {...(hoveredKey !== undefined ? { hoveredKey } : {})}
          onHoverRow={onHoverRow}
          legendPhases={legendPhases}
          footHint={
            <AnnotationToggleRow
              showPeakTicks={showPeakTicks}
              showPeakLabels={showPeakLabels}
              onToggleTicks={onToggleTicks}
              onToggleLabels={onToggleLabels}
            />
          }
          footNote={`offset ×${offset.toFixed(2)} · ${scale === "log" ? "log" : "linear"} q`}
          className="m-4"
        />
        <BuilderRail
          grouping={groupingSummary(series)}
          {...(liveDraft && !confirmBusy ? { onConfirm } : {})}
          onAdjust={ensureDraft}
          orderedBy={series.ordering_variable ?? "—"}
          offset={offset}
          onOffsetChange={onOffsetChange}
          scale={scale}
          onScaleChange={onScaleChange}
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
          onCopyPng={onCopyPng}
        />
      </div>
    </PageFrame>
  );
}

// ── Annotation toggles — local view-prefs (Zustand-backed, never a draft) ──
// Inlined from the carried AnnotationToggles using ui/Button (the design guard
// forbids src/print/** importing src/components/**); the Zustand flags + the
// export-spec read are identical.
function AnnotationToggleRow({
  showPeakTicks,
  showPeakLabels,
  onToggleTicks,
  onToggleLabels,
}: {
  showPeakTicks: boolean;
  showPeakLabels: boolean;
  onToggleTicks: () => void;
  onToggleLabels: () => void;
}): JSX.Element {
  return (
    <span data-testid="annotation-toggles" role="group" aria-label="Annotation toggles" className="inline-flex items-center gap-1">
      <Button variant="ghost" armed={showPeakTicks} onClick={onToggleTicks}>
        Peak ticks
      </Button>
      <Button variant="ghost" armed={showPeakLabels} onClick={onToggleLabels}>
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
  return (
    <div className="flex flex-col gap-0.5" data-testid="builder-recipe">
      {recipe.map((row, i) => {
        const view = recipeRowView(
          { id: row.id, series_id: 0, sample_id: row.sample_id, position: row.position, pinned: row.pinned, excluded: row.excluded },
          sampleNameById,
        );
        return (
          <div
            key={row.id}
            data-testid="builder-recipe-row"
            className="flex items-center gap-2 px-2 py-1.5 rounded border border-transparent hover:bg-plate"
          >
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
