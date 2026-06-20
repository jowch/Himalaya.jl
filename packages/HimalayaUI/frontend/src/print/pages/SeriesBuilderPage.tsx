import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { PageFrame } from "../components/PageFrame";
import { SeriesPlate } from "../components/SeriesPlate";
import { BuilderRail } from "../components/BuilderRail";
import { MemberList } from "../components/MemberList";
import { IconButton, Button, EmptyState, Input, GripHandle } from "../ui";
import { useDragReorder } from "../components/useDragReorder";
import type { DragItemProps } from "../components/useDragReorder";
import { useReorderShortcuts } from "../shell/useReorderShortcuts";
import { useShortcuts } from "../shell/useShortcuts";
import { shortcutLabel } from "../shell/shortcuts";
import {
  useSeries,
  useSeriesTraces,
  useCorpusSamples,
  useCorpusPickerSamples,
  useExperiments,
  useSaveSeries,
  useCommitSeriesPlate,
} from "../../queries";
import { AddSamplePicker } from "../components/AddSamplePicker";
import { useAppState } from "../../state";
import { useDocumentTitle } from "../../hooks/useDocumentTitle";
import * as api from "../../api";
import type { Series, SeriesMember } from "../../api";
import { toWaterfallRows, waterfallQDomain } from "../waterfall/waterfallModel";
import { memberRowLabel } from "../../lib/series/memberRead";
import {
  membersToMemberData,
  legendPhasesOf,
  addableSamples,
  recipeRowView,
} from "./builderAdapters";
import { buildSeriesSaveBody } from "../../lib/series/buildSeriesSaveBody";
import { isSeriesDraftDirty } from "../../lib/series/isSeriesDraftDirty";
import { buildPlateFromRecipe } from "../../lib/series/buildPlateFromRecipe";
import { buildCleanFigureSvg, type FigureTraceKey } from "../export/cleanFigureSvg";
import { buildSeriesFigureKeys } from "../export/seriesFigureKeys";
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
  // Experiments power the add-sample picker's group headers (BU-PICKER).
  const experimentsQ = useExperiments();
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

  // BU-RESCORE3 N3: name the browser tab after the series (was the static
  // "Himalaya"). Committed title — nullish while loading falls back to the app
  // name. Hooks-safe here: above every early return.
  useDocumentTitle(seriesQ.data?.title ?? null);

  // ── Draft (lazy) + view-pref state ──────────────────────────────────────
  const draft = useAppState((s) => s.seriesDraft);
  const startDraft = useAppState((s) => s.startSeriesDraftFromSeries);
  const discardDraft = useAppState((s) => s.discardSeriesDraft);
  const setTitle = useAppState((s) => s.setSeriesDraftTitle);
  const addSample = useAppState((s) => s.addSeriesSample);
  const removeSample = useAppState((s) => s.removeSeriesSample);
  const reorderSample = useAppState((s) => s.reorderSeriesSample);
  const restoreDraft = useAppState((s) => s.restoreSeriesDraft);
  const showPeakTicks = useAppState((s) => s.showPeakTicks);
  const showPeakLabels = useAppState((s) => s.showPeakLabels);
  const setShowPeakTicks = useAppState((s) => s.setShowPeakTicks);
  const setShowPeakLabels = useAppState((s) => s.setShowPeakLabels);

  // ── Local view state (never starts a draft) ─────────────────────────────
  const [offset, setOffset] = useState(1.2);
  const [scale, setScale] = useState<Scale>("log");
  const [hoveredKey, setHoveredKey] = useState<string | undefined>(undefined);

  // ── Draft undo / redo (BU-UNDO, locked decision #3) ─────────────────────
  // Snapshot-based history of the structural recipe edits (reorder / add /
  // remove). The free-text title keeps the field's own native undo (⌘Z there is
  // suppressed by suppressGlobalKeys), so this stack only carries the edits a
  // field cannot recover. `record()` captures the draft slot BEFORE a mutation;
  // undoing the very first edit restores `null` (back to read mode). Updaters
  // stay pure (the snapshot is precomputed) so StrictMode's double-invoke can't
  // double-push. Declared in the top hooks block — above every early return.
  type DraftSlot = ReturnType<typeof useAppState.getState>["seriesDraft"];
  const [undoPast, setUndoPast] = useState<DraftSlot[]>([]);
  const [undoFuture, setUndoFuture] = useState<DraftSlot[]>([]);
  const record = (): void => {
    const snap = useAppState.getState().seriesDraft;
    setUndoPast((p) => [...p, snap]);
    setUndoFuture([]);
  };
  const undo = (): void => {
    if (undoPast.length === 0) return;
    const prev = undoPast[undoPast.length - 1]!;
    const cur = useAppState.getState().seriesDraft;
    restoreDraft(prev);
    setUndoPast(undoPast.slice(0, -1));
    setUndoFuture([cur, ...undoFuture]);
  };
  const redo = (): void => {
    if (undoFuture.length === 0) return;
    const next = undoFuture[0]!;
    const cur = useAppState.getState().seriesDraft;
    restoreDraft(next);
    setUndoFuture(undoFuture.slice(1));
    setUndoPast([...undoPast, cur]);
  };
  useShortcuts({ undo, redo });

  // ── Confirm chain (Save → await fresh cache → Commit → discard) ─────────
  const save = useSaveSeries();
  const commit = useCommitSeriesPlate();
  // BU-PROGRESS (F-ERRSILENT): name the FAILED STEP so a chain failure isn't a
  // generic dead-end. A save failure means nothing was published; a commit
  // failure means the order saved but the figure didn't rebuild — different
  // recovery, different words. Shared by the inline alert + the assertive toast.
  const chainErrorMessage = save.error
    ? "Couldn't save your changes. Check your connection and try Confirm again."
    : commit.error
      ? "Your order was saved, but publishing the figure failed. Try Confirm again."
      : null;
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

  // BU-NAMES: one shared identity token across the builder's member views. The
  // figure plate + reading list label a trace by its exposure (memberRowLabel:
  // label_override, else "exp <id>"); the editable recipe rail labels by SAMPLE
  // name. So in edit mode the recipe ("LL2") and the figure ("exp 42") can't be
  // cross-referenced. Resolve each recipe sample to its committed member's
  // figure label (sample → indexing exposure → member) and surface it as a muted
  // suffix on the row. A draft-only sample not yet in the committed plate has no
  // figure trace, so it gets no token (the lookup misses).
  const figureLabelBySample = useMemo(() => {
    const out = new Map<number, string>();
    const labelByExposure = new Map(
      (series?.members ?? []).map((m) => [m.exposure_id, memberRowLabel(m)]),
    );
    for (const [sampleId, exposureId] of exposureBySample) {
      if (exposureId == null) continue;
      const lab = labelByExposure.get(exposureId);
      if (lab) out.set(sampleId, lab);
    }
    return out;
  }, [series?.members, exposureBySample]);

  // Export-figure per-trace key entries (BU-NAMES figure): the figure labels each
  // trace by SAMPLE name (not "exp <id>") and carries its phase·a·κ. Resolve
  // exposure → sample name (every exposure a sample owns maps to its name, so a
  // member on a non-indexing exposure still resolves) and the q-units (Å/nm, for
  // the lattice + κ units) from the first resolvable member's experiment.
  const figureTraceKeys = useMemo<FigureTraceKey[]>(() => {
    const members = series?.members ?? [];
    const nameByExposure = new Map<number, string>();
    const expByExposure = new Map<number, number>();
    for (const r of pickerQ.data ?? []) {
      const nm = r.sample.name ?? `Sample ${r.sample.id}`;
      if (r.indexing_exposure_id != null) {
        nameByExposure.set(r.indexing_exposure_id, nm);
        expByExposure.set(r.indexing_exposure_id, r.sample.experiment_id);
      }
      for (const e of r.all_exposures ?? []) {
        nameByExposure.set(e.id, nm);
        expByExposure.set(e.id, r.sample.experiment_id);
      }
    }
    const qByExp = new Map((experimentsQ.data ?? []).map((e) => [e.id, e.q_units]));
    let qUnits: string | null = null;
    for (const mem of members) {
      const expId = mem.exposure_id != null ? expByExposure.get(mem.exposure_id) : undefined;
      if (expId != null && qByExp.has(expId)) { qUnits = qByExp.get(expId)!; break; }
    }
    return buildSeriesFigureKeys(members, {
      sampleNameForExposure: (id) => (id != null ? nameByExposure.get(id) ?? null : null),
      fallbackLabel: memberRowLabel,
      qUnits,
    });
  }, [series?.members, pickerQ.data, experimentsQ.data]);

  // BU-NAVAWAY-DRAFT: warn before the tab is closed / reloaded while a draft
  // carries UNSAVED CHANGES. In-app route changes keep the draft alive (it
  // lives in the Zustand store + sessionStorage, recoverable on return), so
  // they need no guard; tab close is the one path that discards the per-tab
  // sessionStorage draft for good. The guard arms only when the draft actually
  // differs from the saved series (controls-don't-lie), never on a pristine
  // fork. `useBlocker` would cover in-app nav too but requires a data router;
  // the app uses <BrowserRouter>, so beforeunload is the available, honest lever.
  const draftDirty = liveDraft !== null && series !== undefined && isSeriesDraftDirty(liveDraft, series);
  useEffect(() => {
    if (!draftDirty) return;
    const handler = (e: BeforeUnloadEvent): void => {
      e.preventDefault();
      e.returnValue = ""; // legacy Chrome/Firefox: a non-empty value arms the prompt
    };
    window.addEventListener("beforeunload", handler);
    return () => window.removeEventListener("beforeunload", handler);
  }, [draftDirty]);

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
    // BU-MODESWITCH-LEAK: same as Cancel — the committed draft is gone, so its
    // undo/redo history must not survive into the next edit of this series.
    setUndoPast([]);
    setUndoFuture([]);
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
      // a toast so the failure reaches a user whose focus is elsewhere. Same
      // cause-named copy (BU-PROGRESS) so the two channels never disagree.
      showToast(
        chainErrorMessage ?? "Couldn't confirm the series. Try again.",
        "error",
      );
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
            as="h1"
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
    record();
    ensureDraft();
    addSample(sampleId);
  };
  // Guard-enforce the draft invariant rather than relying on the recipe editor
  // only rendering under a live draft (layout-enforced). `ensureDraft()` is a
  // no-op when a draft already exists.
  const onRemoveSample = (rowId: number): void => {
    record();
    ensureDraft();
    removeSample(rowId);
  };
  const onReorderSample = (from: number, to: number): void => {
    record();
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
    // BU-MODESWITCH-LEAK: the undo/redo history is scoped to the draft we are
    // discarding. Drop it so the NEXT draft on this series doesn't inherit a
    // stale stack (a stray ⌘Z/⌘⇧Z restoring an edit from the cancelled session).
    setUndoPast([]);
    setUndoFuture([]);
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
  // BU-PROGRESS: the Save→await→Commit chain is multi-step, but the user only
  // saw a single static "Saving…". Name the CURRENT step from the same reactive
  // mutation states so the label advances visibly across the chain (a stuck
  // process looks different from a working one). The awaiting-fresh gap (cache
  // refetch between save success and commit) reads "Confirming…".
  const confirmProgressLabel = !confirmBusy
    ? undefined
    : save.isPending
      ? "Saving order…"
      : commit.isPending
        ? "Publishing the figure…"
        : "Confirming…";

  return (
    <Skeleton
      name="series-builder"
      // Full-height flex column so the builder grid stretches to the bottom of the
      // page even when its content is shorter than the viewport. `min-h-full`
      // fills `main`; `flex flex-col` makes this a flex container — boneyard wraps
      // children in `display:contents`, which UNWRAPS in a flex parent, so the grid
      // becomes a direct flex item and can `grow` to fill (a plain `min-h-full` on
      // the grid can't, because % height won't resolve through a contents box).
      className="flex flex-col min-h-full"
      loading={seriesQ.isLoading || series === undefined}
      fixture={BUILDER_FIXTURE}
      fallback={<div className="p-8 text-sm text-ink-soft">Loading series…</div>}
    >
      {series && (
        <BuilderBody
          series={series}
          tracesById={tracesQ.data ?? {}}
          tracesLoading={tracesQ.isLoading}
          tracesError={tracesQ.isError}
          corpus={corpusQ.data ?? []}
          experiments={experimentsQ.data ?? []}
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
          onUndo={undo}
          canUndo={undoPast.length > 0}
          confirmBusy={confirmBusy}
          confirmLabel={confirmProgressLabel}
          confirmReady={resolverReady}
          resolverError={pickerQ.isError}
          resolverLoading={!resolverReady && !pickerQ.isError}
          chainErrorMessage={chainErrorMessage}
          figureLabelBySample={figureLabelBySample}
          figureTraceKeys={figureTraceKeys}
        />
      )}
    </Skeleton>
  );
}

interface BuilderBodyProps {
  series: Series;
  tracesById: Record<number, api.Trace>;
  /** Trace fetch state — lets a disabled Export name its reason (BU-EXPORT-EMPTYALLOWED). */
  tracesLoading: boolean;
  tracesError: boolean;
  corpus: api.CorpusSample[];
  experiments: api.Experiment[];
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
  /** Step the last structural recipe edit back (BU-UNDO). */
  onUndo: () => void;
  /** Whether there is a recorded edit to undo (drives the visible control). */
  canUndo: boolean;
  confirmBusy: boolean;
  /** BU-PROGRESS: the step-named busy label for the Save action while the
   *  Save→Commit chain is in flight ("Saving order…" / "Publishing the figure…"
   *  / "Confirming…"). Undefined when not busy. */
  confirmLabel: string | undefined;
  /** Picker projection loaded — the recipe→plate resolution source is ready.
   *  Confirm stays disabled until then (an unresolvable recipe would commit
   *  an empty plate). */
  confirmReady: boolean;
  /** Picker projection FAILED — Confirm can never become ready this session. */
  resolverError: boolean;
  resolverLoading: boolean;
  /** BU-PROGRESS: cause-named chain-failure copy, or null when no failure.
   *  Distinguishes a save-step from a commit-step failure. */
  chainErrorMessage: string | null;
  /** BU-NAMES: recipe sample id → the figure trace's label (memberRowLabel),
   *  shown as a muted suffix so a recipe row cross-references its figure trace. */
  figureLabelBySample: Map<number, string>;
  /** Per-trace export-figure key entries (sample name + phase·a·κ + categorical
   *  colour), parallel to the committed members. Built in the outer component
   *  (which holds the picker + experiment data) and threaded into renderSvg. */
  figureTraceKeys: FigureTraceKey[];
}

function BuilderBody({
  series,
  tracesById,
  tracesLoading,
  tracesError,
  corpus,
  experiments,
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
  onUndo,
  canUndo,
  confirmBusy,
  confirmLabel,
  confirmReady,
  resolverError,
  resolverLoading,
  chainErrorMessage,
  figureLabelBySample,
  figureTraceKeys,
}: BuilderBodyProps): JSX.Element {
  // Render model: the figure plate ALWAYS shows the committed plate (members);
  // a draft edits the recipe (membership/title), which only re-resolves the
  // plate after a Save→Commit round-trip. So the waterfall/MemberList read the
  // committed members; the draft drives only the title + the editable rail.
  const members: SeriesMember[] = series.members;
  const rows = toWaterfallRows(members, tracesById);
  // BU-EXPORT-EMPTYALLOWED: Export is WYSIWYG-gated — disabled when no trace
  // carries data — but a bare-disabled button can't say WHY. Distinguish the
  // three causes the predicate folds together so the disabled state names its
  // reason (no-data vs loading vs error) instead of going silently dead.
  const exportDisabled = !rows.some((r) => r.trace.q.length > 0);
  const exportDisabledReason = !exportDisabled
    ? undefined
    : tracesLoading
      ? "Traces are still loading."
      : tracesError
        ? "Traces couldn't load, so there's nothing to export."
        : "This series has no trace data to export.";
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
  // WYSIWYG contract (BU-EXPORTDIVERGE): the export renders the SAME data the
  // plate composes, through the SAME greenfield axis math (makeAxis/axisTicks) —
  // see cleanFigureSvg. It reads the plate's own sources: the waterfall `rows`,
  // the padded q-domain, the log/linear toggle, the offset slider, the
  // annotation flags. The export DIVERGES on identity by design: distinguishable
  // per-trace colours + a sample/phase/a/κ key block (traceKeys), vs the plate's
  // phase-coloured traces. Only the SKIN + that key differ.
  const renderSvg = useCallback(
    () =>
      buildCleanFigureSvg({
        rows,
        traceKeys: figureTraceKeys,
        title: effectiveTitle,
        footer: "q normalized · intensity offset for clarity",
        xType: scale === "log" ? "log" : "linear",
        offsetScale: offset,
        showPeakTicks,
        showPeakLabels,
        qDomain: waterfallQDomain(rows),
      }),
    [rows, figureTraceKeys, effectiveTitle, scale, offset, showPeakTicks, showPeakLabels],
  );
  // Descriptive, product-tagged stem (buildFilename slugifies it): e.g.
  // "himalaya-ll37-titration-2026-06-13.svg". The series title is itself the
  // descriptor (the default is "Series by <var>"), so no redundant "series-".
  const fx = useFigureExport(renderSvg, `himalaya-${effectiveTitle}`, "series figure");

  // ── Traces slot — read MemberList, or the editable recipe in draft ──────
  const tracesSlot = liveDraft ? (
    <RecipeEditor
      recipe={liveDraft.recipe}
      sampleNameById={sampleNameMap(corpus)}
      figureLabelBySample={figureLabelBySample}
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
      // `grow` (flex-grow): as a flex item of the Skeleton's flex-col wrapper this
      // fills the page height so the rail reaches the bottom; `min-h-full` is the
      // fallback when no flex parent is present (stories/tests). It is still a grid
      // for its own [work · rail] columns.
      className="grid grow min-h-full grid-cols-1 lg:grid-cols-[minmax(0,1fr)_336px]"
    >
      {/* work column — full-bleed; inner content capped at the builder width */}
      <div className="min-w-0 px-6 py-6">
        <div className="mx-auto w-full max-w-[1180px]">
          <SeriesPlate
            kicker="Series"
            // BU-EDITMODE: the title is only editable in edit mode. Read mode
            // shows it as a static serif heading (typing can no longer silently
            // start a draft); a live draft swaps in the inline-editable PlateTitle.
            title={
              liveDraft ? (
                <PlateTitle value={effectiveTitle} onChange={onEditTitle} />
              ) : (
                effectiveTitle || "Untitled series"
              )
            }
            // BU-NOHEAD: in edit mode the visible title is an editable Input, so
            // the real page heading is this plain text (rendered sr-only by
            // PlateHeader). In read mode the static title IS the heading, so no
            // headingText is passed (PlateHeader renders the title text as <h1>).
            {...(liveDraft ? { headingText: effectiveTitle || "Untitled series" } : {})}
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
                disabled={exportDisabled}
                {...(exportDisabledReason ? { disabledReason: exportDisabledReason } : {})}
              />
            }
          />
        </div>
      </div>
      <BuilderRail
        {...(liveDraft && !confirmBusy && confirmReady ? { onConfirm } : {})}
        confirmBusy={confirmBusy}
        {...(confirmLabel !== undefined ? { confirmLabel } : {})}
        {...(liveDraft ? { onCancel } : { onAdjust: ensureDraft })}
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
        // Honest section label: rows are only reorderable in draft mode (the
        // editable RecipeEditor); the read-mode MemberList cannot be dragged.
        reorderable={liveDraft != null}
        offset={offset}
        onOffsetChange={onOffsetChange}
        traces={
            <div className="flex flex-col gap-2">
              {chainErrorMessage != null && (
                <div role="alert" className="text-caption text-error">
                  {chainErrorMessage}
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
              {/* BU-EDITMODE: adding samples is an edit, so the picker only
                  shows under a live draft — read mode has no edit affordance, so
                  the Edit door is the sole way in (a real mode toggle). */}
              {liveDraft && addable.length > 0 && (
                <AddSamplePicker
                  options={addable}
                  experiments={experiments}
                  onAdd={onAddSample}
                />
              )}
              {/* Undo is a contextual recipe-edit affordance, kept beside the
                  rows; the draft's commit/discard pair (Save changes / Cancel)
                  are the proper buttons in the rail's compose block, not here. */}
              {liveDraft && canUndo && (
                <div className="flex items-center gap-2 pt-1">
                  <Button
                    variant="ghost"
                    onClick={onUndo}
                    data-testid="builder-undo"
                    title={`Undo last change (${shortcutLabel("undo")})`}
                  >
                    ↺ Undo
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
// BU-TITLE-FIELD: the figure's name is a TITLE, not a form field. The `title`
// variant renders it as Newsreader serif display ink with a dotted underline
// that firms on focus (confident ink, re-openable) — matching the scoping
// plate's title idiom — instead of a permanent bordered sans-serif input box.
function PlateTitle({ value, onChange }: { value: string; onChange: (v: string) => void }): JSX.Element {
  return (
    <Input
      variant="title"
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
  figureLabelBySample,
  onRemove,
  onReorder,
}: {
  recipe: import("../../lib/series/seriesDraft").SeriesRecipeRow[];
  sampleNameById: Record<number, string>;
  figureLabelBySample: Map<number, string>;
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
  // Unified Alt+↑/↓ reorder power-gesture (shared registry). Routes through the
  // focused row's own ▲▼ Move buttons so their boundary focus-retention dance
  // (BU-MOVEFOCUS) and disabled-at-extreme guard are reused verbatim.
  useReorderShortcuts({
    rowSelector: "[data-reorder-row]",
    move: (_index, delta, row) => {
      const btn = row.querySelector<HTMLButtonElement>(
        delta < 0 ? '[data-testid="builder-recipe-up"]' : '[data-testid="builder-recipe-down"]',
      );
      if (btn && !btn.disabled) btn.click();
    },
  });
  return (
    <div className="flex flex-col gap-0.5" data-testid="builder-recipe">
      {visual.map((row, i) => {
        const view = recipeRowView(
          { id: row.id, series_id: 0, sample_id: row.sample_id, position: row.position, pinned: row.pinned, excluded: row.excluded },
          sampleNameById,
        );
        const figureLabel = figureLabelBySample.get(row.sample_id);
        return (
          <RecipeRow
            key={row.id}
            name={view.name}
            {...(figureLabel !== undefined ? { figureLabel } : {})}
            index={i}
            count={visual.length}
            dragItemProps={dragItemProps(i)}
            dropEdge={dropEdge(i)}
            onMove={reorderVisual}
            onRemove={() => onRemove(row.id)}
          />
        );
      })}
    </div>
  );
}

// ── One editable recipe row — owns the focus-retention dance ────────────────
// When a keyboard user activates "Move up"/"Move down" and the move lands the
// row at an extreme, the button they pressed self-disables and can no longer
// hold focus → focus would fall to <body> (WCAG 2.4.3 failure). Mirroring the
// RepresentativeBox LO-REPDROP precedent: a ref flag is set ONLY by THIS row's
// own move activation and consumed in a useLayoutEffect (before paint, no focus
// flash). If the move was ours and the just-activated direction is now
// disabled, focus redirects to the sibling move button (still actionable). A
// re-render from a foreign cause (another row's reorder, SSE) never sets the
// flag, so it never yanks focus.
function RecipeRow({
  name,
  figureLabel,
  index,
  count,
  dragItemProps,
  dropEdge,
  onMove,
  onRemove,
}: {
  name: import("react").ReactNode;
  /** BU-NAMES: the figure trace's label for this sample (e.g. "exp 42"), shown
   *  muted after the name so the row cross-references its plate trace. */
  figureLabel?: string;
  index: number;
  count: number;
  dragItemProps: DragItemProps;
  dropEdge: "top" | "bottom" | null;
  onMove: (fromV: number, toV: number) => void;
  onRemove: () => void;
}): JSX.Element {
  const upRef = useRef<HTMLButtonElement>(null);
  const downRef = useRef<HTMLButtonElement>(null);
  // null = no pending self-move; otherwise the direction the user just pressed.
  const pendingDir = useRef<"up" | "down" | null>(null);
  const atTop = index === 0;
  const atBottom = index === count - 1;

  useLayoutEffect(() => {
    const dir = pendingDir.current;
    pendingDir.current = null; // consume — the flag lives one transition
    if (!dir) return;
    // The pressed direction self-disabled at the extreme → move focus to the
    // sibling that is still actionable. Otherwise the pressed button still
    // holds focus and nothing should move.
    if (dir === "up" && atTop) downRef.current?.focus();
    else if (dir === "down" && atBottom) upRef.current?.focus();
  }, [index, atTop, atBottom]);

  return (
    <div
      data-testid="builder-recipe-row"
      data-reorder-row
      data-reorder-index={index}
      {...dragItemProps}
      className={`group relative cursor-grab flex items-center gap-2 px-2 py-1.5 rounded border border-transparent hover:bg-plate${dragItemProps["data-dragging"] ? " opacity-50" : ""}`}
    >
      {dropEdge && (
        <span
          aria-hidden="true"
          className={`pointer-events-none absolute left-0 right-0 z-10 h-0.5 rounded-full bg-accent ${dropEdge === "top" ? "-top-px" : "-bottom-px"}`}
        />
      )}
      <GripHandle />
      <span className="flex-1 min-w-0 truncate text-meta text-ink">{name}</span>
      {figureLabel !== undefined && (
        <span
          data-testid="builder-recipe-figure-label"
          className="flex-shrink-0 text-caption text-ink-soft tabular-nums"
        >
          {figureLabel}
        </span>
      )}
      {/* BU-REORDER-CRAMP: the reorder/remove controls are ONE tight cluster
          (gap-0.5), set off from the data zone (grip / name / exp) by a left
          hairline + padding, so the ▲▼✕ trio reads as a control group rather than
          jostling the name in a flat equal-gap run. */}
      <div className="flex items-center gap-0.5 flex-shrink-0 ml-1 pl-2 border-l border-hair">
        <IconButton
          ref={upRef}
          label="Move up"
          tone="ghost"
          data-testid="builder-recipe-up"
          disabled={atTop}
          onClick={() => {
            pendingDir.current = "up";
            onMove(index, index - 1);
          }}
        >
          &#9650;
        </IconButton>
        <IconButton
          ref={downRef}
          label="Move down"
          tone="ghost"
          data-testid="builder-recipe-down"
          disabled={atBottom}
          onClick={() => {
            pendingDir.current = "down";
            onMove(index, index + 1);
          }}
        >
          &#9660;
        </IconButton>
        {/* BU-EMPTYREMOVE: a series keeps at least one member — the last row's
            Remove is disabled (with a reason) so a draft can't be emptied to zero
            and then "saved" as an empty series. */}
        <IconButton
          label={count === 1 ? "Remove sample (a series keeps at least one member)" : "Remove sample"}
          tone="ghost"
          dismiss
          data-testid="builder-recipe-remove"
          disabled={count === 1}
          onClick={onRemove}
        />
      </div>
    </div>
  );
}

function sampleNameMap(corpus: api.CorpusSample[]): Record<number, string> {
  const out: Record<number, string> = {};
  for (const s of corpus) out[s.id] = s.name || `Sample ${s.id}`;
  return out;
}
