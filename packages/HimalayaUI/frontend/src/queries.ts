import { useMemo, useRef } from "react";
import { useQuery, useQueries, useMutation, useQueryClient } from "@tanstack/react-query";
import { authOpts } from "./lib/authOpts";
import * as api from "./api";
import { useAppState } from "./state";
import { getClientId } from "./lib/clientId";
import { newClientOpId } from "./lib/clientOpId";
import { useQueueMutation } from "./lib/queue/useQueueMutation";
import { isSampleScreened } from "./lib/sample/screened";
import {
  updateSampleMutator,
  addSampleTagMutator,
  removeSampleTagMutator,
  addCorpusSampleTagMutator,
  removeCorpusSampleTagMutator,
  addExposureTagMutator,
  removeExposureTagMutator,
  postSampleMessageMutator,
  postComparisonMessageMutator,
  setExposureStatusMutator,
  selectExposureMutator,
} from "./lib/queue/mutators/trivial";
import { peakAddMutator } from "./lib/queue/mutators/peakAdd";
import { peakRemoveMutator } from "./lib/queue/mutators/peakRemove";
import {
  peakExcludeMutator, peakUnexcludeMutator,
} from "./lib/queue/mutators/peakSetExcluded";
import {
  addIndexToGroupMutator,
  removeIndexFromGroupMutator,
  deleteIndexMutator,
} from "./lib/queue/mutators/indexGroup";
import { createSpeculativeMutator } from "./lib/queue/mutators/createSpeculative";
import { reanalyzeExposureMutator } from "./lib/queue/mutators/reanalyzeExposure";
import { saveComparisonMutator } from "./lib/queue/mutators/saveComparison";
import { deleteComparisonMutator } from "./lib/queue/mutators/deleteComparison";
import { scopeSeriesMutator } from "./lib/queue/mutators/scopeSeries";
import { saveSeriesMutator } from "./lib/queue/mutators/saveSeries";
import { commitSeriesPlateMutator } from "./lib/queue/mutators/commitSeriesPlate";
import { deleteSeriesMutator } from "./lib/queue/mutators/deleteSeries";
import { useExposureHasPendingPeakOps } from "./lib/queue/hooks";

const CLIENT_ID = getClientId();

// Nullable id slots emit `"none"` instead of the id so prefix invalidations
// (e.g. ["exposure", id]) never accidentally match a disabled query's key.
// Centralising means the "build the key inline" gotcha is mechanically gone:
// every hook + mutator goes through these helpers.
export const queryKeys = {
  experiments: ["experiments"] as const,
  experiment: (id: number) => ["experiment", id] as const,
  samples:    (experimentId: number | undefined) =>
    ["experiment", experimentId ?? "none", "samples"] as const,
  exposures:  (sampleId: number | undefined) =>
    ["sample", sampleId ?? "none", "exposures"] as const,
  trace:      (exposureId: number | undefined) =>
    ["exposure", exposureId ?? "none", "trace"] as const,
  peaks:      (exposureId: number | undefined) =>
    ["exposure", exposureId ?? "none", "peaks"] as const,
  indices:    (exposureId: number | undefined) =>
    ["exposure", exposureId ?? "none", "indices"] as const,
  groups:     (exposureId: number | undefined) =>
    ["exposure", exposureId ?? "none", "groups"] as const,
  messages:   (sampleId: number | undefined) =>
    ["sample", sampleId ?? "none", "messages"] as const,
  speculativeSnap: (
    exposureId: number | undefined,
    phase: string | undefined,
    anchorPeakId: number | undefined,
    anchorRatio: number,
  ) => ["exposure", exposureId ?? "none", "speculative-snap",
        phase ?? "", anchorPeakId ?? -1, anchorRatio] as const,
  // Single-entity keys are namespaced with `-entity` to avoid prefix-matching
  // collisions with the existing collection keys (e.g., a future
  // invalidate(["exposure", id]) would otherwise also blast peaks/indices/groups).
  peak:     (id: number | undefined) => ["peak-entity", id ?? "none"] as const,
  index:    (id: number | undefined) => ["index-entity", id ?? "none"] as const,
  exposure: (id: number | undefined) => ["exposure-entity", id ?? "none"] as const,
  sample:   (id: number | undefined) => ["sample-entity", id ?? "none"] as const,
  // Compare page (Plan §Phase 3). Listing key is parameterized by scope —
  // pass "all" for the global listing, an experimentId for the per-experiment
  // listing. Membership-derived listings can change in either direction when
  // ANY exposure-touching event lands, so the SSE handler invalidates both
  // forms with a prefix `["comparisons"]` invalidation.
  comparisons:        (scope: number | "all") => ["comparisons", scope] as const,
  comparison:         (id: number | undefined) =>
    ["comparison", id ?? "none"] as const,
  comparisonMembers:  (id: number | undefined) =>
    ["comparison", id ?? "none", "members"] as const,
  comparisonForks:    (id: number | undefined) =>
    ["comparison", id ?? "none", "forks"] as const,
  comparisonMessages: (id: number | undefined) =>
    ["comparison", id ?? "none", "messages"] as const,
  // Series (#166/#167/#168). Detail root `["series", id]`, listing root
  // `["series-list"]` — distinct roots so a listing invalidation never
  // clobbers a detail entry (mirrors the comparison/comparisons split). Read
  // hooks (useSeriesList / useSeries) — see below (I3.3).
  series:     (id: number | undefined) => ["series", id ?? "none"] as const,
  seriesList: ["series-list"] as const,
  seriesPins: ["series-pins"] as const,
  // Corpus scoping reads (I3.4 #174): the ordering-variable proposal source
  // and the corpus picker projection. Foreign add_tag replay invalidates both
  // (applyRemoteToCache add_tag/sample branch) so a peer's scoping write
  // refreshes the /series/new proposal.
  corpusSampleTags:    ["corpus-sample-tags"] as const,
  corpusPickerSamples: ["corpus-picker-samples"] as const,
  // Picker support routes (Plan §Phase 5, Task 5.2). Both are read-only —
  // `recentlyPickedExposures` is per-user across all experiments; `sampleTags`
  // is per-experiment (distinct (key, value) pairs).
  recentlyPickedExposures: (userId: number | undefined, limit: number) =>
    ["user", userId ?? "none", "recently-picked-exposures", limit] as const,
  sampleTags: (experimentId: number | undefined) =>
    ["experiment", experimentId ?? "none", "sample-tags"] as const,
  pickerSamples: (experimentId: number | undefined) =>
    ["experiment", experimentId ?? "none", "picker-samples"] as const,
  // Phase 13 — comparison pins, scoped per-user via the X-Username header
  // (no userId in the key — the cache row is implicitly per-tab/per-username).
  comparisonPins: ["comparison-pins"] as const,
  // Corpus-wide sample listing (redesign Phase 1, #156). Distinct
  // ["corpus", ...] namespace so corpus add_tag (#159) patches and
  // invalidations never collide with the per-experiment
  // ["experiment", id, "samples"] entries that `samples` produces.
  corpusSamples: ["corpus", "samples"] as const,
};

export function useExperiments() {
  return useQuery({
    queryKey: queryKeys.experiments,
    queryFn: () => api.listExperiments(),
  });
}

export function useExperiment(id: number) {
  return useQuery({
    queryKey: queryKeys.experiment(id),
    queryFn: () => api.getExperiment(id),
    // Callers pass `id ?? 0` when no experiment is active; gate the fetch so
    // we don't hit GET /api/experiments/0 → 404 (and retries) on every chip
    // mount before the user picks an experiment.
    enabled: id > 0,
  });
}

export function useSamples(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.samples(experimentId),
    queryFn: () => api.listSamples(experimentId),
  });
}

/**
 * Corpus-wide sample list — every sample across all experiments, each with
 * its tags and a q_units from the owning experiment's config. Backs the
 * Phase 1 contact sheet (#160) and sample loupe (#161). Distinct from the
 * per-experiment useSamples(experimentId): no parameter, so no `enabled`
 * gate — the fetch runs unconditionally.
 */
export function useCorpusSamples() {
  return useQuery({
    queryKey: queryKeys.corpusSamples,
    queryFn: () => api.listCorpusSamples(),
  });
}

export function useExposures(sampleId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.exposures(sampleId),
    queryFn: () => api.listExposures(sampleId as number),
    enabled: sampleId !== undefined,
  });
}

/**
 * Page-level screened-progress aggregate for the contact-sheet header (M-1):
 * "N / M samples screened". Runs `useQueries` over the visible samples,
 * sharing the exact `queryKeys.exposures(id)` cache rows each ContactSheetRow
 * fetches — no double-fetch. `screened` is derived per-sample by the shared
 * helper (lib/sample/screened.ts), so it tracks #162's flag when that lands.
 *
 * `total` is the number of samples passed in; `screened` counts those whose
 * derivation resolves true once their exposures have loaded.
 */
export function useScreenedProgress(
  samples: readonly api.CorpusSample[],
): { screened: number; total: number } {
  const queries = useQueries({
    queries: samples.map((s) => ({
      queryKey: queryKeys.exposures(s.id),
      queryFn: () => api.listExposures(s.id),
    })),
  });
  const screened = samples.reduce(
    (n, s, i) => (isSampleScreened(s, queries[i]?.data) ? n + 1 : n),
    0,
  );
  return { screened, total: samples.length };
}

export function useTrace(exposureId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.trace(exposureId),
    queryFn: () => api.getTrace(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

/**
 * Run a parallel `useQueries` over a list of ids and return a stable
 * `Map<id, T>`. Stability is load-bearing: TanStack reuses `query.data` refs
 * when underlying rows don't change, so we hash on `(id, dataRef)` and only
 * mint a new Map when some pair actually changed. Without this, every parent
 * re-render produces a new Map and downstream `useCallback`s that depend on
 * it tear down + replot — wheel/brush smoothness in MultiTracePlot depends
 * on this invariant.
 *
 * Used by the three `useMember*` hooks below (traces, exposures, samples)
 * so any future tweak to the identity-comparison logic is applied once.
 */
function useStableQueryMap<T>(
  ids: number[],
  buildOptions: (id: number) => { queryKey: readonly unknown[]; queryFn: () => Promise<T> },
): Map<number, T> {
  const queries = useQueries({
    queries: ids.map((id) => buildOptions(id)),
  });
  // Collapse (ids, data refs) into a single fixed-length string so `useMemo`'s
  // deps comparison stays well-defined — variable-length deps arrays break
  // memoisation when the previous length was 0 (React's elementwise loop
  // terminates early and returns the stale cached value). A per-hook-instance
  // WeakMap assigns each fresh data object a stable nonce so the signature
  // tracks REF identity (TanStack reuses `q.data` refs when nothing changed),
  // which is what wheel-/brush-smoothness in MultiTracePlot depends on. The
  // WeakMap write is idempotent (`set(k, n)` with the same `n`), so this is
  // StrictMode- and Concurrent-render-safe.
  const refTable = useRef<{ map: WeakMap<object, number>; next: number } | null>(null);
  if (refTable.current === null) refTable.current = { map: new WeakMap(), next: 0 };
  const t = refTable.current;
  const nonceOf = (v: unknown): string => {
    if (v === undefined) return "_";
    if (v === null || typeof v !== "object") return `p${String(v)}`;
    const obj = v as object;
    let n = t.map.get(obj);
    if (n === undefined) {
      n = t.next++;
      t.map.set(obj, n);
    }
    return `o${n}`;
  };
  const signature = ids.map((id, i) => `${id}:${nonceOf(queries[i]?.data)}`).join("|");
  // eslint-disable-next-line react-hooks/exhaustive-deps
  return useMemo(() => {
    const m = new Map<number, T>();
    for (let i = 0; i < ids.length; i++) {
      const d = queries[i]?.data;
      if (d !== undefined) m.set(ids[i]!, d as T);
    }
    return m;
  }, [signature]);
}

/**
 * Fetch live `(q, I)` traces for a variable list of exposure ids in parallel.
 * Used by the Compare page MultiTracePlot — one trace per member. Returns a
 * `Map<exposure_id, Trace>` for compose-time consumption by `MultiTracePlot`.
 *
 * Cache keys mirror `useTrace` exactly so single-exposure pages and Compare
 * share the same cache row (no double-fetching across pages).
 */
export function useMemberTraces(exposureIds: number[]): Map<number, api.Trace> {
  return useStableQueryMap(exposureIds, (id) => ({
    queryKey: queryKeys.trace(id),
    queryFn: () => api.getTrace(id),
  }));
}

/**
 * Sibling of `useMemberTraces` that surfaces a single boolean — true when
 * any underlying trace fetch is in its cold-loading state. Used by the
 * Compare-page skeleton wrappers to gate plot + gutter.
 *
 * Mirrors the boneyard rule (CLAUDE.md): gate on `query.isLoading`, NOT
 * `isPending` — disabled queries (empty `exposureIds`) and background
 * refetches must NOT trigger the skeleton.
 */
export function useMemberTracesLoading(exposureIds: number[]): boolean {
  const queries = useQueries({
    queries: exposureIds.map((id) => ({
      queryKey: queryKeys.trace(id),
      queryFn: () => api.getTrace(id),
    })),
  });
  return queries.some((q) => q.isLoading);
}

/**
 * Sibling of `useMemberTraces` that hydrates per-member EXPOSURE rows. Compare
 * review-mode pages use `sampleIdFor` and label resolution that depend on the
 * exposure cache being populated; without an explicit subscription the cache
 * never warms (only the trace key is fetched), and downstream readers fall
 * back to the orphan path. Issue #61 / #52.
 */
export function useMemberExposures(exposureIds: number[]): Map<number, api.Exposure> {
  return useStableQueryMap(exposureIds, (id) => ({
    queryKey: queryKeys.exposure(id),
    queryFn: () => api.getExposure(id),
  }));
}

/**
 * Hydrates per-member SAMPLE rows. Used together with `useMemberExposures`
 * by the Compare review-mode label resolver — a member's display label is
 * `${sample.display_name || sample.name} · ${exposure.filename}` (issue #52).
 */
export function useMemberSamples(sampleIds: number[]): Map<number, api.Sample> {
  return useStableQueryMap(sampleIds, (id) => ({
    queryKey: queryKeys.sample(id),
    queryFn: () => api.getSample(id),
  }));
}

export function usePeaks(exposureId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.peaks(exposureId),
    queryFn: () => api.listPeaks(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

export function useIndices(exposureId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.indices(exposureId),
    queryFn: () => api.listIndices(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

export function useAddPeak(exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    peakAddMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (q: number) => inner.mutate({ q }),
  };
}

export function useRemovePeak(exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    peakRemoveMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (peakId: number) => inner.mutate({ peakId }),
  };
}

// `useSetPeakExcluded` decomposes into two mutators because the backend uses
// distinct event kinds (`peak_excluded`, `peak_unexcluded`). The hook routes
// the call to the correct mutator based on `excluded` and presents a unified
// `mutate({ peakId, excluded })` surface so existing consumers stay unchanged.
export function useSetPeakExcluded(exposureId: number) {
  const username = useAppState((s) => s.username);
  const exclude = useQueueMutation(
    peakExcludeMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
  const unexclude = useQueueMutation(
    peakUnexcludeMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
  return {
    mutate: (input: { peakId: number; excluded: boolean }) => {
      if (input.excluded) exclude.mutate({ peakId: input.peakId });
      else unexclude.mutate({ peakId: input.peakId });
    },
    isPending: exclude.isPending || unexclude.isPending,
    error: exclude.error ?? unexclude.error,
    reset: () => { exclude.reset(); unexclude.reset(); },
  };
}

export function useReanalyzeExposure(exposureId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(
    reanalyzeExposureMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
}

export function useGroups(exposureId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.groups(exposureId),
    queryFn: () => api.listGroups(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

export function useAddIndexToGroup(exposureId: number, groupId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    addIndexToGroupMutator,
    { exposureId, groupId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (indexId: number) => inner.mutate({ indexId }),
  };
}

export function useRemoveIndexFromGroup(exposureId: number, groupId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    removeIndexFromGroupMutator,
    { exposureId, groupId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (indexId: number) => inner.mutate({ indexId }),
  };
}

// Speculative-snap is a query keyed on (exposureId, phase, anchorPeakId, anchorRatio).
// The hook is enabled-gated because the builder calls it after a phase + anchor
// are both chosen — never on partial input. M2.4 also gates on
// `useExposureHasPendingPeakOps`: while a peak-affecting op is in flight the
// snap response would reflect a stale peak set, so we suspend the query and
// the SpeculativeBuilder keeps showing the last good snap with an "updating"
// indicator until the op settles.
export function useSpeculativeSnap(
  exposureId: number | undefined,
  phase: string | undefined,
  anchorPeakId: number | undefined,
  anchorRatio: number,
) {
  const blocked = useExposureHasPendingPeakOps(exposureId);
  return useQuery({
    queryKey: queryKeys.speculativeSnap(exposureId, phase, anchorPeakId, anchorRatio),
    queryFn: () => api.getSpeculativeSnap(exposureId as number, phase as string, anchorPeakId as number, anchorRatio),
    enabled: exposureId !== undefined && phase !== undefined && anchorPeakId !== undefined && !blocked,
  });
}

export function useCreateSpeculative(exposureId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(
    createSpeculativeMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
}

export function useDeleteIndex(exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    deleteIndexMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (indexId: number) => inner.mutate({ indexId }),
  };
}

export function useUpdateSample(experimentId: number, sampleId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(
    updateSampleMutator,
    { experimentId, sampleId, username, clientId: CLIENT_ID },
  );
}

export function useAddSampleTag(experimentId: number, sampleId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(
    addSampleTagMutator,
    { experimentId, sampleId, username, clientId: CLIENT_ID },
  );
}

export function useSampleMessages(sampleId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.messages(sampleId),
    queryFn: () => api.listSampleMessages(sampleId as number),
    enabled: sampleId !== undefined,
  });
}

// Adapter: useQueueMutation flat-spreads input into the payload, but the
// existing consumers call `postMsg.mutate(body: string)`. Wrap to package the
// string under a `body` key for the mutator.
export function usePostSampleMessage(sampleId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    postSampleMessageMutator,
    { sampleId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (body: string) => inner.mutate({ body }),
  };
}

export function useRemoveSampleTag(experimentId: number, sampleId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    removeSampleTagMutator,
    { experimentId, sampleId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (tagId: number) => inner.mutate({ tagId }),
  };
}

/**
 * Corpus contact-sheet sample tagging (#159). Scope carries `sampleId` only —
 * NO `experimentId` — which is how resolveMutator routes the op to the
 * corpus mutator rather than the experiment-scoped or exposure mutator.
 */
export function useAddCorpusSampleTag(sampleId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(
    addCorpusSampleTagMutator,
    { sampleId, username, clientId: CLIENT_ID },
  );
}

export function useRemoveCorpusSampleTag(sampleId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    removeCorpusSampleTagMutator,
    { sampleId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (tagId: number) => inner.mutate({ tagId }),
  };
}

export function useSetExposureStatus(sampleId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(setExposureStatusMutator, { sampleId, username, clientId: CLIENT_ID });
}

export function useSelectExposure(sampleId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(selectExposureMutator, { sampleId, username, clientId: CLIENT_ID });
  return {
    ...inner,
    mutate: (exposureId: number) => inner.mutate({ exposureId }),
  };
}

export function useAddExposureTag(sampleId: number, exposureId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(
    addExposureTagMutator,
    { sampleId, exposureId, username, clientId: CLIENT_ID },
  );
}

export function useRemoveExposureTag(sampleId: number, exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    removeExposureTagMutator,
    { sampleId, exposureId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (tagId: number) => inner.mutate({ tagId }),
  };
}

export function usePeak(id: number | undefined) {
  return useQuery({
    queryKey: queryKeys.peak(id),
    queryFn: () => api.getPeak(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useIndex(id: number | undefined) {
  return useQuery({
    queryKey: queryKeys.index(id),
    queryFn: () => api.getIndex(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useExposure(id: number | undefined) {
  return useQuery({
    queryKey: queryKeys.exposure(id),
    queryFn: () => api.getExposure(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useSampleById(id: number | undefined) {
  return useQuery({
    queryKey: queryKeys.sample(id),
    queryFn: () => api.getSample(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

// ─── Comparisons (Plan §Phase 4, Task 4.1 Step 2/3) ────────────────────────

/**
 * List comparisons for the given scope.
 * - scope = experiment id → `/api/experiments/:eid/comparisons`
 * - scope = "all"         → `/api/comparisons`
 *
 * Returns the same `ComparisonSummary[]` shape from both routes so the
 * sidebar doesn't have to branch on scope.
 */
export function useComparisons(scope: number | "all") {
  return useQuery({
    queryKey: queryKeys.comparisons(scope),
    queryFn: () => scope === "all"
      ? api.listComparisons()
      : api.listExperimentComparisons(scope),
  });
}

export function useComparison(id: number | undefined) {
  return useQuery({
    queryKey: queryKeys.comparison(id),
    queryFn: () => api.getComparison(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useComparisonForks(id: number | undefined) {
  return useQuery({
    queryKey: queryKeys.comparisonForks(id),
    queryFn: () => api.getComparisonForks(id as number),
    enabled: id !== undefined,
  });
}

/**
 * Corpus-wide series listing — every saved series across all experiments,
 * sorted by recency (backend `last_event_at DESC`). Backs the series folio
 * (#173 / I3.3). No parameter, so no `enabled` gate — fetches unconditionally,
 * exactly like useCorpusSamples().
 */
export function useSeriesList() {
  return useQuery({
    queryKey: queryKeys.seriesList,
    queryFn: () => api.listSeries(),
  });
}

/** Corpus distinct tag pairs — scoping reads these to propose the ordering
 *  variable. May be empty on a cold corpus (accepted by design, #174). */
export function useCorpusSampleTags() {
  return useQuery({
    queryKey: queryKeys.corpusSampleTags,
    queryFn: () => api.getCorpusSampleTags(),
  });
}

/** Corpus picker projection — scoping's candidate member samples. */
export function useCorpusPickerSamples() {
  return useQuery({
    queryKey: queryKeys.corpusPickerSamples,
    queryFn: () => api.getCorpusPickerSamples(),
  });
}

/** Scoping confirm-and-build write: one batch of (key,value) sample_tags with
 *  source='scoping', through the queue (no-op optimistic; SSE-confirmed). */
export function useScopeSeries() {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(scopeSeriesMutator, { username, clientId: CLIENT_ID });
  return {
    ...inner,
    mutate: (input: { key: string; tags: { sampleId: number; value: string }[] }) =>
      inner.mutate(input),
  };
}

/**
 * Series recipe-save (I3.5b). `series_save` → POST /api/series (create) or
 * PATCH /api/series/:id (recipe edit; `payload.id` discriminates). Optimistic
 * surface is the `seriesDraft`; the mutator does no optimistic cache write.
 * Recipe save never 409s (PATCH does not read `expected_content_hash`).
 */
export function useSaveSeries() {
  const username = useAppState((s) => s.username);
  return useQueueMutation(saveSeriesMutator, { username, clientId: CLIENT_ID });
}

/**
 * Series plate-commit (I3.5b). `series_commit` → POST /api/series/:id/commit.
 * Spinner (no optimistic write). On 409 the fetcher throws `ConflictError`,
 * which the App-level `attachConflictBridge` routes to `pendingConflict`
 * (kind `series_commit`) — bridge-free here, single writer, mirroring
 * `useSaveComparison`.
 */
export function useCommitSeriesPlate() {
  const username = useAppState((s) => s.username);
  return useQueueMutation(commitSeriesPlateMutator, { username, clientId: CLIENT_ID });
}

/** Series delete (I3.5b). `series_delete` → DELETE /api/series/:id. */
export function useDeleteSeries() {
  const username = useAppState((s) => s.username);
  return useQueueMutation(deleteSeriesMutator, { username, clientId: CLIENT_ID });
}

/**
 * One series' full nested detail (members + samples). Gated on a defined id
 * so an undefined route param doesn't fire GET /api/series/undefined. Backs
 * the builder (#175) and any folio→detail prefetch.
 */
export function useSeries(id: number | undefined) {
  return useQuery({
    queryKey: queryKeys.series(id),
    queryFn: () => api.getSeries(id as number),
    enabled: id !== undefined,
    // Parity with useComparison (above): a missing/draft series should error
    // fast rather than retry 3× — the I3.5a builder reuses this.
    retry: false,
  });
}

export function useComparisonMessages(id: number | undefined) {
  return useQuery({
    queryKey: queryKeys.comparisonMessages(id),
    queryFn: () => api.listComparisonMessages(id as number),
    enabled: id !== undefined,
  });
}

/**
 * Posts a chat message to the comparison thread. Mirrors
 * `usePostSampleMessage` — the registry discriminates by the presence of
 * `comparisonId` in the payload (vs `sampleId`) to select the right mutator.
 */
export function usePostComparisonMessage(comparisonId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    postComparisonMessageMutator,
    { comparisonId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (body: string) => inner.mutate({ body }),
  };
}

export function useSaveComparison() {
  const username = useAppState((s) => s.username);
  // Phase 12 — bridging ConflictError to Zustand's `pendingConflict` slot
  // happens via a single MutationCache subscriber mounted in App.tsx (see
  // `lib/queue/conflictBridge.ts`). The hook is intentionally bridge-free
  // so multiple mount sites (ComparePageEdit's Save, ConflictModal's
  // Overwrite) cannot race on the slot — there is one subscriber, one
  // writer. Toast suppression on 409 still lives in useQueueMutation's
  // onError (the conflict modal owns that surface).
  return useQueueMutation(
    saveComparisonMutator,
    { username, clientId: CLIENT_ID },
  );
}

export function useDeleteComparison() {
  const username = useAppState((s) => s.username);
  return useQueueMutation(
    deleteComparisonMutator,
    { username, clientId: CLIENT_ID },
  );
}

// ─── Picker support hooks (Plan §Phase 5, Task 5.2) ────────────────────────

/**
 * Fetches the user's most-recently-picked exposures (across all comparisons
 * and experiments). Used by the comparison picker's "Recently used" section.
 * Disabled until `userId` is defined so an empty user state doesn't fire a
 * GET /api/users/undefined/recently-picked-exposures.
 */
export function useRecentlyPickedExposures(
  userId: number | undefined, limit = 20,
) {
  return useQuery({
    queryKey: queryKeys.recentlyPickedExposures(userId, limit),
    queryFn: () => api.getRecentlyPickedExposures(userId as number, limit),
    enabled: userId !== undefined,
  });
}

/**
 * Fetches distinct `(key, value)` sample-tag pairs for an experiment. Used
 * by the picker's tag-filter dropdown. Empty list when no tags exist.
 */
export function useSampleTags(experimentId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.sampleTags(experimentId),
    queryFn: () => api.getSampleTags(experimentId as number),
    enabled: experimentId !== undefined,
  });
}

/**
 * Picker primary list. Returns one row per sample with the resolved
 * indexing-exposure id frozen at fetch time. Spec §PR1 backend.
 *
 * `enabled: experimentId !== undefined` matches `useSampleTags` — picker is
 * always opened from an experiment context, but the hook is shaped to
 * accept `undefined` so render isn't gated on experiment selection.
 */
export function usePickerSamples(experimentId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.pickerSamples(experimentId),
    queryFn: () => api.getPickerSamples(experimentId as number),
    enabled: experimentId !== undefined,
  });
}

// ─── Comparison pins (Plan §Phase 13, Task 13.2) ───────────────────────────
//
// Pin/unpin are trivial idempotent state toggles that don't go through
// `useQueueMutation` (no cross-tab SSE, no `client_op_id` keying — the
// backend uses plain `INSERT OR REPLACE` / `DELETE` with no event-log
// emit). Plain `useMutation` + `invalidateQueries` is the cleanest fit.
//
// The cache key (`queryKeys.comparisonPins`) is global per-tab; if the
// user changes their X-Username mid-session, a manual invalidation is
// needed (out of scope — the username flow already requires a logout
// → re-onboarding round trip).

/** List of comparison ids the current user has pinned (most-recent first). */
export function useComparisonPins() {
  const username = useAppState((s) => s.username);
  return useQuery({
    queryKey: queryKeys.comparisonPins,
    queryFn: () => api.listComparisonPins(authOpts(username, CLIENT_ID)),
    enabled: username !== undefined,
  });
}

export function usePinComparison() {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    // Mint clientOpId per call (not per hook mount) so retries reuse one
    // idempotency key — without it, a network-blip retry would write a
    // duplicate `comparison_pinned` row in user_actions even though the
    // view-table state is already correct.
    mutationFn: (id: number) =>
      api.pinComparison(id, authOpts(username, CLIENT_ID, newClientOpId())),
    onSuccess: () => qc.invalidateQueries({ queryKey: queryKeys.comparisonPins }),
  });
}

export function useUnpinComparison() {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (id: number) =>
      api.unpinComparison(id, authOpts(username, CLIENT_ID, newClientOpId())),
    onSuccess: () => qc.invalidateQueries({ queryKey: queryKeys.comparisonPins }),
  });
}
