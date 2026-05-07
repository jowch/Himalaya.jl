import { useQuery, useQueries } from "@tanstack/react-query";
import * as api from "./api";
import { useAppState } from "./state";
import { getClientId } from "./lib/clientId";
import { useQueueMutation } from "./lib/queue/useQueueMutation";
import {
  updateSampleMutator,
  addSampleTagMutator,
  removeSampleTagMutator,
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
import { useExposureHasPendingPeakOps } from "./lib/queue/hooks";

const CLIENT_ID = getClientId();

export const queryKeys = {
  experiments: ["experiments"] as const,
  experiment: (id: number) => ["experiment", id] as const,
  samples:    (experimentId: number) => ["experiment", experimentId, "samples"] as const,
  exposures:  (sampleId: number) => ["sample", sampleId, "exposures"] as const,
  trace:      (exposureId: number) => ["exposure", exposureId, "trace"] as const,
  peaks:      (exposureId: number) => ["exposure", exposureId, "peaks"] as const,
  indices:    (exposureId: number) => ["exposure", exposureId, "indices"] as const,
  groups:     (exposureId: number) => ["exposure", exposureId, "groups"] as const,
  messages:   (sampleId: number) => ["sample", sampleId, "messages"] as const,
  // Single-entity keys are namespaced with `-entity` to avoid prefix-matching
  // collisions with the existing collection keys (e.g., a future
  // invalidate(["exposure", id]) would otherwise also blast peaks/indices/groups).
  peak:     (id: number) => ["peak-entity", id] as const,
  index:    (id: number) => ["index-entity", id] as const,
  exposure: (id: number) => ["exposure-entity", id] as const,
  sample:   (id: number) => ["sample-entity", id] as const,
  // Compare page (Plan §Phase 3). Listing key is parameterized by scope —
  // pass "all" for the global listing, an experimentId for the per-experiment
  // listing. Membership-derived listings can change in either direction when
  // ANY exposure-touching event lands, so the SSE handler invalidates both
  // forms with a prefix `["comparisons"]` invalidation.
  comparisons:        (scope: number | "all") => ["comparisons", scope] as const,
  comparison:         (id: number) => ["comparison", id] as const,
  comparisonMembers:  (id: number) => ["comparison", id, "members"] as const,
  comparisonForks:    (id: number) => ["comparison", id, "forks"] as const,
  comparisonMessages: (id: number) => ["comparison", id, "messages"] as const,
  // Picker support routes (Plan §Phase 5, Task 5.2). Both are read-only —
  // `recentlyPickedExposures` is per-user across all experiments; `sampleTags`
  // is per-experiment (distinct (key, value) pairs).
  recentlyPickedExposures: (userId: number, limit: number) =>
    ["user", userId, "recently-picked-exposures", limit] as const,
  sampleTags: (experimentId: number) =>
    ["experiment", experimentId, "sample-tags"] as const,
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

export function useExposures(
  sampleId: number | undefined,
  opts?: { excludeRejected?: boolean },
) {
  const excludeRejected = opts?.excludeRejected ?? false;
  return useQuery({
    queryKey: ["sample", sampleId ?? "none", "exposures", { excludeRejected }] as const,
    queryFn: () => api.listExposures(sampleId as number, { excludeRejected }),
    enabled: sampleId !== undefined,
  });
}

export function useTrace(exposureId: number | undefined) {
  return useQuery({
    queryKey: ["exposure", exposureId ?? "none", "trace"] as const,
    queryFn: () => api.getTrace(exposureId as number),
    enabled: exposureId !== undefined,
  });
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
  const queries = useQueries({
    queries: exposureIds.map((id) => ({
      queryKey: ["exposure", id, "trace"] as const,
      queryFn: () => api.getTrace(id),
    })),
  });
  const out = new Map<number, api.Trace>();
  for (let i = 0; i < exposureIds.length; i++) {
    const q = queries[i];
    if (q?.data) out.set(exposureIds[i]!, q.data);
  }
  return out;
}

export function usePeaks(exposureId: number | undefined) {
  return useQuery({
    queryKey: ["exposure", exposureId ?? "none", "peaks"] as const,
    queryFn: () => api.listPeaks(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

export function useIndices(exposureId: number | undefined) {
  return useQuery({
    queryKey: ["exposure", exposureId ?? "none", "indices"] as const,
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
    queryKey: ["exposure", exposureId ?? "none", "groups"] as const,
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
    queryKey: ["exposure", exposureId ?? "none", "speculative-snap", phase ?? "", anchorPeakId ?? -1, anchorRatio] as const,
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
    queryKey: ["sample", sampleId ?? "none", "messages"] as const,
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
    queryKey: id !== undefined ? queryKeys.peak(id) : (["peak-entity", "none"] as const),
    queryFn: () => api.getPeak(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useIndex(id: number | undefined) {
  return useQuery({
    queryKey: id !== undefined ? queryKeys.index(id) : (["index-entity", "none"] as const),
    queryFn: () => api.getIndex(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useExposure(id: number | undefined) {
  return useQuery({
    queryKey: id !== undefined ? queryKeys.exposure(id) : (["exposure-entity", "none"] as const),
    queryFn: () => api.getExposure(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useSampleById(id: number | undefined) {
  return useQuery({
    queryKey: id !== undefined ? queryKeys.sample(id) : (["sample-entity", "none"] as const),
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
    queryKey: id !== undefined ? queryKeys.comparison(id) : (["comparison", "none"] as const),
    queryFn: () => api.getComparison(id as number),
    enabled: id !== undefined,
    retry: false,
  });
}

export function useComparisonForks(id: number | undefined) {
  return useQuery({
    queryKey: id !== undefined ? queryKeys.comparisonForks(id) : (["comparison", "none", "forks"] as const),
    queryFn: () => api.getComparisonForks(id as number),
    enabled: id !== undefined,
  });
}

export function useComparisonMessages(id: number | undefined) {
  return useQuery({
    queryKey: id !== undefined ? queryKeys.comparisonMessages(id) : (["comparison", "none", "messages"] as const),
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
 * and experiments). Used by the ComparisonPicker's "Recently used" section.
 * Disabled until `userId` is defined so an empty user state doesn't fire a
 * GET /api/users/undefined/recently-picked-exposures.
 */
export function useRecentlyPickedExposures(
  userId: number | undefined, limit = 20,
) {
  return useQuery({
    queryKey: userId !== undefined
      ? queryKeys.recentlyPickedExposures(userId, limit)
      : (["user", "none", "recently-picked-exposures", limit] as const),
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
    queryKey: experimentId !== undefined
      ? queryKeys.sampleTags(experimentId)
      : (["experiment", "none", "sample-tags"] as const),
    queryFn: () => api.getSampleTags(experimentId as number),
    enabled: experimentId !== undefined,
  });
}
