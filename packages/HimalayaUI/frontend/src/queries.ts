import { useQuery, useQueries, useMutation, useQueryClient } from "@tanstack/react-query";
import * as api from "./api";
import { useAppState } from "./state";
import { getClientId } from "./lib/clientId";
import { useQueueMutation } from "./lib/queue/useQueueMutation";
import {
  addSampleTagMutator,
  removeSampleTagMutator,
  addCorpusSampleTagMutator,
  removeCorpusSampleTagMutator,
  editCorpusSampleTagMutator,
  addExposureTagMutator,
  removeExposureTagMutator,
  setExposureStatusMutator,
  selectExposureMutator,
} from "./lib/queue/mutators/trivial";
import { peakAddMutator } from "./lib/queue/mutators/peakAdd";
import { peakRemoveMutator } from "./lib/queue/mutators/peakRemove";
import {
  peakExcludeMutator, peakUnexcludeMutator,
} from "./lib/queue/mutators/peakSetExcluded";
import {
  addAssignmentPhaseMutator,
  removeAssignmentPhaseMutator,
  setSampleCallStateMutator,
} from "./lib/queue/mutators/assignment";
import { customIndexMutator } from "./lib/queue/mutators/customIndex";
import {
  moveExposureMutator, renameSampleMutator, mergeSamplesMutator,
  splitSampleMutator, dismissGroupingFlagMutator, undoDismissGroupingFlagMutator,
} from "./lib/queue/mutators/grouping";
import { bridgeSamplesMutator, unbridgeSampleMutator } from "./lib/queue/mutators/bridges";
import { addLabelMutator, removeLabelMutator } from "./lib/queue/mutators/labels";
import { setTagColorMutator } from "./lib/queue/mutators/tagColors";
import { authOpts as buildAuthOpts } from "./lib/authOpts";

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
  bridges:    (experimentId: number | undefined) =>
    ["experiment", experimentId ?? "none", "bridges"] as const,
  labels:     (experimentId: number | undefined) =>
    ["experiment", experimentId ?? "none", "labels"] as const,
  tagColors:  () => ["tag-colors"] as const,
  loads:      (id: number | undefined) =>
    ["experiment", id ?? "none", "loads"] as const,
  exposures:  (sampleId: number | undefined) =>
    ["sample", sampleId ?? "none", "exposures"] as const,
  trace:      (exposureId: number | undefined) =>
    ["exposure", exposureId ?? "none", "trace"] as const,
  peaks:      (exposureId: number | undefined) =>
    ["exposure", exposureId ?? "none", "peaks"] as const,
  indices:    (exposureId: number | undefined) =>
    ["exposure", exposureId ?? "none", "indices"] as const,
  assignment: (exposureId: number | undefined) =>
    ["exposure", exposureId ?? "none", "assignment"] as const,
  sampleCall: (sampleId: number | undefined) =>
    ["sample", sampleId ?? "none", "call"] as const,
  // Single-entity keys are namespaced with `-entity` to avoid prefix-matching
  // collisions with the existing collection keys (e.g., a future
  // invalidate(["exposure", id]) would otherwise also blast peaks/indices).
  exposure: (id: number | undefined) => ["exposure-entity", id ?? "none"] as const,
  sample:   (id: number | undefined) => ["sample-entity", id ?? "none"] as const,
  corpusSampleTags:    ["corpus-sample-tags"] as const,
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
    // Self-heal a stuck scan (8c): if the SSE `ingest_complete` frame is dropped
    // (EventSource reconnect mid-scan), poll the resting `ingest_status` so the
    // surface flips to the post-scan review without a manual reload. Polls ONLY
    // while a scan/analyze is in flight — idle/complete/failed rows don't poll.
    // react-query dedupes by query key, so N mounted observers share one poll.
    refetchInterval: (q) => {
      const s = q.state.data?.ingest_status;
      return s === "scanning" || s === "analyzing" ? 2500 : false;
    },
  });
}

/** The Load ▸ Sample ▸ Exposures roll-up (grouping-review + Configuration
 *  acquisition timeline read off this; spec §9.2/§9.6). Gated id>0 like
 *  useExperiment so a zero/undefined experiment never hits /loads.
 *
 *  `refetchWhileScanning` is the robust floor for live-unfold: during an initial
 *  scan, the SSE `ingest_progress`/`ingest_complete` frames invalidate this query
 *  — but React Query ABSORBS an invalidation that lands while a fetch is already
 *  in flight (the frame burst overlaps the first, empty, mid-scan fetch), and
 *  `staleTime` then keeps that empty result "fresh", so nothing re-fetches until a
 *  manual reload. Polling while the scan is live (mirrors useExperiment's ingest
 *  self-heal) picks up the committed loads regardless of the invalidation race;
 *  it stops the moment the caller's `scanning` flag clears. */
export function useBridges(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.bridges(experimentId),
    queryFn: () => api.listBridges(experimentId),
    enabled: experimentId > 0,
  });
}

export function useBridgeSamples(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(bridgeSamplesMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useUnbridgeSample(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(unbridgeSampleMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useLabels(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.labels(experimentId),
    queryFn: () => api.listLabels(experimentId),
    enabled: experimentId > 0,
  });
}

export function useAddLabel(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(addLabelMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useRemoveLabel(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(removeLabelMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useTagColors() {
  return useQuery({ queryKey: queryKeys.tagColors(), queryFn: () => api.listTagColors() });
}

export function useSetTagColor() {
  const username = useAppState((s) => s.username);
  return useQueueMutation(setTagColorMutator, { username, clientId: CLIENT_ID });
}

export function useLoads(id: number, refetchWhileScanning = false) {
  return useQuery({
    queryKey: queryKeys.loads(id),
    queryFn: () => api.listLoads(id),
    enabled: id > 0,
    refetchInterval: refetchWhileScanning ? 2000 : false,
  });
}

export function useSamples(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.samples(experimentId),
    queryFn: () => api.listSamples(experimentId),
    enabled: experimentId > 0,
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

export function useTrace(exposureId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.trace(exposureId),
    queryFn: () => api.getTrace(exposureId as number),
    enabled: exposureId !== undefined,
    retry: (failures, error) => !(error instanceof api.ApiError && error.status === 404) && failures < 1,
  });
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

// Plan D: the assignment cart data layer. The source of the active phase set;
// deriveActiveIndices(assignment, indices) yields the active indices.
export function useAssignment(exposureId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.assignment(exposureId),
    queryFn: () => api.getAssignment(exposureId as number),
    enabled: exposureId !== undefined,
  });
}

export function useSampleCall(sampleId: number | undefined) {
  return useQuery({
    queryKey: queryKeys.sampleCall(sampleId),
    queryFn: () => api.getSampleCall(sampleId as number),
    enabled: sampleId !== undefined,
  });
}

export function useFrameAssignments(exposureIds: number[]) {
  return useQueries({
    queries: exposureIds.map((id) => ({
      queryKey: queryKeys.assignment(id),
      queryFn: () => api.getAssignment(id),
    })),
  });
}

export function useAddAssignmentPhase(exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    addAssignmentPhaseMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
  return { ...inner, mutate: (indexId: number) => inner.mutate({ indexId }) };
}

export function useRemoveAssignmentPhase(exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    removeAssignmentPhaseMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
  return { ...inner, mutate: (indexId: number) => inner.mutate({ indexId }) };
}

export function useRemoveFramePhase() {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    removeAssignmentPhaseMutator,
    { exposureId: 0, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (exposureId: number, indexId: number) =>
      inner.mutate({ exposureId, indexId }),
  };
}

export function useSetSampleCallState(sampleId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    setSampleCallStateMutator,
    { sampleId, username, clientId: CLIENT_ID },
  );
  return { ...inner, mutate: (state: api.AssignmentState) => inner.mutate({ state }) };
}

export function useCommitCustomIndex(exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    customIndexMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (phase: string, basis: number, ratios: number[]) =>
      inner.mutate({ phase, basis, ratios }),
  };
}


export function useAddSampleTag(experimentId: number, sampleId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(
    addSampleTagMutator,
    { experimentId, sampleId, username, clientId: CLIENT_ID },
  );
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

/**
 * Corpus contact-sheet tag editor (#LO-TAGDUP Slice 2). Scope carries `sampleId`
 * only — same routing as useRemoveCorpusSampleTag — which routes the op to
 * editCorpusSampleTagMutator via resolveMutator (edit_tag, sampleId-only).
 */
export function useEditCorpusSampleTag(sampleId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(
    editCorpusSampleTagMutator,
    { sampleId, username, clientId: CLIENT_ID },
  );
  return {
    ...inner,
    mutate: (a: { tagId: number; key: string; value: string }) => inner.mutate(a),
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
    // Threads the optional per-call callbacks through so consumers can toast
    // on confirmation / terminal failure (see MutateCallbacks).
    mutate: (exposureId: number, callbacks?: Parameters<typeof inner.mutate>[1]) =>
      inner.mutate({ exposureId }, callbacks),
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

/** Corpus distinct tag pairs — scoping reads these to propose the ordering
 *  variable. May be empty on a cold corpus (accepted by design, #174). */
export function useCorpusSampleTags() {
  return useQuery({
    queryKey: queryKeys.corpusSampleTags,
    queryFn: () => api.getCorpusSampleTags(),
  });
}

// ---------------------------------------------------------------------------
// Structural-edit hooks (spec §9.6) — one hook per experiment, not per row.
// Entity ids (sampleId / exposureId) ride in the mutate() input, NOT the scope.
// This lets Task 14 instantiate each hook once and dispatch to any row.
// ---------------------------------------------------------------------------

export function useRenameSample(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(renameSampleMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useMoveExposure(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(moveExposureMutator, { experimentId, username, clientId: CLIENT_ID });
}

// merge: scope = experiment + identity ONLY; the caller passes { loserId, survivorId }
// as the mutate() input (MergeSamplesInput, decided in Task 5) — no `as never`.
export function useMergeSamples(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(mergeSamplesMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useSplitSample(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(splitSampleMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useDismissGroupingFlag(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(dismissGroupingFlagMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useUndoDismissGroupingFlag(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(undoDismissGroupingFlagMutator, { experimentId, username, clientId: CLIENT_ID });
}

/** Trigger a rescan of a single experiment (cheap change-check then additive
 *  ingest of new files). Plain TanStack mutation — not a queue mutator because
 *  the scan is idempotent and emits its own SSE progress frames. Invalidates
 *  the experiment detail on success so the header stats / status chip refresh.
 *
 *  `mutate(true)` bypasses the cheap change-check (force flag) — used when a
 *  pattern field edit changes which files the scan discovers. */
export function useTriggerScan(experimentId: number) {
  const username = useAppState((s) => s.username);
  const qc = useQueryClient();
  const inner = useMutation({
    mutationFn: (force: boolean) =>
      api.triggerScan(experimentId, buildAuthOpts(username, CLIENT_ID), force),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.experiment(experimentId) });
    },
  });
  return {
    ...inner,
    mutate: (force = false) => inner.mutate(force),
  };
}

/** Geometry/name override for the Configuration tab (spec §9.6 — E1 ships only
 *  read hooks for experiments). Wraps E1's `updateExperiment(id, patch)` fetcher;
 *  this is a plain TanStack mutation (not a queue mutator — geometry override is
 *  not part of the structural-edit event log) that invalidates the experiment
 *  detail on success so the provenance chips re-render. */
export function useUpdateExperiment(experimentId: number) {
  const username = useAppState((s) => s.username);
  const qc = useQueryClient();
  return useMutation({
    mutationFn: (patch: api.ExperimentPatch) =>
      api.updateExperiment(experimentId, patch, buildAuthOpts(username, CLIENT_ID)),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.experiment(experimentId) });
    },
  });
}

