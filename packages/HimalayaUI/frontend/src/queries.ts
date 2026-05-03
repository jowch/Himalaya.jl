import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import * as api from "./api";
import { useAppState } from "./state";
import { getClientId } from "./lib/clientId";
import { newClientOpId } from "./lib/clientOpId";
import { authOpts } from "./lib/authOpts";
import { useQueueMutation } from "./lib/queue/useQueueMutation";
import {
  updateSampleMutator,
  addSampleTagMutator,
  removeSampleTagMutator,
  addExposureTagMutator,
  removeExposureTagMutator,
  postSampleMessageMutator,
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

function invalidateExposure(qc: ReturnType<typeof useQueryClient>, exposureId: number): void {
  qc.invalidateQueries({ queryKey: queryKeys.peaks(exposureId) });
  qc.invalidateQueries({ queryKey: queryKeys.indices(exposureId) });
  // `groups` must invalidate too: auto-reanalysis re-attaches custom-group
  // members by semantic identity (phase + basis), so the cached `groups`
  // payload — which carries `members` — is stale until refetched. Without
  // this, the right-rail Active set looks empty after every peak edit even
  // though the backend has the correct membership.
  qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) });
}

export function useAddPeak(exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation<{ q: number }, api.PeakAddResponse>(
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
  const inner = useQueueMutation<{ peakId: number }, void>(
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
  const exclude = useQueueMutation<{ peakId: number }, api.PeakUpdatedResponse>(
    peakExcludeMutator,
    { exposureId, username, clientId: CLIENT_ID },
  );
  const unexclude = useQueueMutation<{ peakId: number }, api.PeakUpdatedResponse>(
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
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: () => api.reanalyzeExposure(exposureId, authOpts(username, CLIENT_ID, newClientOpId())),
    onSuccess: () => invalidateExposure(qc, exposureId),
  });
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
  const inner = useQueueMutation<{ indexId: number }, api.GroupEntry>(
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
  const inner = useQueueMutation<{ indexId: number }, api.GroupEntry>(
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
// are both chosen — never on partial input.
export function useSpeculativeSnap(
  exposureId: number | undefined,
  phase: string | undefined,
  anchorPeakId: number | undefined,
  anchorRatio: number,
) {
  return useQuery({
    queryKey: ["exposure", exposureId ?? "none", "speculative-snap", phase ?? "", anchorPeakId ?? -1, anchorRatio] as const,
    queryFn: () => api.getSpeculativeSnap(exposureId as number, phase as string, anchorPeakId as number, anchorRatio),
    enabled: exposureId !== undefined && phase !== undefined && anchorPeakId !== undefined,
  });
}

export function useCreateSpeculative(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (body: {
      phase: string;
      anchor_peak_id: number;
      anchor_ratio: number;
      additional: api.SpeculativeAdditional[];
      active?: boolean;
    }) => api.createSpeculative(exposureId, body, authOpts(username, CLIENT_ID, newClientOpId())),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.indices(exposureId) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) });
    },
  });
}

export function useDeleteIndex(exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation<{ indexId: number }, { deleted: number }>(
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
  return useQueueMutation<{ name?: string; notes?: string }, api.Sample>(
    updateSampleMutator,
    { experimentId, sampleId, username, clientId: CLIENT_ID },
  );
}

export function useAddSampleTag(experimentId: number, sampleId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation<{ key: string; value: string }, api.SampleTag>(
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
  const inner = useQueueMutation<{ body: string }, api.SampleMessage>(
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
  const inner = useQueueMutation<{ tagId: number }, void>(
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
  return useQueueMutation<
    { exposureId: number; status: "accepted" | "rejected" | null },
    { id: number; status: string | null }
  >(setExposureStatusMutator, { sampleId, username, clientId: CLIENT_ID });
}

export function useSelectExposure(sampleId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation<
    { exposureId: number },
    { id: number; selected: boolean }
  >(selectExposureMutator, { sampleId, username, clientId: CLIENT_ID });
  return {
    ...inner,
    mutate: (exposureId: number) => inner.mutate({ exposureId }),
  };
}

export function useAddExposureTag(sampleId: number, exposureId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation<{ key: string; value: string }, api.ExposureTag>(
    addExposureTagMutator,
    { sampleId, exposureId, username, clientId: CLIENT_ID },
  );
}

export function useRemoveExposureTag(sampleId: number, exposureId: number) {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation<{ tagId: number }, void>(
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
