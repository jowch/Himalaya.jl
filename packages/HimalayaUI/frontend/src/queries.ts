import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import * as api from "./api";
import { useAppState } from "./state";
import { getClientId } from "./lib/clientId";

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

function authOpts(
  username: string | undefined,
  clientId: string | undefined,
): api.AuthOpts {
  const out: api.AuthOpts = {};
  if (username !== undefined) out.username = username;
  if (clientId !== undefined) out.clientId = clientId;
  return out;
}

// After any peak edit we automatically re-run analysis so the indices reflect
// the user's curation immediately — no "stale" intermediate state. Each peak
// mutation chains: peak op → reanalyze → invalidate exposure queries.
async function autoReanalyze(
  exposureId: number,
  username: string | undefined,
): Promise<void> {
  try {
    await api.reanalyzeExposure(exposureId, authOpts(username, CLIENT_ID));
  } catch (e) {
    // Best-effort: surface a console warning but don't block the peak edit.
    // eslint-disable-next-line no-console
    console.warn("auto-reanalyze failed", e);
  }
}

export function useAddPeak(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: async (q: number) => {
      const peak = await api.addPeak(exposureId, q, authOpts(username, CLIENT_ID));
      await autoReanalyze(exposureId, username);
      return peak;
    },
    onSuccess: () => invalidateExposure(qc, exposureId),
  });
}

export function useRemovePeak(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: async (peakId: number) => {
      await api.removePeak(peakId, authOpts(username, CLIENT_ID));
      await autoReanalyze(exposureId, username);
    },
    onSuccess: () => invalidateExposure(qc, exposureId),
  });
}

export function useSetPeakExcluded(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: async ({ peakId, excluded }: { peakId: number; excluded: boolean }) => {
      const out = await api.setPeakExcluded(peakId, excluded, authOpts(username, CLIENT_ID));
      await autoReanalyze(exposureId, username);
      return out;
    },
    onSuccess: () => invalidateExposure(qc, exposureId),
  });
}

export function useReanalyzeExposure(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: () => api.reanalyzeExposure(exposureId, authOpts(username, CLIENT_ID)),
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
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (indexId: number) =>
      api.addIndexToGroup(groupId, indexId, authOpts(username, CLIENT_ID)),
    onSuccess: () => qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) }),
  });
}

export function useRemoveIndexFromGroup(exposureId: number, groupId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (indexId: number) =>
      api.removeIndexFromGroup(groupId, indexId, authOpts(username, CLIENT_ID)),
    onSuccess: () => qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) }),
  });
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
    }) => api.createSpeculative(exposureId, body, authOpts(username, CLIENT_ID)),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.indices(exposureId) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) });
    },
  });
}

export function useDeleteIndex(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (indexId: number) => api.deleteIndex(indexId, authOpts(username, CLIENT_ID)),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.indices(exposureId) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) });
    },
  });
}

export function useUpdateSample(experimentId: number, sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (patch: { name?: string; notes?: string }) =>
      api.updateSample(sampleId, patch, authOpts(username, CLIENT_ID)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) }),
  });
}

export function useAddSampleTag(experimentId: number, sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: ({ key, value }: { key: string; value: string }) =>
      api.addSampleTag(sampleId, key, value, authOpts(username, CLIENT_ID)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) }),
  });
}

export function useSampleMessages(sampleId: number | undefined) {
  return useQuery({
    queryKey: ["sample", sampleId ?? "none", "messages"] as const,
    queryFn: () => api.listSampleMessages(sampleId as number),
    enabled: sampleId !== undefined,
  });
}

export function usePostSampleMessage(sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (body: string) =>
      api.postSampleMessage(sampleId, body, authOpts(username, CLIENT_ID)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.messages(sampleId) }),
  });
}

export function useRemoveSampleTag(experimentId: number, sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (tagId: number) =>
      api.removeSampleTag(sampleId, tagId, authOpts(username, CLIENT_ID)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) }),
  });
}

export function useSetExposureStatus(sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: ({ exposureId, status }: {
      exposureId: number;
      status: "accepted" | "rejected" | null;
    }) => api.setExposureStatus(exposureId, status, authOpts(username, CLIENT_ID)),
    // queryKeys.exposures(id) returns the prefix ["sample", id, "exposures"];
    // TanStack Query invalidateQueries matches by prefix, so this single
    // call covers both excludeRejected variants (and any future variant).
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) }),
  });
}

export function useSelectExposure(sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (exposureId: number) =>
      api.selectExposure(exposureId, authOpts(username, CLIENT_ID)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) }),
  });
}

export function useAddExposureTag(sampleId: number, exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: ({ key, value }: { key: string; value: string }) =>
      api.addExposureTag(exposureId, key, value, authOpts(username, CLIENT_ID)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) }),
  });
}

export function useRemoveExposureTag(sampleId: number, exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (tagId: number) =>
      api.removeExposureTag(exposureId, tagId, authOpts(username, CLIENT_ID)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) }),
  });
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
