import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import * as api from "./api";
import { useAppState } from "./state";

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

function authOpts(username: string | undefined): api.AuthOpts {
  return username !== undefined ? { username } : {};
}

// After any peak edit we automatically re-run analysis so the indices reflect
// the user's curation immediately — no "stale" intermediate state. Each peak
// mutation chains: peak op → reanalyze → invalidate exposure queries.
async function autoReanalyze(
  exposureId: number,
  username: string | undefined,
): Promise<void> {
  try {
    await api.reanalyzeExposure(exposureId, authOpts(username));
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
      const peak = await api.addPeak(exposureId, q, authOpts(username));
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
      await api.removePeak(peakId, authOpts(username));
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
      const out = await api.setPeakExcluded(peakId, excluded, authOpts(username));
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
    mutationFn: () => api.reanalyzeExposure(exposureId, authOpts(username)),
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
      api.addIndexToGroup(groupId, indexId, authOpts(username)),
    onSuccess: () => qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) }),
  });
}

export function useRemoveIndexFromGroup(exposureId: number, groupId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (indexId: number) =>
      api.removeIndexFromGroup(groupId, indexId, authOpts(username)),
    onSuccess: () => qc.invalidateQueries({ queryKey: queryKeys.groups(exposureId) }),
  });
}

export function useUpdateSample(experimentId: number, sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (patch: { name?: string; notes?: string }) =>
      api.updateSample(sampleId, patch, authOpts(username)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) }),
  });
}

export function useAddSampleTag(experimentId: number, sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: ({ key, value }: { key: string; value: string }) =>
      api.addSampleTag(sampleId, key, value, authOpts(username)),
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
      api.postSampleMessage(sampleId, body, authOpts(username)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.messages(sampleId) }),
  });
}

export function useRemoveSampleTag(experimentId: number, sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (tagId: number) =>
      api.removeSampleTag(sampleId, tagId, authOpts(username)),
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
    }) => api.setExposureStatus(exposureId, status, authOpts(username)),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: ["sample", sampleId, "exposures", { excludeRejected: false }] });
      qc.invalidateQueries({ queryKey: ["sample", sampleId, "exposures", { excludeRejected: true }] });
    },
  });
}

export function useSelectExposure(sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (exposureId: number) =>
      api.selectExposure(exposureId, authOpts(username)),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: ["sample", sampleId, "exposures", { excludeRejected: false }] });
      qc.invalidateQueries({ queryKey: ["sample", sampleId, "exposures", { excludeRejected: true }] });
    },
  });
}

export function useAddExposureTag(sampleId: number, exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: ({ key, value }: { key: string; value: string }) =>
      api.addExposureTag(exposureId, key, value, authOpts(username)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) }),
  });
}

export function useRemoveExposureTag(sampleId: number, exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (tagId: number) =>
      api.removeExposureTag(exposureId, tagId, authOpts(username)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) }),
  });
}
