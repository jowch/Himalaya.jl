import { useQuery, useMutation, useQueryClient } from "@tanstack/react-query";
import * as api from "./api";
import { useAppState } from "./state";

export const queryKeys = {
  experiment: (id: number) => ["experiment", id] as const,
  samples:    (experimentId: number) => ["experiment", experimentId, "samples"] as const,
  exposures:  (sampleId: number) => ["sample", sampleId, "exposures"] as const,
  trace:      (exposureId: number) => ["exposure", exposureId, "trace"] as const,
  peaks:      (exposureId: number) => ["exposure", exposureId, "peaks"] as const,
  indices:    (exposureId: number) => ["exposure", exposureId, "indices"] as const,
  groups:     (exposureId: number) => ["exposure", exposureId, "groups"] as const,
};

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

export function useExposures(sampleId: number | undefined) {
  return useQuery({
    queryKey: ["sample", sampleId ?? "none", "exposures"] as const,
    queryFn: () => api.listExposures(sampleId as number),
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
}

function authOpts(username: string | undefined): api.AuthOpts {
  return username !== undefined ? { username } : {};
}

export function useAddPeak(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (q: number) => api.addPeak(exposureId, q, authOpts(username)),
    onSuccess: () => invalidateExposure(qc, exposureId),
  });
}

export function useRemovePeak(exposureId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (peakId: number) => api.removePeak(peakId, authOpts(username)),
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
