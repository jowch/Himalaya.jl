import { useEffect } from "react";
import { useQuery, useQueries, useMutation, useQueryClient, type QueryClient } from "@tanstack/react-query";
import * as api from "./api";
import { useAppState } from "./state";
import { getClientId } from "./lib/clientId";
import { useQueueMutation } from "./lib/queue/useQueueMutation";
import {
  setExposureStatusMutator,
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
import { deleteIndexMutator } from "./lib/queue/mutators/indexGroup";
import { bridgeSamplesMutator, unbridgeSampleMutator } from "./lib/queue/mutators/bridges";
import { patchComparesMutator } from "./lib/queue/mutators/compares";
import { addLabelMutator, removeLabelMutator, restoreLabelSuggestionMutator } from "./lib/queue/mutators/labels";
import { setLabelColorMutator } from "./lib/queue/mutators/labelColors";
import { patchExperimentMutator, type ExperimentPatch } from "./lib/queue/mutators/experiment";
import { patchSetupMutator, type SetupPatch } from "./lib/queue/mutators/setup";
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
  problems:   (experimentId: number) => ["experiment", experimentId, "problems"] as const,
  folderPreview: (path: string) => ["folder-preview", path] as const,
  folderPreviews: ["folder-preview"] as const,
  labelColors:  () => ["label-colors"] as const,
  manifest: (dataDir: string, patterns: { image: string; metadata: string; integration: string }, setupFile: string, analysisDir: string) =>
    ["manifest", dataDir, patterns, setupFile, analysisDir] as const,
  bridges:    (experimentId: number | undefined) =>
    ["experiment", experimentId ?? "none", "bridges"] as const,
  compares:   (experimentId: number | undefined) =>
    ["experiment", experimentId ?? "none", "compares"] as const,
  labels:     (experimentId: number | undefined) =>
    ["experiment", experimentId ?? "none", "labels"] as const,
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
  corpusSamples: ["corpus", "samples"] as const,
  experimentCorpusSamples: (experimentId: number) => ["corpus", "samples", experimentId] as const,
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

export function useProblems(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.problems(experimentId),
    queryFn: () => api.listProblems(experimentId),
    enabled: experimentId > 0,
  });
}

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

export function useCompares(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.compares(experimentId),
    queryFn: () => api.listCompares(experimentId),
    enabled: experimentId > 0,
  });
}

export function usePatchCompares(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(patchComparesMutator, { experimentId, username, clientId: CLIENT_ID });
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

export function useLabelColors() {
  return useQuery({ queryKey: queryKeys.labelColors(), queryFn: () => api.listLabelColors() });
}

export function useSetLabelColor() {
  const username = useAppState((s) => s.username);
  return useQueueMutation(setLabelColorMutator, { username, clientId: CLIENT_ID });
}

export function useCorpusSamples({ enabled = true }: { enabled?: boolean } = {}) {
  return useQuery({
    queryKey: queryKeys.corpusSamples,
    queryFn: () => api.listCorpusSamples(),
    enabled,
  });
}

export function useExperimentCorpusSamples(experimentId: number) {
  return useQuery({
    queryKey: queryKeys.experimentCorpusSamples(experimentId),
    queryFn: () => api.listCorpusSamples(experimentId),
    enabled: Number.isInteger(experimentId) && experimentId > 0,
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
    retry: retryUnlessMissing,
  });
}

export function retryUnlessMissing(failures: number, error: unknown): boolean {
  return !(error instanceof api.ApiError && error.status === 404) && failures < 1;
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


export interface EditReads {
  peaks: (exposureId: number) => Promise<api.Peak[]>;
  indices: (exposureId: number) => Promise<api.IndexEntry[]>;
  knownIndices: (exposureId: number) => Promise<api.IndexEntry[]>;
  assignment: (exposureId: number) => Promise<api.Assignment>;
  labels: () => Promise<api.SampleLabel[]>;
  bridges: () => Promise<api.SampleBridge[]>;
  exposures: (sampleId: number) => Promise<api.Exposure[]>;
  labelColors: () => Promise<api.LabelColorEntry[]>;
  experiment: () => Promise<api.Experiment>;
  compares: () => Promise<api.CompareState>;
}

export interface EditOps {
  addPhase: (exposureId: number, indexId: number) => void;
  removePhase: (exposureId: number, indexId: number) => void;
  setCallState: (sampleId: number, state: api.AssignmentState) => void;
  setStatus: (sampleId: number, exposureId: number, status: api.Exposure["status"]) => void;
  addPeak: (exposureId: number, q: number) => void;
  removePeak: (exposureId: number, peakId: number) => void;
  setPeakExcluded: (exposureId: number, peakId: number, excluded: boolean) => void;
  commitCustom: (exposureId: number, phase: string, basis: number, ratios: number[]) => void;
  deleteIndex: (exposureId: number, indexId: number) => void;
  addLabel: (sampleId: number, label: string) => void;
  removeLabel: (sampleId: number, label: string) => void;
  restoreLabelSuggestion: (sampleId: number, label: string) => void;
  setLabelColor: (label: string, color: api.LabelColor | null) => void;
  patchExperiment: (patch: ExperimentPatch) => void;
  patchSetup: (setupId: number, patch: SetupPatch) => void;
  bridge: (sampleId: number, toSampleId: number) => void;
  unbridge: (sampleId: number) => void;
  patchCompares: (sampleId: number, change: api.CompareChange) => void;
  read: EditReads;
  whenSaved: () => Promise<void>;
  saving: () => boolean;
}

function queueBusy(qc: QueryClient): boolean {
  return qc.getMutationCache().getAll().some((m) => m.state.status === "pending" && m.options.scope?.id === "queue");
}

function whenQueueIdle(qc: QueryClient): Promise<void> {
  if (!queueBusy(qc)) return Promise.resolve();
  return new Promise((resolve) => {
    const stop = qc.getMutationCache().subscribe(() => {
      if (queueBusy(qc)) return;
      stop();
      resolve();
    });
  });
}

export function useEditOps(experimentId: number): EditOps {
  const username = useAppState((s) => s.username);
  const qc = useQueryClient();
  const who = { username, clientId: CLIENT_ID };
  const add = useQueueMutation(addAssignmentPhaseMutator, { ...who, exposureId: 0 });
  const remove = useQueueMutation(removeAssignmentPhaseMutator, { ...who, exposureId: 0 });
  const callState = useQueueMutation(setSampleCallStateMutator, { ...who, sampleId: 0 });
  const status = useQueueMutation(setExposureStatusMutator, { ...who, sampleId: 0 });
  const addPeak = useQueueMutation(peakAddMutator, { ...who, exposureId: 0 });
  const removePeak = useQueueMutation(peakRemoveMutator, { ...who, exposureId: 0 });
  const exclude = useQueueMutation(peakExcludeMutator, { ...who, exposureId: 0 });
  const unexclude = useQueueMutation(peakUnexcludeMutator, { ...who, exposureId: 0 });
  const custom = useQueueMutation(customIndexMutator, { ...who, exposureId: 0 });
  const deleteIndex = useQueueMutation(deleteIndexMutator, { ...who, exposureId: 0 });
  const labelScope = { ...who, experimentId };
  const addLabel = useQueueMutation(addLabelMutator, labelScope);
  const removeLabel = useQueueMutation(removeLabelMutator, labelScope);
  const restoreSuggestion = useQueueMutation(restoreLabelSuggestionMutator, labelScope);
  const labelColor = useQueueMutation(setLabelColorMutator, who);
  const experiment = useQueueMutation(patchExperimentMutator, labelScope);
  const setup = useQueueMutation(patchSetupMutator, labelScope);
  const bridge = useQueueMutation(bridgeSamplesMutator, labelScope);
  const unbridge = useQueueMutation(unbridgeSampleMutator, labelScope);
  const compares = useQueueMutation(patchComparesMutator, labelScope);
  const fresh = <T,>(queryKey: readonly unknown[], queryFn: () => Promise<T>): Promise<T> =>
    qc.fetchQuery({ queryKey, queryFn, staleTime: 0 });
  return {
    addPhase: (exposureId, indexId) => add.mutateIn({ ...who, exposureId }, { indexId }),
    removePhase: (exposureId, indexId) => remove.mutateIn({ ...who, exposureId }, { indexId }),
    setCallState: (sampleId, state) => callState.mutateIn({ ...who, sampleId }, { state }),
    setStatus: (sampleId, exposureId, next) => status.mutateIn({ ...who, sampleId }, { exposureId, status: next }),
    addPeak: (exposureId, q) => addPeak.mutateIn({ ...who, exposureId }, { q }),
    removePeak: (exposureId, peakId) => removePeak.mutateIn({ ...who, exposureId }, { peakId }),
    setPeakExcluded: (exposureId, peakId, excluded) =>
      (excluded ? exclude : unexclude).mutateIn({ ...who, exposureId }, { peakId }),
    commitCustom: (exposureId, phase, basis, ratios) => custom.mutateIn({ ...who, exposureId }, { phase, basis, ratios }),
    deleteIndex: (exposureId, indexId) => deleteIndex.mutateIn({ ...who, exposureId }, { indexId }),
    addLabel: (sampleId, label) => addLabel.mutate({ sampleId, label }),
    removeLabel: (sampleId, label) => removeLabel.mutate({ sampleId, label }),
    restoreLabelSuggestion: (sampleId, label) => restoreSuggestion.mutate({ sampleId, label }),
    setLabelColor: (label, color) => labelColor.mutate({ label, color }),
    patchExperiment: (patch) => experiment.mutate(patch),
    patchSetup: (setupId, patch) => setup.mutate({ setupId, ...patch }),
    bridge: (sampleId, toSampleId) => bridge.mutate({ sampleId, toSampleId }),
    unbridge: (sampleId) => unbridge.mutate({ sampleId }),
    patchCompares: (sampleId, change) => compares.mutate({ sampleId, change }),
    read: {
      peaks: (exposureId) => fresh(queryKeys.peaks(exposureId), () => api.listPeaks(exposureId)),
      indices: (exposureId) => fresh(queryKeys.indices(exposureId), () => api.listIndices(exposureId)),
      knownIndices: (exposureId) =>
        qc.ensureQueryData({ queryKey: queryKeys.indices(exposureId), queryFn: () => api.listIndices(exposureId) }),
      assignment: (exposureId) => fresh(queryKeys.assignment(exposureId), () => api.getAssignment(exposureId)),
      labels: () => fresh(queryKeys.labels(experimentId), () => api.listLabels(experimentId)),
      bridges: () => fresh(queryKeys.bridges(experimentId), () => api.listBridges(experimentId)),
      exposures: (sampleId) => fresh(queryKeys.exposures(sampleId), () => api.listExposures(sampleId)),
      labelColors: () => fresh(queryKeys.labelColors(), () => api.listLabelColors()),
      compares: () => fresh(queryKeys.compares(experimentId), () => api.listCompares(experimentId)),
      experiment: () => fresh(queryKeys.experiment(experimentId), () => api.getExperiment(experimentId)),
    },
    whenSaved: () => whenQueueIdle(qc),
    saving: () => queueBusy(qc),
  };
}

export function useSetExposureStatus(sampleId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(setExposureStatusMutator, { sampleId, username, clientId: CLIENT_ID });
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
    meta: { savesWork: false },
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.experiment(experimentId) });
    },
  });
  return {
    ...inner,
    mutate: (force = false) => inner.mutate(force),
  };
}

export function useDeleteExperiment(experimentId: number) {
  const username = useAppState((s) => s.username);
  const qc = useQueryClient();
  return useMutation({
    mutationFn: () => api.deleteExperiment(experimentId, buildAuthOpts(username, CLIENT_ID)),
    onSuccess: () => forgetExperiment(qc, experimentId),
  });
}

export function forgetExperiment(qc: QueryClient, experimentId: number): void {
  qc.removeQueries({ queryKey: queryKeys.experiment(experimentId) });
  void qc.invalidateQueries({ queryKey: queryKeys.experiments });
  void qc.invalidateQueries({ queryKey: queryKeys.folderPreviews });
}

export function useCheckForNewFilesOnOpen(experimentId: number): void {
  const username = useAppState((s) => s.username);
  const qc = useQueryClient();
  const { mutate } = useMutation({
    mutationFn: () => api.triggerScan(experimentId, buildAuthOpts(username, CLIENT_ID), false),
    meta: { savesWork: false },
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.experiment(experimentId) });
    },
  });
  useEffect(() => {
    if (experimentId > 0) mutate();
  }, [experimentId, mutate]);
}

