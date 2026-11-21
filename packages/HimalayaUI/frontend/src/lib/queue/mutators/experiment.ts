import type { QueryClient } from "@tanstack/react-query";
import * as api from "../../../api";
import type { Experiment } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator } from "../types";

interface ExperimentScope {
  experimentId: number;
  username: string | undefined;
  clientId: string;
}

export type ExperimentPatch = { name?: string; description?: string | null };

export function mergeExperimentPatch(qc: QueryClient, id: number, patch: Partial<Experiment>): void {
  qc.setQueryData<Experiment>(queryKeys.experiment(id), (cur) => (cur ? { ...cur, ...patch } : cur));
  qc.setQueryData<Experiment[]>(queryKeys.experiments, (list) => list?.map((e) => (e.id === id ? { ...e, ...patch } : e)));
}

function patchOf(p: ExperimentPatch): ExperimentPatch {
  const out: ExperimentPatch = {};
  if (p.name !== undefined) out.name = p.name;
  if (p.description !== undefined) out.description = p.description;
  return out;
}

export const patchExperimentMutator: Mutator<ExperimentPatch, ExperimentScope, Experiment> = {
  kind: "patch_experiment",
  onMutate: (p, qc) => {
    const prevOne = qc.getQueryData<Experiment>(queryKeys.experiment(p.experimentId));
    const prevList = qc.getQueryData<Experiment[]>(queryKeys.experiments);
    mergeExperimentPatch(qc, p.experimentId, patchOf(p));
    return {
      restore: () => {
        qc.setQueryData(queryKeys.experiment(p.experimentId), prevOne);
        qc.setQueryData(queryKeys.experiments, prevList);
      },
    };
  },
  request: (p) => api.patchExperiment(p.experimentId, patchOf(p), authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (p, response, qc) => {
    const patch: Partial<Experiment> = {};
    if (response.name !== undefined) patch.name = response.name;
    if (response.description !== undefined) patch.description = response.description;
    mergeExperimentPatch(qc, p.experimentId, patch);
  },
};
