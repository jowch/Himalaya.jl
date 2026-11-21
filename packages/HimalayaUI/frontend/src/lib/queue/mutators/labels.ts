import type { QueryClient } from "@tanstack/react-query";
import * as api from "../../../api";
import type { SampleLabel } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import { applyAddLabel, applyRemoveLabel, applyRestoreLabelSuggestion } from "../../labels";
import type { Mutator, RollbackContext } from "../types";

interface LabelScope {
  username: string | undefined;
  clientId: string;
  experimentId: number;
}

type LabelInput = { sampleId: number; label: string };

function patchLabels(qc: QueryClient, experimentId: number, next: (list: SampleLabel[]) => SampleLabel[]): RollbackContext {
  const key = queryKeys.labels(experimentId);
  const prev = qc.getQueryData<SampleLabel[]>(key);
  qc.setQueryData<SampleLabel[]>(key, next(prev ?? []));
  return { restore: () => qc.setQueryData(key, prev) };
}

function settle(qc: QueryClient, experimentId: number, response: unknown): void {
  if (Array.isArray(response)) qc.setQueryData<SampleLabel[]>(queryKeys.labels(experimentId), response as SampleLabel[]);
  qc.invalidateQueries({ queryKey: queryKeys.labels(experimentId) });
}

export const addLabelMutator: Mutator<LabelInput, LabelScope, SampleLabel[]> = {
  kind: "add_label",
  onMutate: (p, qc) => patchLabels(qc, p.experimentId, (list) => applyAddLabel(list, p.sampleId, p.label)),
  request: (p) => api.addLabel(p.sampleId, p.label, authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (p, response, qc) => settle(qc, p.experimentId, response),
};

export const removeLabelMutator: Mutator<LabelInput, LabelScope, SampleLabel[]> = {
  kind: "remove_label",
  onMutate: (p, qc) => patchLabels(qc, p.experimentId, (list) => applyRemoveLabel(list, p.sampleId, p.label)),
  request: (p) => api.removeLabel(p.sampleId, p.label, authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (p, response, qc) => settle(qc, p.experimentId, response),
};

export const restoreLabelSuggestionMutator: Mutator<LabelInput, LabelScope, SampleLabel[]> = {
  kind: "restore_label_suggestion",
  onMutate: (p, qc) => patchLabels(qc, p.experimentId, (list) => applyRestoreLabelSuggestion(list, p.sampleId, p.label)),
  request: (p) => api.restoreLabelSuggestion(p.sampleId, p.label, authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (p, response, qc) => settle(qc, p.experimentId, response),
};
