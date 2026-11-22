import type { QueryClient } from "@tanstack/react-query";
import * as api from "../../../api";
import type { SampleBridge } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import { applyBridge, applyUnbridge } from "../../bridges";
import type { KeepCall } from "../../bridgeCall";
import { invalidateAssignmentsExcept, invalidateSampleCalls } from "../sampleCall";
import type { Mutator, RollbackContext } from "../types";

interface BridgeScope {
  username: string | undefined;
  clientId: string;
  experimentId: number;
}

function patchBridges(qc: QueryClient, experimentId: number, next: (list: SampleBridge[]) => SampleBridge[]): RollbackContext {
  const key = queryKeys.bridges(experimentId);
  const prev = qc.getQueryData<SampleBridge[]>(key);
  qc.setQueryData<SampleBridge[]>(key, next(prev ?? []));
  return { restore: () => qc.setQueryData(key, prev) };
}

function settle(qc: QueryClient, experimentId: number, response: unknown): void {
  if (Array.isArray(response)) qc.setQueryData<SampleBridge[]>(queryKeys.bridges(experimentId), response as SampleBridge[]);
  qc.invalidateQueries({ queryKey: queryKeys.bridges(experimentId) });
  invalidateSampleCalls(qc);
  invalidateAssignmentsExcept(qc, []);
  qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
}

type BridgeInput = { sampleId: number; toSampleId: number; keepCall?: KeepCall };

export const bridgeSamplesMutator: Mutator<BridgeInput, BridgeScope, SampleBridge[]> = {
  kind: "bridge_samples",
  onMutate: (p, qc) => patchBridges(qc, p.experimentId, (list) => applyBridge(list, p.sampleId, p.toSampleId)),
  request: (p) => api.bridgeSamples(p.sampleId, p.toSampleId, p.keepCall, authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (p, response, qc) => settle(qc, p.experimentId, response),
};

type UnbridgeInput = { sampleId: number };

export const unbridgeSampleMutator: Mutator<UnbridgeInput, BridgeScope, SampleBridge[]> = {
  kind: "unbridge_sample",
  onMutate: (p, qc) => patchBridges(qc, p.experimentId, (list) => applyUnbridge(list, p.sampleId)),
  request: (p) => api.unbridgeSample(p.sampleId, authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (p, response, qc) => settle(qc, p.experimentId, response),
};
