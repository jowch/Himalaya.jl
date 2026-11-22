import type { QueryClient } from "@tanstack/react-query";
import * as api from "../../../api";
import type { CompareChange, CompareState } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import { applyCompareChange, EMPTY_COMPARES } from "../../compares";
import type { Mutator, RollbackContext } from "../types";

interface CompareScope {
  username: string | undefined;
  clientId: string;
  experimentId: number;
}

type CompareInput = { sampleId: number; change: CompareChange };

function isCompareState(v: unknown): v is CompareState {
  return typeof v === "object" && v !== null && Array.isArray((v as CompareState).compares) && Array.isArray((v as CompareState).subtractions);
}

export const patchComparesMutator: Mutator<CompareInput, CompareScope, CompareState> = {
  kind: "patch_compares",
  onMutate: (p, qc: QueryClient): RollbackContext => {
    const key = queryKeys.compares(p.experimentId);
    const prev = qc.getQueryData<CompareState>(key);
    qc.setQueryData<CompareState>(key, applyCompareChange(prev ?? EMPTY_COMPARES, p.sampleId, p.change, p.username ?? null));
    return { restore: () => qc.setQueryData(key, prev) };
  },
  request: (p) => api.patchCompares(p.sampleId, p.change, authOpts(p.username, p.clientId, p.clientOpId)),
  onSuccess: (p, response, qc) => {
    if (isCompareState(response)) qc.setQueryData<CompareState>(queryKeys.compares(p.experimentId), response);
    qc.invalidateQueries({ queryKey: queryKeys.compares(p.experimentId) });
  },
};
