import * as api from "../../../api";
import type { IndexEntry, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";
import { stripQueueMetadata } from "../queueMeta";
type CustomIndexScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};
export type CustomIndexInput = { phase: string; basis: number; ratios: number[] };
function buildAuth(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}
export const customIndexMutator: Mutator<CustomIndexInput, CustomIndexScope, api.CustomIndexResponse> = {
  kind: "custom_index_commit",
  onMutate: (): RollbackContext => {
    return { restore: () => {} };
  },
  request: (p) => api.createCustomIndex(p.exposureId, p.phase, p.basis, p.ratios, buildAuth(p)),
  onSuccess: (p, response, qc) => {
    const { payload: row } = stripQueueMetadata(response);
    const hasFullShape = (row as Partial<IndexEntry>).phase !== undefined;
    if (hasFullShape) {
      qc.setQueryData<IndexEntry[]>(queryKeys.indices(p.exposureId), (old = []) => {
        const ix = row as IndexEntry;
        if (old.some((e) => e.id === ix.id)) return old;
        return [...old, ix];
      });
    } else {
      qc.invalidateQueries({ queryKey: queryKeys.indices(p.exposureId) });
    }
    qc.invalidateQueries({ queryKey: queryKeys.assignment(p.exposureId) });
  },
};