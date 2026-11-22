import * as api from "../../../api";
import type { IndexEntry, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";
export type CreateSpeculativeInput = {
  phase: string;
  anchor_peak_id: number;
  anchor_ratio: number;
  additional: api.SpeculativeAdditional[];
  active?: boolean;
};
type CreateSpeculativeScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};
function buildAuth(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}
export const createSpeculativeMutator: Mutator<
  CreateSpeculativeInput,
  CreateSpeculativeScope,
  IndexEntry
> = {
  kind: "speculative_created",
  onMutate: (p, qc): RollbackContext => {
    const indicesKey = queryKeys.indices(p.exposureId);
    const prevIndices = qc.getQueryData<IndexEntry[]>(indicesKey);
    return {
      restore: () => {
        if (prevIndices !== undefined) qc.setQueryData(indicesKey, prevIndices);
      },
    };
  },
  request: (p) => api.createSpeculative(p.exposureId, {
    phase:          p.phase,
    anchor_peak_id: p.anchor_peak_id,
    anchor_ratio:   p.anchor_ratio,
    additional:     p.additional,
    ...(p.active !== undefined ? { active: p.active } : {}),
  }, buildAuth(p)),
  onSuccess: (p, response, qc) => {
    const indicesKey = queryKeys.indices(p.exposureId);
    const hasFullShape = (response as Partial<IndexEntry>).phase !== undefined;
    if (hasFullShape) {
      qc.setQueryData<IndexEntry[]>(indicesKey, (old) => {
        const list = old ?? [];
        if (list.some((ix) => ix.id === response.id)) return list;
        return [...list, response];
      });
    } else {
      qc.invalidateQueries({ queryKey: indicesKey });
    }
  },
};