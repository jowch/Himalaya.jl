/**
 * speculative_created mutator (M2.4).
 *
 * Speculative-create can affect the indices ordering globally (the new index
 * lands in `indices`, possibly altering score-sort positions, and may opt into
 * the active group). Reconstructing the exact post-create cache shape
 * optimistically would duplicate backend score logic, so we snapshot the prior
 * cache for rollback and let SSE invalidation (in `applyRemoteToCache`) +
 * `onSuccess` produce the final state.
 *
 * The backend POST returns the freshly-built index (same shape as
 * GET /api/indices/:id). On success we splice it into the indices cache and
 * leave groups invalidation to the SSE post-commit broadcast — the response
 * doesn't carry the canonical groups payload.
 */
import * as api from "../../../api";
import type { IndexEntry, GroupEntry, AuthOpts } from "../../../api";
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
    const groupsKey  = queryKeys.groups(p.exposureId);
    const prevIndices = qc.getQueryData<IndexEntry[]>(indicesKey);
    const prevGroups  = qc.getQueryData<GroupEntry[]>(groupsKey);
    return {
      restore: () => {
        if (prevIndices !== undefined) qc.setQueryData(indicesKey, prevIndices);
        if (prevGroups  !== undefined) qc.setQueryData(groupsKey,  prevGroups);
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
    // SSE-wins synthesis only carries `index_id` (the SSE payload for
    // speculative_created omits the IndexEntry shape — see routes_analysis.jl).
    // Detect by checking for `phase`: HTTP response has it, synth doesn't.
    // Without this guard, the splice below pushes a phantom row with
    // id=undefined/phase=undefined/score=undefined into the indices cache.
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
    // Groups membership lives on the backend; invalidate so the active set
    // reflects any auto-add when `active: true`.
    qc.invalidateQueries({ queryKey: queryKeys.groups(p.exposureId) });
  },
  affectsExposurePeaks: () => false,
};
