/**
 * Index deletion mutator:
 *
 * - deleteIndexMutator       (kind: "delete_index")
 *
 * Deleting an index doesn't change the effective peak set, so it declares
 * `affectsExposurePeaks: () => false` — it doesn't gate the stale-indices
 * banner or the speculative-snap query.
 *
 * The backend emits `speculative_deleted` over SSE, which `applyRemoteToCache`
 * already handles via cache invalidation.
 */
import * as api from "../../../api";
import type { IndexEntry, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

// ---------------------------------------------------------------------------
// Shared scope + auth
// ---------------------------------------------------------------------------

type ExposureOnlyScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};

export type DeleteIndexInput = { indexId: number };

function buildAuth(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

// ---------------------------------------------------------------------------
// deleteIndexMutator
// ---------------------------------------------------------------------------
//
// One cache write on optimistic apply: filter the index out of
// `queryKeys.indices(exposureId)`. Rollback restores it. onSuccess is a no-op —
// optimistic state is already correct; SSE post_state may follow if the delete
// cascades.

export const deleteIndexMutator: Mutator<DeleteIndexInput, ExposureOnlyScope, { deleted: number }> = {
  kind: "delete_index",
  onMutate: (p, qc): RollbackContext => {
    const indicesKey = queryKeys.indices(p.exposureId);
    const prevIndices = qc.getQueryData<IndexEntry[]>(indicesKey);
    if (prevIndices) {
      qc.setQueryData<IndexEntry[]>(indicesKey,
        prevIndices.filter((ix) => ix.id !== p.indexId));
    }
    return {
      restore: () => {
        if (prevIndices !== undefined) qc.setQueryData(indicesKey, prevIndices);
      },
    };
  },
  request: (p) => api.deleteIndex(p.indexId, buildAuth(p)),
  onSuccess: () => {
    // No-op: optimistic effect already reflects the post-delete state. The
    // backend emits `speculative_deleted` over SSE; applyRemoteToCache
    // invalidates indices for any concurrent reconciliation.
  },
  affectsExposurePeaks: () => false,
  // 404 = the index is already gone → desired end state.
  treats404AsSuccess: true,
};
