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
import type { IndexEntry, AuthOpts, Assignment } from "../../../api";
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
// TWO cache writes on optimistic apply: filter the index out of
// `queryKeys.indices(exposureId)` AND out of the assignment's members.
// `assignment_members.index_id` is ON DELETE CASCADE (db.jl:223), so deleting
// an index that is in the call removes it from the durable assignment too —
// dropping only from `indices` leaves the cart holding a member id that
// resolves to nothing. Rollback restores both. onSuccess is a no-op: the
// optimistic state already equals the post-delete state.

export const deleteIndexMutator: Mutator<DeleteIndexInput, ExposureOnlyScope, { deleted: number }> = {
  kind: "delete_index",
  onMutate: (p, qc): RollbackContext => {
    const indicesKey = queryKeys.indices(p.exposureId);
    const assignmentKey = queryKeys.assignment(p.exposureId);
    const prevIndices = qc.getQueryData<IndexEntry[]>(indicesKey);
    const prevAssignment = qc.getQueryData<Assignment>(assignmentKey);
    if (prevIndices) {
      qc.setQueryData<IndexEntry[]>(indicesKey,
        prevIndices.filter((ix) => ix.id !== p.indexId));
    }
    // Only patch an assignment the cache already holds — writing one here
    // would seed a synthetic entry and suppress the next real fetch.
    if (prevAssignment) {
      qc.setQueryData<Assignment>(assignmentKey, {
        ...prevAssignment,
        members: prevAssignment.members.filter((m) => m !== p.indexId),
      });
    }
    return {
      restore: () => {
        if (prevIndices !== undefined) qc.setQueryData(indicesKey, prevIndices);
        if (prevAssignment !== undefined) qc.setQueryData(assignmentKey, prevAssignment);
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
