/**
 * Index/group mutators (M2.3):
 *
 * - addIndexToGroupMutator   (kind: "index_confirmed")
 * - removeIndexFromGroupMutator (kind: "index_unconfirmed")
 * - deleteIndexMutator       (kind: "delete_index")
 *
 * Index/group operations don't change the effective peak set, so all three
 * declare `affectsExposurePeaks: () => false` — they don't gate the stale-
 * indices banner or the speculative-snap query.
 *
 * The delete mutator touches two caches (indices list + group memberships)
 * so its rollback must restore both. The backend emits `speculative_deleted`
 * over SSE, which `applyRemoteToCache` already handles via cache invalidation.
 */
import * as api from "../../../api";
import type { GroupEntry, IndexEntry, AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

// ---------------------------------------------------------------------------
// Shared scope + auth
// ---------------------------------------------------------------------------

type GroupScope = {
  exposureId: number;
  groupId: number;
  username: string | undefined;
  clientId: string;
};
type ExposureOnlyScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};

export type AddIndexToGroupInput = { indexId: number };
export type RemoveIndexFromGroupInput = { indexId: number };
export type DeleteIndexInput = { indexId: number };

function buildAuth(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

// ---------------------------------------------------------------------------
// addIndexToGroupMutator
// ---------------------------------------------------------------------------

export const addIndexToGroupMutator: Mutator<AddIndexToGroupInput, GroupScope, GroupEntry> = {
  kind: "index_confirmed",
  onMutate: (p, qc): RollbackContext => {
    const groupsKey = queryKeys.groups(p.exposureId);
    const prev = qc.getQueryData<GroupEntry[]>(groupsKey);
    if (prev) {
      qc.setQueryData<GroupEntry[]>(groupsKey, prev.map((g) =>
        g.id === p.groupId && !g.members.includes(p.indexId)
          ? { ...g, members: [...g.members, p.indexId] }
          : g,
      ));
    }
    return {
      restore: () => {
        if (prev !== undefined) qc.setQueryData(groupsKey, prev);
      },
    };
  },
  request: (p) => api.addIndexToGroup(p.groupId, p.indexId, buildAuth(p)),
  onSuccess: (p, response, qc) => {
    const groupsKey = queryKeys.groups(p.exposureId);
    qc.setQueryData<GroupEntry[]>(groupsKey, (old) =>
      (old ?? []).map((g) => (g.id === response.id ? response : g)));
  },
  affectsExposurePeaks: () => false,
};

// ---------------------------------------------------------------------------
// removeIndexFromGroupMutator
// ---------------------------------------------------------------------------

export const removeIndexFromGroupMutator: Mutator<RemoveIndexFromGroupInput, GroupScope, GroupEntry> = {
  kind: "index_unconfirmed",
  onMutate: (p, qc): RollbackContext => {
    const groupsKey = queryKeys.groups(p.exposureId);
    const prev = qc.getQueryData<GroupEntry[]>(groupsKey);
    if (prev) {
      qc.setQueryData<GroupEntry[]>(groupsKey, prev.map((g) =>
        g.id === p.groupId
          ? { ...g, members: g.members.filter((m) => m !== p.indexId) }
          : g,
      ));
    }
    return {
      restore: () => {
        if (prev !== undefined) qc.setQueryData(groupsKey, prev);
      },
    };
  },
  request: (p) => api.removeIndexFromGroup(p.groupId, p.indexId, buildAuth(p)),
  onSuccess: (p, response, qc) => {
    const groupsKey = queryKeys.groups(p.exposureId);
    qc.setQueryData<GroupEntry[]>(groupsKey, (old) =>
      (old ?? []).map((g) => (g.id === response.id ? response : g)));
  },
  affectsExposurePeaks: () => false,
};

// ---------------------------------------------------------------------------
// deleteIndexMutator
// ---------------------------------------------------------------------------
//
// Two cache writes on optimistic apply:
//   1. Filter the index out of `queryKeys.indices(exposureId)`.
//   2. Filter `indexId` out of any group's `members` in
//      `queryKeys.groups(exposureId)`.
// Rollback restores both. onSuccess is a no-op — optimistic state is already
// correct; SSE post_state may follow if the delete cascades.

export const deleteIndexMutator: Mutator<DeleteIndexInput, ExposureOnlyScope, { deleted: number }> = {
  kind: "delete_index",
  onMutate: (p, qc): RollbackContext => {
    const indicesKey = queryKeys.indices(p.exposureId);
    const groupsKey = queryKeys.groups(p.exposureId);
    const prevIndices = qc.getQueryData<IndexEntry[]>(indicesKey);
    const prevGroups = qc.getQueryData<GroupEntry[]>(groupsKey);
    if (prevIndices) {
      qc.setQueryData<IndexEntry[]>(indicesKey,
        prevIndices.filter((ix) => ix.id !== p.indexId));
    }
    if (prevGroups) {
      qc.setQueryData<GroupEntry[]>(groupsKey, prevGroups.map((g) => ({
        ...g,
        members: g.members.filter((m) => m !== p.indexId),
      })));
    }
    return {
      restore: () => {
        if (prevIndices !== undefined) qc.setQueryData(indicesKey, prevIndices);
        if (prevGroups !== undefined) qc.setQueryData(groupsKey, prevGroups);
      },
    };
  },
  request: (p) => api.deleteIndex(p.indexId, buildAuth(p)),
  onSuccess: () => {
    // No-op: optimistic effect already reflects the post-delete state. The
    // backend emits `speculative_deleted` over SSE; applyRemoteToCache
    // invalidates indices + groups for any concurrent reconciliation.
  },
  affectsExposurePeaks: () => false,
};
