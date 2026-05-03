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
import type { Mutator, OpPayload, RollbackContext } from "../types";

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

type AddFlat = OpPayload<AddIndexToGroupInput> & GroupScope & AddIndexToGroupInput;
type RemoveFlat = OpPayload<RemoveIndexFromGroupInput> & GroupScope & RemoveIndexFromGroupInput;
type DeleteFlat = OpPayload<DeleteIndexInput> & ExposureOnlyScope & DeleteIndexInput;

const flatAdd = (p: OpPayload<AddIndexToGroupInput>): AddFlat => p as unknown as AddFlat;
const flatRemove = (p: OpPayload<RemoveIndexFromGroupInput>): RemoveFlat => p as unknown as RemoveFlat;
const flatDelete = (p: OpPayload<DeleteIndexInput>): DeleteFlat => p as unknown as DeleteFlat;

function buildAuth(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

// ---------------------------------------------------------------------------
// addIndexToGroupMutator
// ---------------------------------------------------------------------------

export const addIndexToGroupMutator: Mutator<OpPayload<AddIndexToGroupInput>, GroupEntry> = {
  kind: "index_confirmed",
  onMutate: (raw, qc): RollbackContext => {
    const p = flatAdd(raw);
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
  request: (raw) => {
    const p = flatAdd(raw);
    return api.addIndexToGroup(p.groupId, p.indexId, buildAuth(p));
  },
  onSuccess: (raw, response, qc) => {
    const p = flatAdd(raw);
    const groupsKey = queryKeys.groups(p.exposureId);
    qc.setQueryData<GroupEntry[]>(groupsKey, (old) =>
      (old ?? []).map((g) => (g.id === response.id ? response : g)));
  },
  affectsExposurePeaks: () => false,
};

// ---------------------------------------------------------------------------
// removeIndexFromGroupMutator
// ---------------------------------------------------------------------------

export const removeIndexFromGroupMutator: Mutator<OpPayload<RemoveIndexFromGroupInput>, GroupEntry> = {
  kind: "index_unconfirmed",
  onMutate: (raw, qc): RollbackContext => {
    const p = flatRemove(raw);
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
  request: (raw) => {
    const p = flatRemove(raw);
    return api.removeIndexFromGroup(p.groupId, p.indexId, buildAuth(p));
  },
  onSuccess: (raw, response, qc) => {
    const p = flatRemove(raw);
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

export const deleteIndexMutator: Mutator<OpPayload<DeleteIndexInput>, { deleted: number }> = {
  kind: "delete_index",
  onMutate: (raw, qc): RollbackContext => {
    const p = flatDelete(raw);
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
  request: (raw) => {
    const p = flatDelete(raw);
    return api.deleteIndex(p.indexId, buildAuth(p));
  },
  onSuccess: () => {
    // No-op: optimistic effect already reflects the post-delete state. The
    // backend emits `speculative_deleted` over SSE; applyRemoteToCache
    // invalidates indices + groups for any concurrent reconciliation.
  },
  affectsExposurePeaks: () => false,
};
