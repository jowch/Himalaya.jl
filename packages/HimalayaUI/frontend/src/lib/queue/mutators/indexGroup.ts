import * as api from "../../../api";
import type { IndexEntry, AuthOpts, Assignment } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";
type ExposureOnlyScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};
export type DeleteIndexInput = { indexId: number };
function buildAuth(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}
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
  },
  treats404AsSuccess: true,
};