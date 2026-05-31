/**
 * customIndexMutator (Plan D Task D-9, OpKind "custom_index_commit").
 *
 * Commits a client-fitted custom index: POST /api/exposures/:id/custom-index
 * with {phase, basis} (basis = q₁ slope from the modal's physics). The backend
 * persists a speculative index AND adds it to the assignment, emitting two SSE
 * frames (speculative_created + assignment_add).
 *
 * Optimistic policy (review finding #4 / issue-#37 pattern): the new index id
 * does not exist client-side until the server assigns it, so we do NOT splice a
 * negative placeholder into the assignment cache. onMutate is a no-op; onSuccess
 * writes the real IndexEntry into the indices cache and INVALIDATES
 * queryKeys.assignment on the fresh id (the next read picks up the canonical
 * cart with the new member). affectsExposurePeaks stays false.
 */
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

export type CustomIndexInput = { phase: string; basis: number };

function buildAuth(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

export const customIndexMutator: Mutator<CustomIndexInput, CustomIndexScope, api.CustomIndexResponse> = {
  kind: "custom_index_commit",
  onMutate: (): RollbackContext => {
    // No optimistic write: the index id is server-assigned, and splicing a
    // negative placeholder into the assignment would desync on confirm.
    return { restore: () => {} };
  },
  request: (p) => api.createCustomIndex(p.exposureId, p.phase, p.basis, buildAuth(p)),
  onSuccess: (p, response, qc) => {
    const { payload: row } = stripQueueMetadata(response);
    // On the SSE-wins own-tab race the deferred resolves off the FIRST frame
    // (speculative_created), whose mutator has no synthesizeFromSse — so `row`
    // is the generic synth `{event_id, client_op_id, index_id}`, NOT a full
    // IndexEntry. Splicing it would land a phantom `{id:undefined, phase:…}`
    // row (issue-#37). Only splice when the response carries the full shape;
    // otherwise invalidate so the indices cache converges via refetch. (Mirrors
    // createSpeculativeMutator's `hasFullShape` guard.)
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
    // Fresh custom-index id is not yet in the assignment cache — invalidate so
    // the next read reflects the new member (never splice a placeholder).
    qc.invalidateQueries({ queryKey: queryKeys.assignment(p.exposureId) });
  },
  affectsExposurePeaks: () => false,
};
