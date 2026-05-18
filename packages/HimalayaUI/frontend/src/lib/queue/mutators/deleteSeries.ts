/**
 * deleteSeries mutator (#166). Mirrors `deleteComparison`: no optimistic
 * write; `onSuccess` removes the detail cache (NOT invalidate — refetching a
 * deleted resource 404s) and filters the id out of the listing cache.
 */
import * as api from "../../../api";
import type { AuthOpts } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export interface DeleteSeriesInput {
  id: number;
}

interface DeleteSeriesScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

type DeleteResponse = { id: number; deleted: boolean; event_id: number };

export const deleteSeriesMutator: Mutator<DeleteSeriesInput, DeleteSeriesScope, DeleteResponse> = {
  kind: "series_delete",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => api.deleteSeries(p.id, buildAuthOpts(p)),
  onSuccess: (p, _response, qc) => {
    qc.removeQueries({ queryKey: queryKeys.series(p.id) });
    // Filter the id out of the cached listing. Typed structurally (`{id}`) —
    // the cluster does not own the listing summary type (I3.3 does).
    qc.setQueriesData<{ id: number }[]>(
      { queryKey: queryKeys.seriesList },
      (old) => (old ? old.filter((s) => s.id !== p.id) : old),
    );
  },
  // 404 = "already gone" → desired end state (idempotent delete under retry).
  treats404AsSuccess: true,
  synthesizeFromSse: (remote, base) => ({
    ...base,
    id: remote.entity_id,
    deleted: true,
  } as DeleteResponse),
};
