/**
 * deleteComparison mutator (Plan §Phase 3, Task 3.2).
 *
 * Mirrors `saveComparison`'s shape but for `DELETE /api/comparisons/:id`.
 * Per spec §"Mutation queue integration": no optimistic write. The
 * comparison entity + members caches are removed (NOT invalidated — refetching
 * a deleted resource 404s and leaves stale error state) and the id is
 * filtered out of every listing key.
 */
import * as api from "../../../api";
import type { AuthOpts, ComparisonSummary } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export interface DeleteComparisonInput {
  id: number;
}

interface DeleteComparisonScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

type DeleteResponse = { id: number; deleted: boolean; event_id: number };

export const deleteComparisonMutator: Mutator<DeleteComparisonInput, DeleteComparisonScope, DeleteResponse> = {
  kind: "comparison_delete",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => api.deleteComparison(p.id, buildAuthOpts(p)),
  onSuccess: (p, _response, qc) => {
    // Remove the entity caches — refetching a deleted resource 404s, which
    // would leave stale error state in cache (`isError: true`) and break
    // any background "is this comparison still around?" check. Use
    // `removeQueries` instead of `invalidateQueries` per spec line 345.
    qc.removeQueries({ queryKey: queryKeys.comparison(p.id) });
    qc.removeQueries({ queryKey: queryKeys.comparisonMembers(p.id) });
    qc.removeQueries({ queryKey: queryKeys.comparisonMessages(p.id) });
    qc.removeQueries({ queryKey: queryKeys.comparisonForks(p.id) });
    // Filter the id out of every cached listing — `setQueriesData` walks the
    // prefix, so `["comparisons", "all"]` and `["comparisons", expId]` both
    // get pruned. This mirrors what the SSE handler does on a foreign delete.
    qc.setQueriesData<ComparisonSummary[]>(
      { queryKey: ["comparisons"] },
      (old) => (old ? old.filter((c) => c.id !== p.id) : old),
    );
  },
  // 404 = "already gone" → desired end state. Without this, a 5xx-then-
  // retry can land a 404 on a comparison the first attempt already deleted,
  // and rollback (which is a no-op here anyway) doesn't matter — but the
  // toast would fire. Suppressing it keeps the UX clean.
  treats404AsSuccess: true,
};
