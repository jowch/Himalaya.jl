/**
 * saveComparison mutator (Plan §Phase 3, Task 3.1).
 *
 * One mutator handles BOTH create and update. The discriminator is
 * `payload.id`: absent ⇒ create (POST /api/comparisons); present ⇒ submit
 * (POST /api/comparisons/:id/submit). The OpKind name describes the *user
 * gesture* (Save), not the wire-side event — the backend emits
 * `comparison_created` for the create path and `comparison_submitted` for
 * the submit path.
 *
 * Per spec §"Mutation queue integration": no optimistic write — submission
 * shows a spinner. `onSuccess` writes the canonical comparison + members
 * into the cache and invalidates listing keys (membership-derived listings
 * can change in either direction).
 *
 * On 409, `api.saveComparison` throws a typed `ConflictError` (not a generic
 * ApiError). The error has `status = 409`, so the queue's failure-class
 * router treats it as a validation error (no retry). The conflict modal is
 * driven off the typed throw lifted via `useMutation`'s `error` field — it
 * does NOT route through the toast (validation copy is suppressed for 409
 * via `errors.ts`'s `isValidationError` semantics; see commentary there).
 */
import * as api from "../../../api";
import type {
  AuthOpts, Comparison, ComparisonMemberInput, SaveComparisonBody,
} from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type { Mutator, RollbackContext } from "../types";

export interface SaveComparisonInput {
  /** Absent ⇒ create; present ⇒ submit existing comparison. */
  id?: number;
  title: string;
  description?: string | null;
  members: ComparisonMemberInput[];
  /** Required on submit; absent on create. The frontend captures `baseHash`
   *  at edit-mode entry and rides it as `expected_content_hash`. */
  expected_content_hash?: string;
  forked_from_id?: number | null;
  forked_at_hash?: string | null;
}

interface SaveComparisonScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

function buildBody(p: SaveComparisonInput): SaveComparisonBody {
  const out: SaveComparisonBody = {
    title: p.title,
    members: p.members,
  };
  if (p.description !== undefined) out.description = p.description;
  if (p.expected_content_hash !== undefined) {
    out.expected_content_hash = p.expected_content_hash;
  }
  if (p.forked_from_id !== undefined) out.forked_from_id = p.forked_from_id;
  if (p.forked_at_hash !== undefined) out.forked_at_hash = p.forked_at_hash;
  return out;
}

export const saveComparisonMutator: Mutator<SaveComparisonInput, SaveComparisonScope, Comparison> = {
  kind: "comparison_save",
  // No optimistic write per spec. The user has been seeing the local draft
  // (Zustand) throughout edit mode; the submit shows a spinner instead of
  // an optimistic membership swap. Restore is a no-op.
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) => api.saveComparison(buildBody(p), p.id, buildAuthOpts(p)),
  onSuccess: (_p, response, qc) => {
    // SSE-wins path: `synthesizeResponseFromSse` builds a partial shape
    // (no `forked_from_title`, no per-member `is_stale`/`created_by`, etc.).
    // Detect that case by probing for the full Comparison shape — if the
    // response doesn't carry `members` as an array, fall back to invalidate
    // and let the next read fetch the canonical state.
    const looksFull = Array.isArray((response as { members?: unknown }).members)
      && typeof (response as { content_hash?: unknown }).content_hash === "string";
    if (looksFull) {
      qc.setQueryData(queryKeys.comparison(response.id), response);
      qc.setQueryData(queryKeys.comparisonMembers(response.id), response.members);
    } else if (typeof response?.id === "number") {
      qc.invalidateQueries({ queryKey: queryKeys.comparison(response.id) });
      qc.invalidateQueries({ queryKey: queryKeys.comparisonMembers(response.id) });
    }
    // Invalidate every comparisons listing — comparisons are membership-
    // derived, so the response.members can touch any experiment. Prefix
    // invalidation matches both `["comparisons", "all"]` and any per-
    // experiment scope.
    qc.invalidateQueries({ queryKey: ["comparisons"] });
  },
};
