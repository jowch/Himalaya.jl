import type { Mutator, RollbackContext } from "../types";
import type { AuthOpts } from "../../../api";
import { authOpts } from "../../authOpts";
import { queryKeys } from "../../../queries";
import * as api from "../../../api";

export interface ScopeSeriesInput {
  key: string;
  tags: { sampleId: number; value: string }[];
}
export interface ScopeSeriesScope { username: string | undefined; clientId: string }

// Local auth-opts wrapper — same pattern every mutator file uses (saveSeries,
// trivial, …): wrap the shared authOpts(username, clientId,
// clientOpId) so request() reads the op's metadata off the flat payload.
function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

/**
 * scopeSeriesMutator — series scoping's batch sample_tags write (#174).
 * Queue-routed for `client_op_id` idempotency + multiplayer SSE, but with a
 * NO-OP onMutate: on success the user navigates away (to the folio), so there
 * is no on-surface tag view to optimistically patch/roll back (mirrors
 * saveSeriesMutator). The batch route emits N `add_tag` frames sharing this
 * op's single client_op_id; own-op confirmation clears the deferred after the
 * first frame, the rest self-echo (dropped) — verified by queue-reviewer
 * against idx_events_unique_op (N distinct entity_id rows admitted) + the
 * clear-after-first/self-echo guard. onSuccess invalidates only the two caches
 * this write actually touches.
 *
 * No `synthesizeFromSse`: it is never invoked for this mutator. Foreign
 * `add_tag`(sample) frames route through `resolveMutatorForEvent` to
 * `addSampleTagMutator.synthesizeFromSse` (the registry maps no event kind to
 * scopeSeriesMutator), and `applyRemoteToCache`'s add_tag(sample) branch
 * drives the foreign-tab invalidation. Defining it here would be dead code +
 * a needless `as unknown as` cast.
 */
export const scopeSeriesMutator: Mutator<ScopeSeriesInput, ScopeSeriesScope, api.BatchSampleTagResult[]> = {
  kind: "add_tag",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: (p) =>
    api.batchSampleTags(
      p.key,
      p.tags.map((t) => ({ sample_id: t.sampleId, value: t.value })),
      "scoping",
      buildAuthOpts(p),
    ),
  onSuccess: (_p, _response, qc) => {
    // Narrow: invalidate ONLY the caches this write touches — the corpus tag
    // proposal source and the corpus picker projection (which carries each
    // sample's tags). No broad `["experiment"]` blast (queue anti-pattern);
    // the per-experiment samples list is refreshed by the foreign-replay path
    // in applyRemoteToCache for any tab that has it cached.
    qc.invalidateQueries({ queryKey: queryKeys.corpusSampleTags });
    qc.invalidateQueries({ queryKey: queryKeys.corpusPickerSamples });
  },
};
