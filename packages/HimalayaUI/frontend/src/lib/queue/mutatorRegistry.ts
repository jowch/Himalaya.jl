import type { OpKind, Mutator } from "./types";
import {
  updateSampleMutator,
  addSampleTagMutator,
  removeSampleTagMutator,
  addExposureTagMutator,
  removeExposureTagMutator,
  addCorpusSampleTagMutator,
  removeCorpusSampleTagMutator,
  postSampleMessageMutator,
  postComparisonMessageMutator,
  setExposureStatusMutator,
  selectExposureMutator,
} from "./mutators/trivial";
import { peakAddMutator } from "./mutators/peakAdd";
import { peakRemoveMutator } from "./mutators/peakRemove";
import {
  peakExcludeMutator,
  peakUnexcludeMutator,
} from "./mutators/peakSetExcluded";
import {
  addIndexToGroupMutator,
  removeIndexFromGroupMutator,
  deleteIndexMutator,
} from "./mutators/indexGroup";
import { createSpeculativeMutator } from "./mutators/createSpeculative";
import { reanalyzeExposureMutator } from "./mutators/reanalyzeExposure";
import { saveComparisonMutator } from "./mutators/saveComparison";
import { deleteComparisonMutator } from "./mutators/deleteComparison";
import { saveSeriesMutator } from "./mutators/saveSeries";
import { deleteSeriesMutator } from "./mutators/deleteSeries";

/**
 * Minimal shape required by the resolver: just enough of a persisted op to
 * dispatch. Defined locally so callers (rehydrate, tests) don't have to
 * import the persistence module's internal type.
 */
export interface PersistedOpForResolution {
  kind: OpKind;
  payload: unknown;
}

/**
 * Map a persisted op to the right Mutator. Most kinds have a 1:1 mapping;
 * `add_tag` and `remove_tag` are dual-scope kinds whose dispatch depends on
 * the payload shape — sample-scoped ops carry an `experimentId`, while
 * exposure-scoped ops carry only `sampleId` + `exposureId`. The discriminator
 * here mirrors how `useQueueMutation` flat-spreads `scope` into the payload
 * (see `queries.ts::useAddSampleTag` vs `useAddExposureTag`).
 *
 * Returns `undefined` when the kind has no outbound mutator (currently only
 * `speculative_deleted`, which is server-driven via SSE — the frontend never
 * persists an outbound op of that kind), in which case `rehydrate` will count
 * it as `dropped`.
 */
export function resolveMutator(
  op: PersistedOpForResolution,
): Mutator<any, any, any> | undefined {
  const p = op.payload as
    | {
        experimentId?: number; sampleId?: number;
        exposureId?: number; comparisonId?: number;
      }
    | undefined;
  switch (op.kind) {
    case "update_sample":
      return updateSampleMutator;
    case "add_tag":
      // Tri-scope (#159). Peel exposure off FIRST: an exposure-tag op carries
      // {sampleId, exposureId} and — like a corpus-tag op — has no
      // experimentId, so the two collide unless exposureId is tested first.
      if (p?.exposureId !== undefined) return addExposureTagMutator;
      return p?.experimentId !== undefined
        ? addSampleTagMutator
        : addCorpusSampleTagMutator;
    case "remove_tag":
      if (p?.exposureId !== undefined) return removeExposureTagMutator;
      return p?.experimentId !== undefined
        ? removeSampleTagMutator
        : removeCorpusSampleTagMutator;
    case "post_message":
      return p?.comparisonId !== undefined
        ? postComparisonMessageMutator
        : postSampleMessageMutator;
    case "set_exposure_status":
      return setExposureStatusMutator;
    case "select_exposure":
      return selectExposureMutator;
    case "peak_added":
      return peakAddMutator;
    case "peak_removed":
      return peakRemoveMutator;
    case "peak_excluded":
      return peakExcludeMutator;
    case "peak_unexcluded":
      return peakUnexcludeMutator;
    case "index_confirmed":
      return addIndexToGroupMutator;
    case "index_unconfirmed":
      return removeIndexFromGroupMutator;
    case "delete_index":
      return deleteIndexMutator;
    case "speculative_created":
      return createSpeculativeMutator;
    case "speculative_deleted":
      return undefined; // server-driven; no outbound op of this kind
    case "reanalyze_exposure":
      return reanalyzeExposureMutator;
    case "comparison_save":
      return saveComparisonMutator;
    case "comparison_delete":
      return deleteComparisonMutator;
    case "series_save":
      return saveSeriesMutator;
    case "series_delete":
      return deleteSeriesMutator;
    default:
      return undefined;
  }
}

/**
 * Resolve the mutator that owns a given SSE event kind. Used by the replay
 * coordinator to dispatch `synthesizeFromSse` per event kind.
 *
 * `entity_type` is required because three event kinds (`add_tag`,
 * `remove_tag`, `post_message`) are shared across mutators that differ only
 * by the entity scope they operate on.
 *
 * Returns undefined for unknown event kinds — replayCoordinator will fall
 * back to the generic `{ ...base, ...payload }` shape.
 */
export function resolveMutatorForEvent(
  eventKind: string,
  entityType: string,
): Mutator<any, any, any> | undefined {
  switch (eventKind) {
    case "peak_added":          return peakAddMutator;
    case "peak_removed":        return peakRemoveMutator;
    case "peak_excluded":       return peakExcludeMutator;
    case "peak_unexcluded":     return peakUnexcludeMutator;
    case "index_confirmed":     return addIndexToGroupMutator;
    case "index_unconfirmed":   return removeIndexFromGroupMutator;
    case "speculative_created": return createSpeculativeMutator;
    // event-kind speculative_deleted is the SSE counterpart of op-kind delete_index
    case "speculative_deleted": return deleteIndexMutator;
    case "analyze_run":         return reanalyzeExposureMutator;
    case "comparison_created":
    case "comparison_submitted":
      return saveComparisonMutator;
    case "comparison_deleted":  return deleteComparisonMutator;
    case "series_created":
    case "series_recipe_updated":
      return saveSeriesMutator;
    case "series_deleted":
      return deleteSeriesMutator;
    case "update_sample":       return updateSampleMutator;
    case "set_exposure_status": return setExposureStatusMutator;
    case "select_exposure":     return selectExposureMutator;
    // add_tag / remove_tag stay 2-arm here — deliberately NOT tri-scope like
    // resolveMutator above. Verified (#159): a corpus-originated and an
    // experiment-originated sample tag are byte-identical on the SSE wire
    // (same entity_type="sample", same payload — the route always resolves
    // experiment_id from the sample row). There is no corpus-vs-experiment
    // discriminator to key on, and none is needed: this function only picks a
    // synthesizeFromSse for the own-op SSE-confirmation *response shape*.
    // addSampleTagMutator.synthesizeFromSse yields a SampleTag, which the
    // corpus mutator's own onSuccess consumes; the cache patch is done by the
    // pending mutation's own onSuccess (the corpus mutator). A literal third
    // arm would be unreachable dead code. This is a conscious, plan-aware
    // deviation from the literal "convert both to tri-scope" wording of #159
    // / master plan §11. See
    // docs/superpowers/specs/2026-05-18-corpus-sample-tag-mutations-design.md.
    case "add_tag":
      return entityType === "sample" ? addSampleTagMutator : addExposureTagMutator;
    case "remove_tag":
      return entityType === "sample" ? removeSampleTagMutator : removeExposureTagMutator;
    case "post_message":
      // entity_type is the wire string ("sample_message" / "comparison_message"),
      // matching applyRemoteToCache.ts's discriminator. Default the
      // unknown branch to sample (mirroring applyRemoteToCache's default).
      return entityType === "comparison_message"
        ? postComparisonMessageMutator
        : postSampleMessageMutator;
    // Event kinds with no queue mutator fall through here. This includes
    // invalidate-only foreign events like comparison_pinned / comparison_unpinned
    // (handled entirely by applyRemoteToCache — no optimistic outbound op exists).
    // replayCoordinator treats undefined as "use the generic {...base,...payload} shape".
    default:
      return undefined;
  }
}
