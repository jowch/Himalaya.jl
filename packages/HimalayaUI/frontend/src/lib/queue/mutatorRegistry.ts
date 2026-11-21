import type { OpKind, Mutator } from "./types";
import {
  updateSampleMutator,
  addSampleTagMutator,
  removeSampleTagMutator,
  addExposureTagMutator,
  removeExposureTagMutator,
  addCorpusSampleTagMutator,
  removeCorpusSampleTagMutator,
  editCorpusSampleTagMutator,
  setExposureStatusMutator,
  selectExposureMutator,
} from "./mutators/trivial";
import { peakAddMutator } from "./mutators/peakAdd";
import { peakRemoveMutator } from "./mutators/peakRemove";
import {
  peakExcludeMutator,
  peakUnexcludeMutator,
} from "./mutators/peakSetExcluded";
import { deleteIndexMutator } from "./mutators/indexGroup";
import { createSpeculativeMutator } from "./mutators/createSpeculative";
import { reanalyzeExposureMutator } from "./mutators/reanalyzeExposure";
import {
  addAssignmentPhaseMutator,
  removeAssignmentPhaseMutator,
  setSampleCallStateMutator,
} from "./mutators/assignment";
import { customIndexMutator } from "./mutators/customIndex";
import { bridgeSamplesMutator, unbridgeSampleMutator } from "./mutators/bridges";
import { moveExposureMutator, renameSampleMutator, mergeSamplesMutator, splitSampleMutator, dismissGroupingFlagMutator, undoDismissGroupingFlagMutator } from "./mutators/grouping";

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
        exposureId?: number;
      }
    | undefined;
  switch (op.kind) {
    case "update_sample":
      return updateSampleMutator;
    case "add_tag":
      // An exposure-tag op carries {sampleId, exposureId} and —
      // like a corpus-tag op — has no experimentId, so the two collide unless
      // exposureId is tested before the experimentId split.
      if (p?.exposureId !== undefined) return addExposureTagMutator;
      return p?.experimentId !== undefined
        ? addSampleTagMutator
        : addCorpusSampleTagMutator;
    case "remove_tag":
      if (p?.exposureId !== undefined) return removeExposureTagMutator;
      return p?.experimentId !== undefined
        ? removeSampleTagMutator
        : removeCorpusSampleTagMutator;
    case "edit_tag":
      // edit_tag is corpus-only (no experiment-scoped or exposure-scoped path).
      // sampleId-only scope routes here, mirroring the corpus add_tag/remove_tag fallthrough.
      return editCorpusSampleTagMutator;
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
    case "delete_index":
      return deleteIndexMutator;
    case "speculative_created":
      return createSpeculativeMutator;
    case "speculative_deleted":
      return undefined; // server-driven; no outbound op of this kind
    case "reanalyze_exposure":
      return reanalyzeExposureMutator;
    case "assignment_add":
      return addAssignmentPhaseMutator;
    case "assignment_remove":
      return removeAssignmentPhaseMutator;
    case "sample_call_set_state":
      return setSampleCallStateMutator;
    case "custom_index_commit":
      return customIndexMutator;
    case "move_exposure": return moveExposureMutator;
    case "rename_sample": return renameSampleMutator;
    case "merge_samples":          return mergeSamplesMutator;
    case "split_sample":           return splitSampleMutator;
    case "dismiss_grouping_flag":  return dismissGroupingFlagMutator;
    case "undo_dismiss_grouping_flag": return undoDismissGroupingFlagMutator;
    case "bridge_samples":         return bridgeSamplesMutator;
    case "unbridge_sample":        return unbridgeSampleMutator;
    default:
      return undefined;
  }
}

/**
 * Resolve the mutator that owns a given SSE event kind. Used by the replay
 * coordinator to dispatch `synthesizeFromSse` per event kind.
 *
 * `entity_type` is required because two event kinds (`add_tag`,
 * `remove_tag`) are shared across mutators that differ only
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
    case "speculative_created": return createSpeculativeMutator;
    // event-kind speculative_deleted is the SSE counterpart of op-kind delete_index
    case "speculative_deleted": return deleteIndexMutator;
    case "analyze_run":         return reanalyzeExposureMutator;
    case "assignment_add":         return addAssignmentPhaseMutator;
    case "assignment_remove":      return removeAssignmentPhaseMutator;
    case "sample_call_set_state":  return setSampleCallStateMutator;
    // custom_index_commit has no own event kind — its route emits
    // speculative_created + assignment_add, both resolved above. Nothing to add.
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
    // docs/mutation-queue.md.
    case "add_tag":
      return entityType === "sample" ? addSampleTagMutator : addExposureTagMutator;
    case "remove_tag":
      return entityType === "sample" ? removeSampleTagMutator : removeExposureTagMutator;
    case "edit_tag":
      // edit_tag is sample-scoped only. The corpus mutator owns synthesizeFromSse.
      return editCorpusSampleTagMutator;
    // Ingestion grouping-review structural edits (Phase E2). A merge fans out as
    // exposure_moved frames — so exposure_moved is the only SSE kind move AND merge
    // share on the wire. See brief §Step 4 note.
    case "exposure_moved":  return moveExposureMutator;
    case "sample_renamed":  return renameSampleMutator;
    // split emits sample_created (new sample) + sample_split (source update).
    // The own-tab deferred resolves on sample_created (the first-emitted frame).
    // sample_split has no own-op mutator — it routes through applyRemoteToCache.
    case "sample_split":             return splitSampleMutator;
    case "grouping_flag_dismissed":  return dismissGroupingFlagMutator;
    case "sample_bridged":           return bridgeSamplesMutator;
    case "sample_unbridged":         return unbridgeSampleMutator;
    // merge_samples fans out as exposure_moved frames (no sample_merged kind exists);
    // own-op confirmation uses the exposure_moved arm above.
    // Event kinds with no queue mutator fall through here (handled entirely by
    // applyRemoteToCache — no optimistic outbound op exists).
    // replayCoordinator treats undefined as "use the generic {...base,...payload} shape".
    default:
      return undefined;
  }
}
