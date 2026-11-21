import type { OpKind, Mutator } from "./types";
import {
  updateSampleMutator,
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
import { addLabelMutator, removeLabelMutator } from "./mutators/labels";
import { setLabelColorMutator } from "./mutators/labelColors";

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
 * Map a persisted op to the right Mutator.
 *
 * Returns `undefined` when the kind has no outbound mutator (currently only
 * `speculative_deleted`, which is server-driven via SSE — the frontend never
 * persists an outbound op of that kind), in which case `rehydrate` will count
 * it as `dropped`.
 */
export function resolveMutator(
  op: PersistedOpForResolution,
): Mutator<any, any, any> | undefined {
  switch (op.kind) {
    case "update_sample":
      return updateSampleMutator;
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
    case "bridge_samples":         return bridgeSamplesMutator;
    case "unbridge_sample":        return unbridgeSampleMutator;
    case "add_label":              return addLabelMutator;
    case "remove_label":           return removeLabelMutator;
    case "set_tag_color":          return setLabelColorMutator;
    default:
      return undefined;
  }
}

/**
 * Resolve the mutator that owns a given SSE event kind. Used by the replay
 * coordinator to dispatch `synthesizeFromSse` per event kind.
 *
 * Returns undefined for unknown event kinds — replayCoordinator will fall
 * back to the generic `{ ...base, ...payload }` shape.
 */
export function resolveMutatorForEvent(
  eventKind: string,
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
    case "sample_bridged":           return bridgeSamplesMutator;
    case "sample_unbridged":         return unbridgeSampleMutator;
    case "label_added":              return addLabelMutator;
    case "label_removed":            return removeLabelMutator;
    case "tag_color_set":            return setLabelColorMutator;
    // Event kinds with no queue mutator fall through here (handled entirely by
    // applyRemoteToCache — no optimistic outbound op exists).
    // replayCoordinator treats undefined as "use the generic {...base,...payload} shape".
    default:
      return undefined;
  }
}
