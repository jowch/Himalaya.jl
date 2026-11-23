import type { OpKind, Mutator } from "./types";
import {
  updateSampleMutator,
  setExposureStatusMutator,
  setSelectionMutator,
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
import { patchComparesMutator } from "./mutators/compares";
import { addLabelMutator, removeLabelMutator, restoreLabelSuggestionMutator } from "./mutators/labels";
import { setLabelColorMutator } from "./mutators/labelColors";
import { patchExperimentMutator } from "./mutators/experiment";
import { patchSetupMutator } from "./mutators/setup";

export interface PersistedOpForResolution {
  kind: OpKind;
  payload: unknown;
}
export function resolveMutator(
  op: PersistedOpForResolution,
): Mutator<any, any, any> | undefined {
  switch (op.kind) {
    case "update_sample":
      return updateSampleMutator;
    case "set_exposure_status":
      return setExposureStatusMutator;
    case "set_selection":
      return setSelectionMutator;
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
      return undefined;
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
    case "patch_compares":         return patchComparesMutator;
    case "add_label":              return addLabelMutator;
    case "remove_label":           return removeLabelMutator;
    case "restore_label_suggestion": return restoreLabelSuggestionMutator;
    case "set_tag_color":          return setLabelColorMutator;
    case "patch_experiment":       return patchExperimentMutator;
    case "patch_setup":            return patchSetupMutator;
    default:
      return undefined;
  }
}
export function resolveMutatorForEvent(
  eventKind: string,
): Mutator<any, any, any> | undefined {
  switch (eventKind) {
    case "peak_added":          return peakAddMutator;
    case "peak_removed":        return peakRemoveMutator;
    case "peak_excluded":       return peakExcludeMutator;
    case "peak_unexcluded":     return peakUnexcludeMutator;
    case "speculative_created": return createSpeculativeMutator;
    case "speculative_deleted": return deleteIndexMutator;
    case "analyze_run":         return reanalyzeExposureMutator;
    case "assignment_add":         return addAssignmentPhaseMutator;
    case "assignment_remove":      return removeAssignmentPhaseMutator;
    case "sample_call_set_state":  return setSampleCallStateMutator;
    case "update_sample":       return updateSampleMutator;
    case "set_exposure_status": return setExposureStatusMutator;
    case "set_selection":       return setSelectionMutator;
    case "sample_bridged":           return bridgeSamplesMutator;
    case "sample_unbridged":         return unbridgeSampleMutator;
    case "compares_changed":         return patchComparesMutator;
    case "label_added":              return addLabelMutator;
    case "label_removed":            return removeLabelMutator;
    case "label_suggestion_restored": return restoreLabelSuggestionMutator;
    case "tag_color_set":            return setLabelColorMutator;
    case "experiment_updated":       return patchExperimentMutator;
    case "setup_updated":            return patchSetupMutator;
    default:
      return undefined;
  }
}