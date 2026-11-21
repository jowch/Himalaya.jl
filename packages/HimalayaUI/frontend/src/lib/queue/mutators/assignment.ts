/**
 * Assignment-cart mutators (Plan D Task D-3):
 *
 * - addAssignmentPhaseMutator    (kind: "assignment_add")
 * - removeAssignmentPhaseMutator (kind: "assignment_remove")
 * - setSampleCallStateMutator    (kind: "sample_call_set_state")
 *
 * All three write `queryKeys.assignment` (the third for every frame of the sample)
 * and declare `affectsExposurePeaks: () => false` (they never touch the peak set).
 *
 * Two load-bearing invariants from the queue review:
 *
 *   1. onSuccess writes ONLY queryKeys.assignment — NEVER queryKeys.exposure.
 *      The framework's synthesizeResponseFromSse lifts analysis_inputs_hash
 *      from post_state, which is `undefined` for an assignment frame (the
 *      {assignment} post_state carries no hash). Writing the exposure cache
 *      here would wipe its analysis_inputs_hash with undefined.
 *
 *   2. The rollback snapshots the WHOLE Assignment object, so it is symmetric
 *      across both `members` and `state` mutations — members are existing
 *      positive index ids (no negative-optimistic-id concern).
 *
 * The SSE-confirmation arm lives in applyRemoteToCache (`assignment_*` case),
 * keyed on the DISTINCT `{assignment}` post_state shape.
 */
import * as api from "../../../api";
import type {
  Assignment, AssignmentState, AssignmentMutationResponse, AuthOpts, SampleCall, SampleCallMutationResponse,
} from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type {
  Mutator, RollbackContext, SseEvent, QueueResponseMeta, AssignmentPostState, SampleCallPostState,
} from "../types";
import { sampleFrameIds, writeSampleCall, invalidateOtherFrameAssignments } from "../sampleCall";
import { stripQueueMetadata } from "../queueMeta";

type AssignmentScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};

export type AddPhaseInput = { indexId: number };
export type RemovePhaseInput = { indexId: number };
export type SetStateInput = { state: AssignmentState };

function buildAuth(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

/**
 * Build the SSE-wins synthetic response from an assignment frame's distinct
 * `{assignment:{state,members}}` post_state. Shared by all three mutators —
 * the response shape (an Assignment + queue metadata) is identical regardless
 * of which gesture produced it. Falls back to `undefined` (generic shape) when
 * the frame lacks the {assignment} envelope.
 */
function synthAssignment(remote: SseEvent, base: QueueResponseMeta): AssignmentMutationResponse | undefined {
  const ps = remote.post_state as AssignmentPostState | undefined;
  if (!ps?.assignment) return undefined;
  return {
    ...base,
    view_row_id: base.view_row_id ?? null,
    exposure_id: remote.entity_id,
    state: ps.assignment.state,
    members: ps.assignment.members,
  };
}

// ---------------------------------------------------------------------------
// addAssignmentPhaseMutator
// ---------------------------------------------------------------------------

export const addAssignmentPhaseMutator: Mutator<AddPhaseInput, AssignmentScope, AssignmentMutationResponse> = {
  kind: "assignment_add",
  onMutate: (p, qc): RollbackContext => {
    const key = queryKeys.assignment(p.exposureId);
    const prev = qc.getQueryData<Assignment>(key);
    if (prev) {
      qc.setQueryData<Assignment>(key, {
        ...prev,
        state: "indexed",
        members: prev.members.includes(p.indexId)
          ? prev.members
          : [...prev.members, p.indexId],
      });
    }
    return {
      restore: () => { if (prev !== undefined) qc.setQueryData(key, prev); },
    };
  },
  request: (p) => api.addAssignmentPhase(p.exposureId, p.indexId, buildAuth(p)),
  onSuccess: (p, response, qc) => {
    // Strip queue metadata; write ONLY the assignment cache (never exposure).
    const { payload: row } = stripQueueMetadata(response);
    qc.setQueryData<Assignment>(queryKeys.assignment(p.exposureId), row as Assignment);
    invalidateOtherFrameAssignments(qc, p.exposureId);
    qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
  },
  synthesizeFromSse: synthAssignment,
  affectsExposurePeaks: () => false,
};

// ---------------------------------------------------------------------------
// removeAssignmentPhaseMutator
// ---------------------------------------------------------------------------

export const removeAssignmentPhaseMutator: Mutator<RemovePhaseInput, AssignmentScope, AssignmentMutationResponse> = {
  kind: "assignment_remove",
  onMutate: (p, qc): RollbackContext => {
    const key = queryKeys.assignment(p.exposureId);
    const prev = qc.getQueryData<Assignment>(key);
    if (prev) {
      qc.setQueryData<Assignment>(key, {
        ...prev,
        members: prev.members.filter((m) => m !== p.indexId),
      });
    }
    return {
      restore: () => { if (prev !== undefined) qc.setQueryData(key, prev); },
    };
  },
  request: (p) => api.removeAssignmentPhase(p.exposureId, p.indexId, buildAuth(p)),
  onSuccess: (p, response, qc) => {
    const { payload: row } = stripQueueMetadata(response);
    qc.setQueryData<Assignment>(queryKeys.assignment(p.exposureId), row as Assignment);
    qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
  },
  synthesizeFromSse: synthAssignment,
  affectsExposurePeaks: () => false,
  // 404 = the member isn't in the assignment → desired end state.
  treats404AsSuccess: true,
};

// ---------------------------------------------------------------------------
// setSampleCallStateMutator
// ---------------------------------------------------------------------------

type SampleCallScope = {
  sampleId: number;
  username: string | undefined;
  clientId: string;
};

function synthSampleCall(remote: SseEvent, base: QueueResponseMeta): SampleCallMutationResponse | undefined {
  const ps = remote.post_state as SampleCallPostState | undefined;
  if (!ps?.call) return undefined;
  return {
    ...base,
    view_row_id: base.view_row_id ?? null,
    sample_id: remote.entity_id,
    state: ps.call.state,
    phase: null,
    call_phases: [],
    members: ps.call.members,
  };
}

export const setSampleCallStateMutator: Mutator<SetStateInput, SampleCallScope, SampleCallMutationResponse> = {
  kind: "sample_call_set_state",
  onMutate: (p, qc): RollbackContext => {
    const prev = sampleFrameIds(qc, p.sampleId).map((id) => {
      const key = queryKeys.assignment(id);
      return [key, qc.getQueryData<Assignment>(key)] as const;
    });
    for (const [key, a] of prev) {
      if (a) qc.setQueryData<Assignment>(key, { ...a, state: p.state, members: p.state === "indexed" ? a.members : [] });
    }
    return {
      restore: () => { for (const [key, a] of prev) if (a !== undefined) qc.setQueryData(key, a); },
    };
  },
  request: (p) => api.setSampleCallState(p.sampleId, p.state, buildAuth(p)),
  onSuccess: (p, response, qc) => {
    const { payload: call } = stripQueueMetadata(response);
    writeSampleCall(qc, p.sampleId, call as SampleCall);
    qc.invalidateQueries({ queryKey: queryKeys.corpusSamples });
  },
  synthesizeFromSse: synthSampleCall,
  affectsExposurePeaks: () => false,
};
