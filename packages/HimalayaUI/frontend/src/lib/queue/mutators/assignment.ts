import * as api from "../../../api";
import type {
  Assignment, AssignmentState, AssignmentMutationResponse, AuthOpts, SampleCall, SampleCallMutationResponse,
} from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";
import type {
  Mutator, RollbackContext, SseEvent, QueueResponseMeta, AssignmentPostState, SampleCallPostState,
} from "../types";
import { sampleFrameIds, writeSampleCall, invalidateAssignmentsExcept, invalidateSampleCalls, invalidateReviewProgress } from "../sampleCall";
import { stripQueueMetadata } from "../queueMeta";

type AssignmentScope = {
  exposureId: number;
  username: string | undefined;
  clientId: string;
};

export type AddPhaseInput = { indexId: number };
export type RemovePhaseInput = { indexId: number; exposureId?: number };
export type SetStateInput = { state: AssignmentState };

function buildAuth(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

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
    const { payload: row } = stripQueueMetadata(response);
    qc.setQueryData<Assignment>(queryKeys.assignment(p.exposureId), row as Assignment);
    invalidateAssignmentsExcept(qc, [p.exposureId]);
    invalidateReviewProgress(qc, "own");
    invalidateSampleCalls(qc);
  },
  synthesizeFromSse: synthAssignment,
};


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
    invalidateReviewProgress(qc, "own");
    invalidateSampleCalls(qc);
  },
  synthesizeFromSse: synthAssignment,
  treats404AsSuccess: true,
};


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
    invalidateAssignmentsExcept(qc, [...sampleFrameIds(qc, p.sampleId), ...(call as SampleCall).members.map((m) => m.exposure_id)]);
    invalidateReviewProgress(qc, "own");
    invalidateSampleCalls(qc);
  },
  synthesizeFromSse: synthSampleCall,
};
