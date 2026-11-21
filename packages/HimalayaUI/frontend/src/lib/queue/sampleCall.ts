import type { QueryClient } from "@tanstack/react-query";
import type { Assignment, AssignmentState, Exposure, SampleCallMember } from "../../api";
import { queryKeys } from "../../queries";

export function sampleFrameIds(qc: QueryClient, sampleId: number): number[] {
  return (qc.getQueryData<Exposure[]>(queryKeys.exposures(sampleId)) ?? []).map((e) => e.id);
}

export function writeSampleCall(
  qc: QueryClient,
  sampleId: number,
  call: { state: AssignmentState; members: SampleCallMember[] },
): void {
  const known = qc.getQueryData<Exposure[]>(queryKeys.exposures(sampleId)) !== undefined;
  const frames = new Set([...sampleFrameIds(qc, sampleId), ...call.members.map((m) => m.exposure_id)]);
  for (const id of frames) {
    qc.setQueryData<Assignment>(queryKeys.assignment(id), {
      exposure_id: id,
      state: call.state,
      members: call.members.filter((m) => m.exposure_id === id).map((m) => m.index_id).sort((a, b) => a - b),
    });
  }
  if (!known) {
    qc.invalidateQueries({
      predicate: (q) => q.queryKey[0] === "exposure" && q.queryKey[2] === "assignment" && !frames.has(q.queryKey[1] as number),
    });
  }
  qc.invalidateQueries({ queryKey: queryKeys.sampleCall(sampleId) });
}

function siblingFrames(qc: QueryClient, exposureId: number): number[] | undefined {
  for (const [, list] of qc.getQueriesData<Exposure[]>({
    predicate: (q) => q.queryKey[0] === "sample" && q.queryKey[2] === "exposures",
  })) {
    if (list?.some((e) => e.id === exposureId)) return list.map((e) => e.id);
  }
  return undefined;
}

export function invalidateOtherFrameAssignments(qc: QueryClient, exposureId: number): void {
  const siblings = siblingFrames(qc, exposureId);
  qc.invalidateQueries({
    predicate: (q) => q.queryKey[0] === "exposure" && q.queryKey[2] === "assignment" && q.queryKey[1] !== exposureId
      && (siblings === undefined || siblings.includes(q.queryKey[1] as number)),
  });
}

export function invalidateSampleCalls(qc: QueryClient): void {
  qc.invalidateQueries({ predicate: (q) => q.queryKey[0] === "sample" && q.queryKey[2] === "call" });
}
