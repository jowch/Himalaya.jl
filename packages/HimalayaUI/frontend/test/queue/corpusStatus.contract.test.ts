import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import {
  addAssignmentPhaseMutator,
  removeAssignmentPhaseMutator,
  setSampleCallStateMutator,
} from "../../src/lib/queue/mutators/assignment";
import { setExposureStatusMutator } from "../../src/lib/queue/mutators/trivial";
import type { SseEvent } from "../../src/lib/queue/types";
import type { CorpusSample } from "../../src/api";
import { queryKeys } from "../../src/queries";

let qc: QueryClient;

function corpusStale(): boolean {
  return qc.getQueryState(queryKeys.corpusSamples)?.isInvalidated ?? false;
}

beforeEach(() => {
  qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  qc.setQueryData<CorpusSample[]>(queryKeys.corpusSamples, []);
  expect(corpusStale()).toBe(false);
});

const scope = { exposureId: 5, sampleId: 1, username: "alice", clientId: "tab-1", clientOpId: "op-1" };
const assignmentRow = { exposure_id: 5, state: "indexed", members: [42], event_id: 1, view_row_id: 1 };

describe("corpus sample status follows the calls that derive it", () => {
  describe("own op: mutator onSuccess", () => {
    it.each([
      ["assignment_add", () => addAssignmentPhaseMutator.onSuccess({ ...scope, kind: "assignment_add", indexId: 42 }, assignmentRow as never, qc)],
      ["assignment_remove", () => removeAssignmentPhaseMutator.onSuccess({ ...scope, kind: "assignment_remove", indexId: 42 }, assignmentRow as never, qc)],
      ["sample_call_set_state", () => setSampleCallStateMutator.onSuccess({ ...scope, kind: "sample_call_set_state", state: "form_factor" }, { sample_id: 1, state: "form_factor", members: [] } as never, qc)],
      ["set_exposure_status", () => setExposureStatusMutator.onSuccess({ ...scope, status: "rejected" } as never, { id: 5, status: "rejected" }, qc)],
    ])("%s invalidates the corpus samples", (_kind, run) => {
      run();
      expect(corpusStale()).toBe(true);
    });
  });

  describe("foreign op: applyRemoteToCache", () => {
    const frames: SseEvent[] = [
      { id: 9, kind: "assignment_add", entity_type: "exposure", entity_id: 5, payload: { index_id: 42 },
        post_state: { assignment: { state: "indexed", members: [42] } } },
      { id: 9, kind: "assignment_remove", entity_type: "exposure", entity_id: 5, payload: { index_id: 42 },
        post_state: { assignment: { state: "indexed", members: [] } } },
      { id: 9, kind: "sample_call_set_state", entity_type: "sample", entity_id: 1, payload: { state: "form_factor" },
        post_state: { call: { state: "form_factor", members: [] } } },
      { id: 9, kind: "set_exposure_status", entity_type: "exposure", entity_id: 5, payload: { status: "rejected" } },
    ];
    it.each(frames.map((f) => [f.kind, f] as const))("%s invalidates the corpus samples", (_kind, frame) => {
      applyRemoteToCache(frame, qc);
      expect(corpusStale()).toBe(true);
    });
  });
});
