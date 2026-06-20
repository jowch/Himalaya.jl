import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../src/queries";
import type { Load } from "../src/api";
import { dismissGroupingFlagMutator } from "../src/lib/queue/mutators/grouping";

function seed(qc: QueryClient) {
  const loads: Load[] = [{
    load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [{ sample_id: 20, name: "B", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
      flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "A" }, exposures: [] }],
  }];
  qc.setQueryData(queryKeys.loads(7), loads);
}

describe("dismissGroupingFlagMutator (durable 'Keep separate')", () => {
  it("onMutate clears sample.flag optimistically; restore reverts", () => {
    const qc = new QueryClient(); seed(qc);
    // sampleId is in the Input (not scope) — one hook per experiment, not per row
    const ctx = dismissGroupingFlagMutator.onMutate(
      { kind: "dismiss_grouping_flag", clientOpId: "op",
        payload: { sampleId: 20, flagKind: "merge", mergeWithSampleId: 10 },
        experimentId: 7, sampleId: 20, flagKind: "merge", mergeWithSampleId: 10,
        username: "a", clientId: "c" } as never, qc);
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples[0]!.flag).toBeNull();
    ctx.restore();
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples[0]!.flag).not.toBeNull();
  });
  it("onSuccess invalidates loads(experimentId)", () => {
    const qc = new QueryClient(); seed(qc);
    const inv = vi.spyOn(qc, "invalidateQueries");
    dismissGroupingFlagMutator.onSuccess({ experimentId: 7, sampleId: 20 } as never, undefined as never, qc);
    const keys = inv.mock.calls.map((c) => JSON.stringify((c[0] as { queryKey: unknown }).queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
  });
});
