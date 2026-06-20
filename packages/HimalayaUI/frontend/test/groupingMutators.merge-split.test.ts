import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../src/queries";
import type { Load } from "../src/api";
import { mergeSamplesMutator, splitSampleMutator } from "../src/lib/queue/mutators/grouping";

function seed(qc: QueryClient) {
  const loads: Load[] = [{
    load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [
      { sample_id: 10, name: "survivor", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
        exposures: [{ id: 100, filename: "a.tif", horizontal_position: 8, timestamp: null }] },
      { sample_id: 20, name: "loser", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
        flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "survivor" },
        exposures: [{ id: 200, filename: "b.tif", horizontal_position: 12, timestamp: null }] },
    ],
  }];
  qc.setQueryData(queryKeys.loads(7), loads);
}

describe("mergeSamplesMutator", () => {
  it("onMutate re-points loser exposures onto survivor and removes the loser; restore reverts", () => {
    const qc = new QueryClient(); seed(qc);
    const ctx = mergeSamplesMutator.onMutate(
      { kind: "merge_samples", clientOpId: "op", payload: { loserId: 20, survivorId: 10 },
        experimentId: 7, loserId: 20, survivorId: 10, username: "a", clientId: "c" } as never, qc);
    const after = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    expect(after[0]!.samples.map((s) => s.sample_id)).toEqual([10]);
    expect(after[0]!.samples[0]!.exposures.map((e) => e.id)).toEqual([100, 200]);
    ctx.restore();
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples.map((s) => s.sample_id)).toEqual([10, 20]);
  });
  it("onSuccess invalidates loads(experimentId)", () => {
    const qc = new QueryClient(); seed(qc);
    const inv = vi.spyOn(qc, "invalidateQueries");
    mergeSamplesMutator.onSuccess({ experimentId: 7, loserId: 20, survivorId: 10 } as never, { loser_id: 20, survivor_id: 10 } as never, qc);
    const keys = inv.mock.calls.map((c) => JSON.stringify((c[0] as { queryKey: unknown }).queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
  });
});

describe("splitSampleMutator", () => {
  it("onMutate creates a new placeholder sample with the chosen exposures; restore reverts", () => {
    const qc = new QueryClient(); seed(qc);
    // add a second exposure to sample 10 so we can split it
    const loads = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    loads[0]!.samples[0]!.exposures.push({ id: 101, filename: "a2.tif", horizontal_position: 30, timestamp: null });
    qc.setQueryData(queryKeys.loads(7), loads);

    // sampleId is in the Input (not scope) — one hook per experiment, not per row
    const ctx = splitSampleMutator.onMutate(
      { kind: "split_sample", clientOpId: "op", payload: { sampleId: 10, exposureIds: [101], name: "survivorb" },
        experimentId: 7, sampleId: 10, exposureIds: [101], name: "survivorb",
        username: "a", clientId: "c" } as never, qc);
    const after = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    const src = after[0]!.samples.find((s) => s.sample_id === 10)!;
    expect(src.exposures.map((e) => e.id)).toEqual([100]);
    const created = after[0]!.samples.find((s) => s.sample_id < 0)!;
    expect(created.name).toBe("survivorb");
    expect(created.exposures.map((e) => e.id)).toEqual([101]);
    ctx.restore();
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples.every((s) => s.sample_id > 0)).toBe(true);
  });
});
