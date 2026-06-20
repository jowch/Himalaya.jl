import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../src/queries";
import type { Load } from "../src/api";
import { moveExposureMutator, renameSampleMutator } from "../src/lib/queue/mutators/grouping";

function seedLoads(qc: QueryClient, experimentId: number) {
  // §8.8 shape: exposures keyed `.id`; samples carry `merged_into_id`; load keyed `load_id`.
  const loads: Load[] = [{
    load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 8, note: null,
    samples: [
      { sample_id: 10, name: "A", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
        exposures: [{ id: 100, filename: "a1.tif", horizontal_position: 8, timestamp: null }] },
      { sample_id: 20, name: "B", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
        exposures: [{ id: 200, filename: "b1.tif", horizontal_position: 12, timestamp: null }] },
    ],
  }];
  qc.setQueryData(queryKeys.loads(experimentId), loads);
}

describe("moveExposureMutator", () => {
  it("onMutate moves the exposure between samples in the loads cache; restore reverts", () => {
    const qc = new QueryClient();
    seedLoads(qc, 7);
    // exposureId + sampleId are in the Input (not scope) — one hook per experiment
    const ctx = moveExposureMutator.onMutate(
      { kind: "move_exposure", clientOpId: "op", payload: { exposureId: 100, sampleId: 20 },
        experimentId: 7, exposureId: 100, sampleId: 20,
        username: "alice", clientId: "c" } as never, qc);
    const after = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    const a = after[0]!.samples.find((s) => s.sample_id === 10)!;
    const b = after[0]!.samples.find((s) => s.sample_id === 20)!;
    expect(a.exposures.map((e) => e.id)).toEqual([]);
    expect(b.exposures.map((e) => e.id)).toEqual([200, 100]);
    ctx.restore();
    const reverted = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    expect(reverted[0]!.samples.find((s) => s.sample_id === 10)!.exposures.map((e) => e.id)).toEqual([100]);
  });

  it("onSuccess invalidates loads(experimentId) (own-op reconcile — replay never calls applyRemoteToCache)", () => {
    const qc = new QueryClient();
    seedLoads(qc, 7);
    const inv = vi.spyOn(qc, "invalidateQueries");
    moveExposureMutator.onSuccess(
      { experimentId: 7, exposureId: 100, sampleId: 20 } as never,
      { id: 100 } as never, qc);
    const keys = inv.mock.calls.map((c) => JSON.stringify((c[0] as { queryKey: unknown }).queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
  });
});

describe("renameSampleMutator", () => {
  it("onMutate rewrites the sample name + sets name_source=user; restore reverts", () => {
    const qc = new QueryClient();
    seedLoads(qc, 7);
    // sampleId is in the Input (not scope) — one hook per experiment
    const ctx = renameSampleMutator.onMutate(
      { kind: "rename_sample", clientOpId: "op", payload: { sampleId: 10, name: "Renamed" },
        experimentId: 7, sampleId: 10, name: "Renamed",
        username: "alice", clientId: "c" } as never, qc);
    const after = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    const a = after[0]!.samples.find((s) => s.sample_id === 10)!;
    expect(a.name).toBe("Renamed");
    expect(a.name_source).toBe("user");
    ctx.restore();
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples.find((s) => s.sample_id === 10)!.name).toBe("A");
  });
});
