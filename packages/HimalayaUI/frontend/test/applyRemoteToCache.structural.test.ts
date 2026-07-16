// test/applyRemoteToCache.structural.test.ts
import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../src/lib/queue/applyRemoteToCache";
import { queryKeys } from "../src/queries";
import type { SseEvent } from "../src/lib/queue/types";
import type { Sample } from "../src/api";

function frame(kind: string, entity_id: number, payload: Record<string, unknown>): SseEvent {
  return { kind, entity_id, entity_type: "sample", payload } as unknown as SseEvent;
}

describe("structural grouping SSE receive (Phase E1)", () => {
  it("sample_renamed splices the new name into the sample cache + invalidates listings", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    qc.setQueryData<Sample>(queryKeys.sample(5), {
      id: 5, experiment_id: 7, name: "Old", notes: null, tags: [],
    });
    applyRemoteToCache(frame("sample_renamed", 5, { name: "HA85 (S01P15)", experiment_id: 7 }), qc);
    expect(qc.getQueryData<Sample>(queryKeys.sample(5))?.name).toBe("HA85 (S01P15)");
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
  });

  it("exposure_moved invalidates both samples' exposures + the loads roll-up", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    // Contract payload: sample_id = destination, from_sample_id = source.
    applyRemoteToCache(
      frame("exposure_moved", 99, { sample_id: 2, from_sample_id: 1, experiment_id: 7 }),
      qc,
    );
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.exposures(1)));
    expect(keys).toContain(JSON.stringify(queryKeys.exposures(2)));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
  });

  it("sample_created invalidates loads + sample listings, NEVER peaks/indices (precedes default)", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(frame("sample_created", 42, { experiment_id: 7 }), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.corpusSamples));
    // entity_id 42 is a sample id — the default arm would poison these:
    expect(keys).not.toContain(JSON.stringify(queryKeys.peaks(42)));
    expect(keys).not.toContain(JSON.stringify(queryKeys.indices(42)));
  });

  it("sample_split / grouping_flag_dismissed refetch loads + sample listings (invalidate-only)", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(frame("sample_split", 2, { new_sample_id: 9, exposure_ids: [1, 2], experiment_id: 7 }), qc);
    applyRemoteToCache(frame("grouping_flag_dismissed", 2, { flag_kind: "merge", merge_with_sample_id: 4, experiment_id: 7 }), qc);
    const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.samples(7)));
    expect(keys).toContain(JSON.stringify(queryKeys.corpusSamples));
  });
});
