import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import { bridgeSamplesMutator, unbridgeSampleMutator } from "../../src/lib/queue/mutators/bridges";
import { queryKeys } from "../../src/queries";
import type { SseEvent } from "../../src/lib/queue/types";
import type { SampleBridge } from "../../src/api";

function frame(kind: string, entity_id: number, payload: Record<string, unknown>): SseEvent {
  return { kind, entity_id, entity_type: "sample", payload } as unknown as SseEvent;
}

const scope = { username: "jc", clientId: "c", clientOpId: "op", kind: "bridge_samples", experimentId: 7 } as const;

describe("bridge mutations", () => {
  it("a remote bridge event refetches that experiment's bridges and nothing exposure-scoped", () => {
    for (const kind of ["sample_bridged", "sample_unbridged"]) {
      const qc = new QueryClient();
      const spy = vi.spyOn(qc, "invalidateQueries");
      applyRemoteToCache(frame(kind, 42, { group_id: 3, experiment_id: 7 }), qc);
      const keys = spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey));
      expect(keys).toEqual([JSON.stringify(queryKeys.bridges(7))]);
    }
  });

  it("routes both the outbound op and the echo back to the bridge mutators", () => {
    expect(resolveMutator({ kind: "bridge_samples", payload: {} })).toBe(bridgeSamplesMutator);
    expect(resolveMutator({ kind: "unbridge_sample", payload: {} })).toBe(unbridgeSampleMutator);
    expect(resolveMutatorForEvent("sample_bridged")).toBe(bridgeSamplesMutator);
    expect(resolveMutatorForEvent("sample_unbridged")).toBe(unbridgeSampleMutator);
  });

  it("bridging shows at once and rolls back cleanly", () => {
    const qc = new QueryClient();
    const before: SampleBridge[] = [{ sample_id: 1, group_id: 5 }, { sample_id: 2, group_id: 5 }];
    qc.setQueryData(queryKeys.bridges(7), before);
    const ctx = bridgeSamplesMutator.onMutate({ ...scope, sampleId: 3, toSampleId: 1 }, qc);
    expect(qc.getQueryData<SampleBridge[]>(queryKeys.bridges(7))?.map((r) => r.sample_id).sort()).toEqual([1, 2, 3]);
    ctx.restore();
    expect(qc.getQueryData(queryKeys.bridges(7))).toEqual(before);
  });

  it("takes the server's list when the response is one, and refetches when it is a partial from the live stream", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    const server: SampleBridge[] = [{ sample_id: 1, group_id: 9 }, { sample_id: 4, group_id: 9 }];
    unbridgeSampleMutator.onSuccess({ ...scope, kind: "unbridge_sample", sampleId: 2 }, server, qc);
    expect(qc.getQueryData(queryKeys.bridges(7))).toEqual(server);
    unbridgeSampleMutator.onSuccess({ ...scope, kind: "unbridge_sample", sampleId: 2 }, { id: 1 } as unknown as SampleBridge[], qc);
    expect(qc.getQueryData(queryKeys.bridges(7))).toEqual(server);
    expect(spy).toHaveBeenCalledWith({ queryKey: queryKeys.bridges(7) });
  });
});
