import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { resolveMutator, resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import { patchComparesMutator } from "../../src/lib/queue/mutators/compares";
import { applyCompareChange, comparesOf, EMPTY_COMPARES } from "../../src/lib/compares";
import { queryKeys } from "../../src/queries";
import type { SseEvent } from "../../src/lib/queue/types";
import type { CompareState } from "../../src/api";

const scope = { username: "jc", clientId: "c", clientOpId: "op", kind: "patch_compares", payload: {}, experimentId: 7 } as const;
const base: CompareState = {
  compares: [
    { sample_id: 1, ref_sample_id: 2, added_by: "ana" },
    { sample_id: 1, ref_sample_id: 3, added_by: "jc" },
    { sample_id: 4, ref_sample_id: 2, added_by: "jc" },
  ],
  subtractions: [{ sample_id: 1, ref_sample_id: 2 }],
};

describe("compare changes", () => {
  it("adds after the existing ones, skips repeats, and stamps who added them", () => {
    const next = applyCompareChange(base, 1, { add: [3, 5, 5] }, "jc");
    expect(comparesOf(next, 1).refs).toEqual([2, 3, 5]);
    expect(comparesOf(next, 1).addedBy.get(5)).toBe("jc");
    expect(comparesOf(next, 4).refs).toEqual([2]);
  });

  it("removing the subtracted comparison clears the subtraction unless the change sets one", () => {
    expect(comparesOf(applyCompareChange(base, 1, { remove: [2] }, "jc"), 1)).toMatchObject({ refs: [3], subtraction: undefined });
    expect(comparesOf(applyCompareChange(base, 1, { remove: [3] }, "jc"), 1).subtraction).toBe(2);
    expect(comparesOf(applyCompareChange(base, 1, { remove: [2], subtraction: 3 }, "jc"), 1).subtraction).toBe(3);
    expect(comparesOf(applyCompareChange(base, 1, { subtraction: null }, "jc"), 1).subtraction).toBeUndefined();
    expect(comparesOf(applyCompareChange(EMPTY_COMPARES, 9, { subtraction: 1 }, null), 9).subtraction).toBe(1);
  });

  it("a saved subtraction that is not a comparison is cleared by the next change, so adding it never makes it the background", () => {
    const stale = { compares: [{ sample_id: 1, ref_sample_id: 3, added_by: "jc" }], subtractions: [{ sample_id: 1, ref_sample_id: 2 }] };
    expect(comparesOf(applyCompareChange(stale, 1, { add: [2] }, "jc"), 1)).toMatchObject({ refs: [3, 2], subtraction: undefined });
    expect(comparesOf(applyCompareChange(stale, 1, { remove: [3] }, "jc"), 1).subtraction).toBeUndefined();
    expect(comparesOf(applyCompareChange(stale, 1, { remove: [9] }, "jc"), 1).subtraction).toBe(2);
    expect(comparesOf(applyCompareChange(stale, 1, { add: [2], subtraction: 2 }, "jc"), 1).subtraction).toBe(2);
  });
});

describe("compare mutations", () => {
  it("a remote change refetches that experiment's comparisons only", () => {
    const qc = new QueryClient();
    const spy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache({ kind: "compares_changed", entity_id: 1, entity_type: "sample", payload: { add: [2], experiment_id: 7 } } as unknown as SseEvent, qc);
    expect(spy.mock.calls.map((c) => JSON.stringify(c[0]?.queryKey))).toEqual([JSON.stringify(queryKeys.compares(7))]);
  });

  it("routes the outbound op and its echo to the compare mutator", () => {
    expect(resolveMutator({ kind: "patch_compares", payload: {} })).toBe(patchComparesMutator);
    expect(resolveMutatorForEvent("compares_changed")).toBe(patchComparesMutator);
  });

  it("shows at once, rolls back cleanly, and takes the server's state when the response is one", () => {
    const qc = new QueryClient();
    qc.setQueryData(queryKeys.compares(7), base);
    const ctx = patchComparesMutator.onMutate({ ...scope, sampleId: 1, change: { add: [9] } }, qc);
    expect(comparesOf(qc.getQueryData<CompareState>(queryKeys.compares(7))!, 1).refs).toEqual([2, 3, 9]);
    ctx.restore();
    expect(qc.getQueryData(queryKeys.compares(7))).toEqual(base);
    const server: CompareState = { compares: [], subtractions: [] };
    patchComparesMutator.onSuccess({ ...scope, sampleId: 1, change: {} }, server, qc);
    expect(qc.getQueryData(queryKeys.compares(7))).toEqual(server);
    patchComparesMutator.onSuccess({ ...scope, sampleId: 1, change: {} }, { id: 3 } as unknown as CompareState, qc);
    expect(qc.getQueryData(queryKeys.compares(7))).toEqual(server);
  });
});
