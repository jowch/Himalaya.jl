import { describe, it, expect, beforeEach, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { handleRemoteEvent } from "../../src/lib/queue/replayCoordinator";
import { makeDeferred, getDeferred, pendingDeferreds } from "../../src/lib/queue/deferred";
import { makeFakeMutation, remoteForeignEvent } from "./helpers";

describe("handleRemoteEvent", () => {
  let qc: QueryClient;
  beforeEach(() => {
    pendingDeferreds.clear();
    qc = new QueryClient();
  });

  it("SSE with matching client_op_id resolves the pending deferred and clears the registry", async () => {
    const d = makeDeferred<{ event_id: number; client_op_id: string }>("op-1");
    handleRemoteEvent({
      id: 7,
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: "op-1",
      payload: { q: 1.0 },
    }, qc, qc.getMutationCache());
    const result = await d.promise;
    expect(result).toMatchObject({ event_id: 7, client_op_id: "op-1" });
    expect(getDeferred("op-1")).toBeUndefined();
  });

  it("SSE without matching client_op_id triggers rollback-then-rerun on pending mutations", () => {
    const restore = vi.fn();
    const onMutate = vi.fn();
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      context: { restore },
      variables: { kind: "peak_excluded", exposureId: 42, q: 2.0, clientOpId: "mine-1" },
      onMutate,
    }));
    handleRemoteEvent(remoteForeignEvent({ entity_id: 42 }), qc, qc.getMutationCache());
    expect(restore).toHaveBeenCalledTimes(1);
    expect(onMutate).toHaveBeenCalledTimes(1);
  });

  it("Rollback iterates pending mutations in REVERSE order", () => {
    const order: number[] = [];
    for (let i = 0; i < 3; i++) {
      qc.getMutationCache().add(makeFakeMutation({
        status: "pending",
        context: { restore: () => order.push(i) },
        onMutate: () => {},
      }));
    }
    handleRemoteEvent(remoteForeignEvent(), qc, qc.getMutationCache());
    expect(order).toEqual([2, 1, 0]);
  });

  it("Re-run iterates pending mutations in INSERTION order", () => {
    const order: number[] = [];
    for (let i = 0; i < 3; i++) {
      qc.getMutationCache().add(makeFakeMutation({
        status: "pending",
        context: { restore: () => {} },
        onMutate: () => order.push(i),
      }));
    }
    handleRemoteEvent(remoteForeignEvent(), qc, qc.getMutationCache());
    expect(order).toEqual([0, 1, 2]);
  });

  it("MutationCache.getAll() preserves insertion order (load-bearing TanStack invariant)", () => {
    const ids = ["a", "b", "c", "d", "e"];
    for (const id of ids) {
      qc.getMutationCache().add(makeFakeMutation({
        status: "pending",
        context: {},
        mutationKey: [id],
      }));
    }
    const observed = qc.getMutationCache().getAll().map(
      m => (m.options.mutationKey as readonly unknown[])[0]
    );
    expect(observed).toEqual(ids);
  });

  it("Two foreign events trigger restore twice (one per event), not four times", () => {
    const restore = vi.fn();
    qc.getMutationCache().add(makeFakeMutation({
      status: "pending",
      context: { restore },
      onMutate: () => {},
    }));
    handleRemoteEvent(remoteForeignEvent(), qc, qc.getMutationCache());
    handleRemoteEvent(remoteForeignEvent(), qc, qc.getMutationCache());
    expect(restore).toHaveBeenCalledTimes(2);
  });

  it("Non-pending mutations (success/error) are skipped during rollback and re-run", () => {
    const restore = vi.fn();
    const onMutate = vi.fn();
    qc.getMutationCache().add(makeFakeMutation({
      status: "success",
      context: { restore },
      onMutate,
    }));
    handleRemoteEvent(remoteForeignEvent(), qc, qc.getMutationCache());
    expect(restore).not.toHaveBeenCalled();
    expect(onMutate).not.toHaveBeenCalled();
  });

  it("applyRemoteToCache for peak_added appends an optimistic-id row to the peaks cache", () => {
    qc.setQueryData(["exposure", 42, "peaks"], [
      { id: 1, q: 0.5, kind: "auto", excluded: false },
    ]);
    handleRemoteEvent({
      id: 99,
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: "foreign-add",
      payload: { q: 1.7 },
    }, qc, qc.getMutationCache());
    const peaks = qc.getQueryData(["exposure", 42, "peaks"]) as any[];
    expect(peaks).toHaveLength(2);
    expect(peaks[1]).toMatchObject({ q: 1.7, kind: "add", excluded: false });
    expect(peaks[1].id).toBeLessThan(0);
  });

  it("applyRemoteToCache for peak_excluded sets excluded=true on the matching peak", () => {
    qc.setQueryData(["exposure", 42, "peaks"], [
      { id: 1, q: 0.5, kind: "auto", excluded: false },
      { id: 2, q: 1.7, kind: "auto", excluded: false },
    ]);
    handleRemoteEvent({
      id: 100,
      kind: "peak_excluded",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: null,
      payload: { q: 1.7 },
    }, qc, qc.getMutationCache());
    const peaks = qc.getQueryData(["exposure", 42, "peaks"]) as any[];
    expect(peaks.find(p => p.q === 1.7)?.excluded).toBe(true);
    expect(peaks.find(p => p.q === 0.5)?.excluded).toBe(false);
  });

  it("applyRemoteToCache for analyze_run with post_state writes the indices array and updates exposure.analysis_inputs_hash", () => {
    qc.setQueryData(["exposure-entity", 42], { id: 42, analysis_inputs_hash: "old-hash" });
    qc.setQueryData(["exposure", 42, "indices"], [{ id: 1 }]);
    handleRemoteEvent({
      id: 200,
      kind: "analyze_run",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: null,
      post_state: { analysis_inputs_hash: "new-hash", indices: [{ id: 5 }, { id: 6 }] },
    }, qc, qc.getMutationCache());
    expect(qc.getQueryData(["exposure", 42, "indices"])).toEqual([{ id: 5 }, { id: 6 }]);
    expect((qc.getQueryData(["exposure-entity", 42]) as any).analysis_inputs_hash).toBe("new-hash");
  });

  it("applyRemoteToCache default branch invalidates peaks/indices/groups for unknown kinds", () => {
    const spy = vi.spyOn(qc, "invalidateQueries");
    handleRemoteEvent({
      id: 300,
      kind: "weird_unknown_kind",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: null,
    }, qc, qc.getMutationCache());
    expect(spy).toHaveBeenCalledTimes(3);
  });
});
