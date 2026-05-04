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

  it("applyRemoteToCache for peak_added inserts the row using server-assigned peak_curation_id", () => {
    qc.setQueryData(["exposure", 42, "peaks"], [
      { id: 1, q: 0.5, kind: "auto", excluded: false },
    ]);
    handleRemoteEvent({
      id: 99, // event id (NOT used as peak id post issue #2 fix)
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: "foreign-add",
      payload: { q: 1.7, peak_curation_id: 314 },
    }, qc, qc.getMutationCache());
    const peaks = qc.getQueryData(["exposure", 42, "peaks"]) as any[];
    expect(peaks).toHaveLength(2);
    expect(peaks[1]).toMatchObject({
      id: 314,
      q: 1.7,
      source: "manual",
      excluded: false,
    });
  });

  it("peak_added without peak_curation_id falls back to invalidation rather than inserting a phantom row", () => {
    qc.setQueryData(["exposure", 42, "peaks"], [
      { id: 1, q: 0.5, kind: "auto", excluded: false },
    ]);
    handleRemoteEvent({
      id: 99,
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: "foreign-add-legacy",
      payload: { q: 1.7 }, // no peak_curation_id
    }, qc, qc.getMutationCache());
    // Cache untouched; invalidation has been queued (no phantom insert).
    const peaks = qc.getQueryData(["exposure", 42, "peaks"]) as any[];
    expect(peaks).toHaveLength(1);
  });

  it("peak_removed filters by peak_curation_id (not q) — preserves untargeted peaks even at the same q", () => {
    qc.setQueryData(["exposure", 42, "peaks"], [
      { id: 1, q: 0.5, source: "auto", excluded: false },
      { id: 2, q: 1.7, source: "manual", excluded: false },
      { id: 3, q: 1.7, source: "manual", excluded: false },
    ]);
    handleRemoteEvent({
      id: 100,
      kind: "peak_removed",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: null,
      payload: { peak_curation_id: 2 },
    }, qc, qc.getMutationCache());
    const peaks = qc.getQueryData(["exposure", 42, "peaks"]) as any[];
    // Only peak id=2 removed; the same-q peak id=3 survives because we filter
    // by id, not q (issue #1 from PR review — pre-fix behavior wiped both).
    expect(peaks.map(p => p.id).sort()).toEqual([1, 3]);
  });

  it("peak_excluded prefers auto_peak_id over q tolerance (suggestion #1)", () => {
    qc.setQueryData(["exposure", 42, "peaks"], [
      { id: 10, q: 1.7000, source: "auto", excluded: false },
      { id: 11, q: 1.7005, source: "auto", excluded: false }, // within tol of 1.7
    ]);
    handleRemoteEvent({
      id: 200,
      kind: "peak_excluded",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: null,
      payload: { q: 1.7000, auto_peak_id: 11 },
    }, qc, qc.getMutationCache());
    const peaks = qc.getQueryData(["exposure", 42, "peaks"]) as any[];
    expect(peaks.find(p => p.id === 10)?.excluded).toBe(false);
    expect(peaks.find(p => p.id === 11)?.excluded).toBe(true);
  });

  it("self-echo SSE for an own op with no matching deferred is dropped (issue #8)", () => {
    qc.setQueryData(["exposure", 42, "peaks"], [
      { id: 1, q: 0.5, source: "auto", excluded: false },
    ]);
    // Read the per-tab client_id the way the runtime does so we can spoof it.
    const ourClientId = (() => {
      const k = "himalaya.client_id";
      let v = sessionStorage.getItem(k);
      if (!v) { v = "test-client"; sessionStorage.setItem(k, v); }
      return v;
    })();
    handleRemoteEvent({
      id: 99,
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_id: ourClientId,
      client_op_id: "own-op-no-deferred",
      payload: { q: 1.7, peak_curation_id: 999 },
    }, qc, qc.getMutationCache());
    const peaks = qc.getQueryData(["exposure", 42, "peaks"]) as any[];
    // Self-echo guard short-circuits Case 2; cache is unchanged.
    expect(peaks).toHaveLength(1);
  });

  it("aborts the HTTP request when SSE resolves the deferred first", () => {
    const controller = new AbortController();
    const d = makeDeferred<unknown>("op-abort-test");
    d.controller = controller;
    handleRemoteEvent({
      id: 1,
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: "op-abort-test",
    }, qc, qc.getMutationCache());
    expect(controller.signal.aborted).toBe(true);
  });

  it("applyRemoteToCache for peak_excluded sets excluded=true on the matching peak", () => {
    qc.setQueryData(["exposure", 42, "peaks"], [
      { id: 1, q: 0.5, source: "auto", excluded: false },
      { id: 2, q: 1.7, source: "auto", excluded: false },
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

  it("own-op SSE (matching deferred) applies post_state to indices and exposure caches", async () => {
    qc.setQueryData(["exposure-entity", 42], { id: 42, analysis_inputs_hash: "old-hash" });
    qc.setQueryData(["exposure", 42, "indices"], [{ id: 1, inputs_hash: "old-hash" }]);
    const d = makeDeferred<{ id?: number }>("own-op-with-postState");
    handleRemoteEvent({
      id: 7,
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: "own-op-with-postState",
      payload: { q: 1.7, peak_curation_id: 314 },
      post_state: { analysis_inputs_hash: "new-hash", indices: [{ id: 5, inputs_hash: "new-hash" }] },
    }, qc, qc.getMutationCache());
    await d.promise;
    // Without this propagation the indices stay frozen at "old-hash" and the
    // StaleIndicesBanner sticks until a hard refetch.
    expect(qc.getQueryData(["exposure", 42, "indices"])).toEqual([{ id: 5, inputs_hash: "new-hash" }]);
    expect((qc.getQueryData(["exposure-entity", 42]) as any).analysis_inputs_hash).toBe("new-hash");
  });

  it("self-echo SSE (own client_id, no deferred) still applies post_state to indices/exposure", () => {
    qc.setQueryData(["exposure-entity", 42], { id: 42, analysis_inputs_hash: "old-hash" });
    qc.setQueryData(["exposure", 42, "indices"], [{ id: 1, inputs_hash: "old-hash" }]);
    const ourClientId = (() => {
      const k = "himalaya.client_id";
      let v = sessionStorage.getItem(k);
      if (!v) { v = "test-client"; sessionStorage.setItem(k, v); }
      return v;
    })();
    handleRemoteEvent({
      id: 99,
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_id: ourClientId,
      client_op_id: "self-echo-no-deferred",
      payload: { q: 1.7, peak_curation_id: 999 },
      post_state: { analysis_inputs_hash: "new-hash", indices: [{ id: 5, inputs_hash: "new-hash" }] },
    }, qc, qc.getMutationCache());
    // Self-echo guard short-circuits the per-kind body (mutator already wrote
    // peaks + exposure hash), but post_state must still propagate to indices —
    // mutator onSuccess paths don't write indices, so without this the banner
    // sticks. Regression for "stale banner persists after own peak add."
    expect(qc.getQueryData(["exposure", 42, "indices"])).toEqual([{ id: 5, inputs_hash: "new-hash" }]);
    expect((qc.getQueryData(["exposure-entity", 42]) as any).analysis_inputs_hash).toBe("new-hash");
  });

  it("synthesizeResponseFromSse for peak_added produces a properly-shaped Peak (id, source, exposure_id)", async () => {
    const d = makeDeferred<any>("op-synth-shape");
    handleRemoteEvent({
      id: 7,
      kind: "peak_added",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: "op-synth-shape",
      payload: { q: 1.7, peak_curation_id: 314 },
      post_state: { analysis_inputs_hash: "new-hash", indices: [] },
    }, qc, qc.getMutationCache());
    const result = await d.promise;
    // Without these fields, peakAdd.onSuccess pushes a Peak with id=undefined
    // into the cache; TraceViewer then renders a permanent halo because
    // hoveredPeakId (undefined) === peak.id (undefined). Regression for
    // "circle around added peak doesn't clear without refresh."
    expect(result).toMatchObject({
      id: 314,
      exposure_id: 42,
      q: 1.7,
      source: "manual",
      excluded: false,
      analysis_inputs_hash: "new-hash",
    });
    expect(result.intensity).toBeNull();
    expect(result.prominence).toBeNull();
    expect(result.sharpness).toBeNull();
  });

  it("synthesizeResponseFromSse for peak_excluded maps auto_peak_id → id and derives source/excluded from kind", async () => {
    const d = makeDeferred<any>("op-pe-shape");
    handleRemoteEvent({
      id: 8,
      kind: "peak_excluded",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: "op-pe-shape",
      payload: { q: 0.5, auto_peak_id: 7 },
      post_state: { analysis_inputs_hash: "h-pe", indices: [] },
    }, qc, qc.getMutationCache());
    const result = await d.promise;
    // Without these fields, peakSetExcluded.onSuccess's `pk.id === peakOnly.id`
    // map matches no row (undefined === undefined fails on real Peak rows),
    // and the canonical state never replaces the optimistic one.
    expect(result).toMatchObject({
      id: 7,
      q: 0.5,
      source: "auto",
      excluded: true,
      analysis_inputs_hash: "h-pe",
    });
    // auto_peak_id MUST NOT leak into the response (it's not a Peak field).
    expect(result.auto_peak_id).toBeUndefined();
  });

  it("synthesizeResponseFromSse for peak_unexcluded derives excluded=false from kind", async () => {
    const d = makeDeferred<any>("op-pue-shape");
    handleRemoteEvent({
      id: 9,
      kind: "peak_unexcluded",
      entity_type: "exposure",
      entity_id: 42,
      client_op_id: "op-pue-shape",
      payload: { q: 0.5, auto_peak_id: 7 },
      post_state: { analysis_inputs_hash: "h-pue", indices: [] },
    }, qc, qc.getMutationCache());
    const result = await d.promise;
    expect(result).toMatchObject({
      id: 7,
      q: 0.5,
      source: "auto",
      excluded: false,
      analysis_inputs_hash: "h-pue",
    });
    expect(result.auto_peak_id).toBeUndefined();
  });

  it("synthesizeResponseFromSse for add_tag maps tag_id → id and adds source: 'manual'", async () => {
    const d = makeDeferred<any>("op-add-tag-shape");
    handleRemoteEvent({
      id: 50,
      kind: "add_tag",
      entity_type: "sample",
      entity_id: 10,
      client_op_id: "op-add-tag-shape",
      payload: { key: "category", value: "control", tag_id: 77, experiment_id: 1 },
    }, qc, qc.getMutationCache());
    const result = await d.promise;
    // Without this synthesis, addSampleTagMutator.onSuccess would read
    // response.id → undefined and response.source → undefined, landing a
    // malformed (undeletable, mis-classified) tag in the cache.
    expect(result).toMatchObject({
      id: 77,
      key: "category",
      value: "control",
      source: "manual",
    });
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
