/**
 * F-WIPE W2 — peak_* frames consume the W1 assignment post_state envelope.
 *
 * The W1 backend contract: ALL reanalyzing SSE frames (peak_added /
 * peak_removed / peak_excluded / peak_unexcluded AND analyze_run — the
 * manual-reanalyze route was the last producer to gain the envelope) carry
 *   post_state = { analysis_inputs_hash, indices,
 *                  assignment: { state, members },        // same serializer
 *                  assignment_dropped?: string[] }        // per-member phases
 *
 * Frontend obligations pinned here:
 *   (a) assignment envelope present  → assignment cache WRITTEN (no invalidate)
 *   (b) envelope absent (old backend) → assignment cache UNTOUCHED
 *   (c) assignment_dropped non-empty → exactly ONE warning toast (aggregated)
 *   (d) dropped absent / empty       → no toast
 *   (e) exactly-once across the three handleRemoteEvent paths — the
 *       announcement lives in applyPostStateOnly, the single per-frame
 *       chokepoint (Case 1 deferred-match, self-echo guard, foreign
 *       applyRemoteToCache→applyPostState are mutually exclusive per frame),
 *       and no mutator onSuccess re-processes the frame.
 */
import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache, applyPostStateOnly } from "../../src/lib/queue/applyRemoteToCache";
import { handleRemoteEvent } from "../../src/lib/queue/replayCoordinator";
import { makeDeferred, pendingDeferreds } from "../../src/lib/queue/deferred";
import { peakAddMutator } from "../../src/lib/queue/mutators/peakAdd";
import type { SseEvent } from "../../src/lib/queue/types";
import type { Assignment, Exposure, PeakAddResponse } from "../../src/api";
import { queryKeys } from "../../src/queries";
import { getClientId } from "../../src/lib/clientId";
import { setToastImpl } from "../../src/lib/toast";

const EXPOSURE: Exposure = {
  id: 5, sample_id: 1, filename: null, kind: "file", selected: true,
  status: null, image_path: null, image_version: "",
  trace_hash: null, analysis_inputs_hash: "h0",
  tags: [], sources: [],
};

const SEEDED_ASSIGNMENT: Assignment = { exposure_id: 5, state: "indexed", members: [10, 42] };

/** A peak_excluded frame carrying the full W1 curation post_state. */
function w1Frame(overrides: Partial<SseEvent> = {}): SseEvent {
  return {
    id: 99, kind: "peak_excluded", entity_type: "exposure", entity_id: 5,
    payload: { q: 0.5, auto_peak_id: 7 },
    post_state: {
      analysis_inputs_hash: "h-new",
      indices: [{ id: 1, exposure_id: 5, phase: "Im3m", basis: 0.1 }],
      assignment: { state: "indexed", members: [10] },
    },
    ...overrides,
  };
}

describe("peak_* frames with the W1 assignment envelope", () => {
  let qc: QueryClient;
  let toastCalls: Array<{ msg: string; kind: string }>;

  beforeEach(() => {
    pendingDeferreds.clear();
    qc = new QueryClient();
    qc.setQueryData<Exposure>(queryKeys.exposure(5), EXPOSURE);
    qc.setQueryData(queryKeys.indices(5), []);
    qc.setQueryData<Assignment>(queryKeys.assignment(5), SEEDED_ASSIGNMENT);
    toastCalls = [];
    setToastImpl((msg, kind) => { toastCalls.push({ msg, kind }); });
  });
  afterEach(() => { setToastImpl(null); });

  // ---- (a) envelope present → assignment cache written, no invalidation ----

  it("foreign peak frame writes the assignment cache from post_state.assignment", () => {
    const invalidateSpy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache(w1Frame(), qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(
      { exposure_id: 5, state: "indexed", members: [10] });
    // Splice, not invalidate — replay-without-refetch.
    expect(invalidateSpy).not.toHaveBeenCalledWith(
      { queryKey: queryKeys.assignment(5) });
    // Curation post_state still lands as before.
    expect(qc.getQueryData<Exposure>(queryKeys.exposure(5))!.analysis_inputs_hash)
      .toBe("h-new");
    expect(qc.getQueryData<unknown[]>(queryKeys.indices(5))).toHaveLength(1);
  });

  it("own-op path (applyPostStateOnly) writes the assignment cache too", () => {
    applyPostStateOnly(w1Frame(), qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(
      { exposure_id: 5, state: "indexed", members: [10] });
  });

  // ---- (b) envelope absent (old backend) → assignment cache untouched ------

  it("peak frame WITHOUT the assignment envelope leaves the assignment cache untouched", () => {
    const invalidateSpy = vi.spyOn(qc, "invalidateQueries");
    const frame = w1Frame({
      post_state: { analysis_inputs_hash: "h-new", indices: [] },
    });
    applyRemoteToCache(frame, qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(SEEDED_ASSIGNMENT);
    expect(invalidateSpy).not.toHaveBeenCalledWith(
      { queryKey: queryKeys.assignment(5) });
  });

  it("applyPostStateOnly without the envelope does not touch the assignment cache", () => {
    applyPostStateOnly(w1Frame({
      post_state: { analysis_inputs_hash: "h-new", indices: [] },
    }), qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(SEEDED_ASSIGNMENT);
  });

  // ---- (c)/(d) assignment_dropped toast gating ------------------------------

  it("assignment_dropped non-empty fires exactly one warning toast with aggregated copy", () => {
    applyRemoteToCache(w1Frame({
      post_state: {
        analysis_inputs_hash: "h-new", indices: [],
        assignment: { state: "null", members: [] },
        assignment_dropped: ["Pn3m", "Lamellar"],
      },
    }), qc);
    expect(toastCalls).toEqual([{
      msg: "Pn3m and Lamellar indices dropped from the call. They no longer fit the peaks.",
      kind: "warning",
    }]);
  });

  it("no assignment_dropped → no toast", () => {
    applyRemoteToCache(w1Frame(), qc);
    expect(toastCalls).toEqual([]);
  });

  it("empty assignment_dropped → no toast", () => {
    applyRemoteToCache(w1Frame({
      post_state: {
        analysis_inputs_hash: "h-new", indices: [],
        assignment: { state: "indexed", members: [10] },
        assignment_dropped: [],
      },
    }), qc);
    expect(toastCalls).toEqual([]);
  });

  // ---- (e) exactly-once across the three handleRemoteEvent paths -----------

  const droppedFrame = (overrides: Partial<SseEvent> = {}): SseEvent => w1Frame({
    post_state: {
      analysis_inputs_hash: "h-new", indices: [],
      assignment: { state: "null", members: [] },
      assignment_dropped: ["Pn3m"],
    },
    ...overrides,
  });

  it("own-op Case 1 (deferred match): one toast, and the resolved mutator onSuccess does not re-fire", async () => {
    const d = makeDeferred<PeakAddResponse>("op-1");
    handleRemoteEvent(
      droppedFrame({
        kind: "peak_added",
        client_op_id: "op-1",
        payload: { q: 0.42, peak_curation_id: 100 },
      }),
      qc, qc.getMutationCache());
    expect(toastCalls).toHaveLength(1);
    // The deferred resolution feeds the mutator's onSuccess (the HTTP request
    // is aborted) — run it exactly as useQueueMutation would and pin that no
    // second toast appears. This is the own-op double-fire guard: the HTTP
    // confirmation path never re-processes post_state.
    const response = await d.promise;
    peakAddMutator.onSuccess(
      { kind: "peak_added", payload: { q: 0.42 }, clientOpId: "op-1",
        exposureId: 5, q: 0.42, username: undefined, clientId: getClientId() },
      response, qc);
    expect(toastCalls).toHaveLength(1);
  });

  it("own-op self-echo guard (HTTP won, deferred already cleared): exactly one toast", () => {
    handleRemoteEvent(
      droppedFrame({ client_id: getClientId(), client_op_id: "already-cleared" }),
      qc, qc.getMutationCache());
    expect(toastCalls).toHaveLength(1);
  });

  it("foreign frame (replay-as-rerun path): exactly one toast", () => {
    handleRemoteEvent(
      droppedFrame({ client_id: "some-other-tab", client_op_id: "not-ours" }),
      qc, qc.getMutationCache());
    expect(toastCalls).toHaveLength(1);
  });

  // ---- analyze_run frames (last producer to gain the W1 envelope) ----------
  // Pin tests: the analyze_run branch already funnels through the same
  // applyPostState() → applyPostStateOnly chokepoint as the peak_* branches,
  // so the envelope is consumed identically. These rows make that a contract
  // rather than a coincidence.

  it("analyze_run frame with assignment + assignment_dropped writes the assignment cache and fires exactly one toast", () => {
    const invalidateSpy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache({
      id: 99, kind: "analyze_run", entity_type: "exposure", entity_id: 5,
      payload: {},
      post_state: {
        analysis_inputs_hash: "h-new", indices: [],
        assignment: { state: "null", members: [] },
        assignment_dropped: ["Pn3m"],
      },
    }, qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(
      { exposure_id: 5, state: "null", members: [] });
    expect(invalidateSpy).not.toHaveBeenCalledWith(
      { queryKey: queryKeys.assignment(5) });
    expect(toastCalls).toEqual([{
      msg: "Pn3m index dropped from the call. Its index no longer fits the peaks.",
      kind: "warning",
    }]);
    expect(qc.getQueryData<Exposure>(queryKeys.exposure(5))!.analysis_inputs_hash)
      .toBe("h-new");
  });

  it("analyze_run frame without the envelope leaves the assignment cache untouched (no toast)", () => {
    const invalidateSpy = vi.spyOn(qc, "invalidateQueries");
    applyRemoteToCache({
      id: 99, kind: "analyze_run", entity_type: "exposure", entity_id: 5,
      payload: {},
      post_state: { analysis_inputs_hash: "h-new", indices: [] },
    }, qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(SEEDED_ASSIGNMENT);
    expect(invalidateSpy).not.toHaveBeenCalledWith(
      { queryKey: queryKeys.assignment(5) });
    expect(toastCalls).toEqual([]);
  });
});
