/**
 * Cache-shape integrity for mutators that spread a route response into a
 * typed query cache.
 *
 * The bug class this catches: a mutator's onSuccess writes the full
 * response into a cache typed as `Foo`, but the response has extra
 * (event_id, view_row_id, sample_id, etc.) or missing fields. TypeScript
 * can't see this at runtime — issue #16 (group response pollution) and
 * sixth-pass issue #17 (Peak missing intensity/etc.) were both caught
 * only by re-reading the code or re-running the route.
 *
 * Each test runs a mutator end-to-end against a route-shaped mock and
 * asserts the cache row's keys match the type's exactly (both directions
 * strict — extras mean pollution, missing mean the route omitted a
 * type-required field).
 *
 * Discipline: the mocked response shape MUST be derived from the actual
 * route's emit, not from the TypeScript interface (since the type might
 * be wrong). Each mock fixture below is annotated with the file:line of
 * the route handler it mirrors.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { peakAddMutator } from "../../src/lib/queue/mutators/peakAdd";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import {
  peakExcludeMutator, peakUnexcludeMutator,
} from "../../src/lib/queue/mutators/peakSetExcluded";
import { createSpeculativeMutator } from "../../src/lib/queue/mutators/createSpeculative";
import {
  addAssignmentPhaseMutator,
  removeAssignmentPhaseMutator,
  setSampleCallStateMutator,
} from "../../src/lib/queue/mutators/assignment";
import { customIndexMutator } from "../../src/lib/queue/mutators/customIndex";
import {
  updateSampleMutator,
} from "../../src/lib/queue/mutators/trivial";
import { queryKeys } from "../../src/queries";
import { pendingDeferreds } from "../../src/lib/queue/deferred";

const PEAK_KEYS = new Set([
  "id", "exposure_id", "q", "intensity", "prominence", "sharpness",
  "source", "excluded",
]);
const ASSIGNMENT_KEYS = new Set([
  "exposure_id", "state", "members",
]);
const SAMPLE_KEYS = new Set([
  "id", "experiment_id", "name", "display_name", "notes",
]);
const INDEX_ENTRY_KEYS = new Set([
  "id", "exposure_id", "phase", "basis", "score", "r_squared",
  "lattice_d", "ngc", "status", "kind", "inputs_hash",
  "peaks", "predicted_q",
]);
function mockFetchOnce(body: unknown, status = 200): void {
  const original = globalThis.fetch;
  globalThis.fetch = (async () => {
    globalThis.fetch = original;
    return new Response(JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    });
  }) as typeof fetch;
}

async function runMutator<R>(
  qc: QueryClient,
  m: { onMutate: (p: any, q: QueryClient) => any;
       request: (p: any, s: AbortSignal) => Promise<R>;
       onSuccess: (p: any, r: R, q: QueryClient) => void; },
  flat: any,
): Promise<R> {
  const ctx = m.onMutate(flat, qc);
  try {
    const response = await m.request(flat, new AbortController().signal);
    m.onSuccess(flat, response, qc);
    return response;
  } catch (e) {
    ctx?.restore?.();
    throw e;
  }
}

function assertKeys(obj: unknown, expected: Set<string>, label: string): void {
  expect(obj).toBeTypeOf("object");
  expect(obj).not.toBeNull();
  const actual = new Set(Object.keys(obj as object));
  expect({ label, actual: [...actual].sort() })
    .toEqual({ label, actual: [...expected].sort() });
}

describe("Cache-shape integrity (mutator onSuccess writes type-shaped rows)", () => {
  let qc: QueryClient;
  beforeEach(() => {
    qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    pendingDeferreds.clear();
  });

  // -------------------------------------------------------------------------
  // Peak mutators
  // -------------------------------------------------------------------------

  it("peakAdd writes a Peak with exactly the 8 declared keys (issue #17 + deep-scan #1)", async () => {
    qc.setQueryData(queryKeys.peaks(5), []);
    // Mock derived from routes_peaks.jl:156-186 (POST response Dict).
    mockFetchOnce({
      id: 100, exposure_id: 5, q: 0.42,
      intensity: null, prominence: null, sharpness: null,
      source: "manual", excluded: false,
      event_id: 9, view_row_id: 100,
      analysis_inputs_hash: "h1",
    }, 201);
    await runMutator(qc, peakAddMutator, {
      kind: "peak_added",
      clientOpId: "op-shape-1",
      exposureId: 5, username: "alice", clientId: "tab-1",
      q: 0.42, payload: { q: 0.42 },
    });
    const list = qc.getQueryData<unknown[]>(queryKeys.peaks(5));
    expect(list).toHaveLength(1);
    assertKeys(list![0], PEAK_KEYS, "peakAdd cache row");
  });

  it("peakExclude writes a Peak with exactly 8 keys (no event metadata leak)", async () => {
    const initial = {
      id: 7, exposure_id: 5, q: 0.5, intensity: 1.2, prominence: 0.8,
      sharpness: 30.0, source: "auto", excluded: false,
    };
    qc.setQueryData(queryKeys.peaks(5), [initial]);
    // Mock derived from routes_peaks.jl PATCH /peaks/auto/:id response.
    mockFetchOnce({
      ...initial, excluded: true,
      event_id: 10, view_row_id: 11,
      analysis_inputs_hash: "h1",
    }, 200);
    await runMutator(qc, peakExcludeMutator, {
      kind: "peak_excluded",
      clientOpId: "op-shape-2",
      exposureId: 5, username: "alice", clientId: "tab-1",
      peakId: 7, q: 0.5,
      payload: { peakId: 7, q: 0.5 },
    });
    const list = qc.getQueryData<{ id: number }[]>(queryKeys.peaks(5));
    assertKeys(list!.find(p => p.id === 7)!, PEAK_KEYS,
      "peakExclude cache row");
  });

  it("peakUnexclude writes a Peak with exactly 8 keys", async () => {
    const initial = {
      id: 7, exposure_id: 5, q: 0.5, intensity: 1.2, prominence: 0.8,
      sharpness: 30.0, source: "auto", excluded: true,
    };
    qc.setQueryData(queryKeys.peaks(5), [initial]);
    mockFetchOnce({
      ...initial, excluded: false,
      event_id: 11, view_row_id: 12,
      analysis_inputs_hash: "h1",
    }, 200);
    await runMutator(qc, peakUnexcludeMutator, {
      kind: "peak_unexcluded",
      clientOpId: "op-shape-3",
      exposureId: 5, username: "alice", clientId: "tab-1",
      peakId: 7, q: 0.5,
      payload: { peakId: 7, q: 0.5 },
    });
    const list = qc.getQueryData<{ id: number }[]>(queryKeys.peaks(5));
    assertKeys(list!.find(p => p.id === 7)!, PEAK_KEYS,
      "peakUnexclude cache row");
  });

  it("peakRemove drops the manual peak from the peaks cache", async () => {
    // Only manual peaks are removable. Seed a manual peak so the optimistic
    // remove actually drops it.
    const initial = {
      id: 7, exposure_id: 5, q: 0.5, intensity: null, prominence: null,
      sharpness: null, source: "manual", excluded: false,
    };
    qc.setQueryData(queryKeys.peaks(5), [initial]);
    // Mock derived from routes_peaks.jl DELETE /peaks/manual/:id response.
    mockFetchOnce({
      event_id: 12, view_row_id: null,
      analysis_inputs_hash: "h2",
    }, 200);
    await runMutator(qc, peakRemoveMutator, {
      kind: "peak_removed",
      clientOpId: "op-shape-4",
      exposureId: 5, username: "alice", clientId: "tab-1",
      peakId: 7, payload: { peakId: 7 },
    });
    const peaks = qc.getQueryData<unknown[]>(queryKeys.peaks(5));
    expect(peaks).toHaveLength(0);
  });

  // -------------------------------------------------------------------------
  // Assignment mutators (Plan D-3) — onSuccess writes ONLY the assignment
  // cache (the 3-key Assignment shape); NEVER the exposure cache.
  // -------------------------------------------------------------------------

  it("addAssignmentPhase writes an Assignment with exactly 3 keys", async () => {
    qc.setQueryData(queryKeys.assignment(5), { exposure_id: 5, state: "indexed", members: [] });
    // Mock derived from routes_analysis.jl POST /assignment/members response.
    mockFetchOnce({
      exposure_id: 5, state: "indexed", members: [42],
      event_id: 21, view_row_id: 7,
    }, 200);
    await runMutator(qc, addAssignmentPhaseMutator, {
      kind: "assignment_add",
      clientOpId: "op-asg-1",
      exposureId: 5, username: "alice", clientId: "tab-1",
      indexId: 42, payload: { indexId: 42 },
    });
    const a = qc.getQueryData<unknown>(queryKeys.assignment(5));
    assertKeys(a, ASSIGNMENT_KEYS, "addAssignmentPhase cache row");
  });

  it("removeAssignmentPhase writes an Assignment with exactly 3 keys", async () => {
    qc.setQueryData(queryKeys.assignment(5), { exposure_id: 5, state: "indexed", members: [42] });
    mockFetchOnce({
      exposure_id: 5, state: "indexed", members: [],
      event_id: 22, view_row_id: 8,
    }, 200);
    await runMutator(qc, removeAssignmentPhaseMutator, {
      kind: "assignment_remove",
      clientOpId: "op-asg-2",
      exposureId: 5, username: "alice", clientId: "tab-1",
      indexId: 42, payload: { indexId: 42 },
    });
    assertKeys(qc.getQueryData(queryKeys.assignment(5)), ASSIGNMENT_KEYS, "removeAssignmentPhase cache row");
  });

  it("setSampleCallState writes an Assignment with exactly 3 keys on every frame of the sample", async () => {
    qc.setQueryData(queryKeys.exposures(1), [{ id: 5 }, { id: 6 }]);
    qc.setQueryData(queryKeys.assignment(5), { exposure_id: 5, state: "indexed", members: [42] });
    qc.setQueryData(queryKeys.assignment(6), { exposure_id: 6, state: "indexed", members: [43] });
    mockFetchOnce({
      sample_id: 1, state: "form_factor", phase: null, call_phases: [], members: [],
      event_id: 23, view_row_id: null,
    }, 200);
    await runMutator(qc, setSampleCallStateMutator, {
      kind: "sample_call_set_state",
      clientOpId: "op-asg-3",
      sampleId: 1, username: "alice", clientId: "tab-1",
      state: "form_factor", payload: { state: "form_factor" },
    });
    for (const id of [5, 6]) {
      const a = qc.getQueryData<{ state: string; members: number[] }>(queryKeys.assignment(id));
      assertKeys(a, ASSIGNMENT_KEYS, "setSampleCallState cache row");
      expect(a!.state).toBe("form_factor");
      expect(a!.members).toEqual([]);
    }
  });

  it("two pending assignment_add ops converge (reverse-rollback + insertion-replay)", () => {
    // Finding #5b: simulate two queued optimistic adds, then roll BOTH back in
    // reverse order and replay in insertion order — both members must end up
    // present (the rollback snapshots the whole Assignment, so it is symmetric).
    qc.setQueryData(queryKeys.assignment(5), { exposure_id: 5, state: "indexed", members: [] });
    const op1 = { kind: "assignment_add" as const, clientOpId: "c1", exposureId: 5,
      username: "a", clientId: "t", indexId: 10, payload: { indexId: 10 } };
    const op2 = { kind: "assignment_add" as const, clientOpId: "c2", exposureId: 5,
      username: "a", clientId: "t", indexId: 11, payload: { indexId: 11 } };
    const ctx1 = addAssignmentPhaseMutator.onMutate(op1, qc);
    const ctx2 = addAssignmentPhaseMutator.onMutate(op2, qc);
    // reverse-rollback
    ctx2.restore();
    ctx1.restore();
    expect(qc.getQueryData<{ members: number[] }>(queryKeys.assignment(5))!.members).toEqual([]);
    // insertion-replay
    addAssignmentPhaseMutator.onMutate(op1, qc);
    addAssignmentPhaseMutator.onMutate(op2, qc);
    expect(qc.getQueryData<{ members: number[] }>(queryKeys.assignment(5))!.members).toEqual([10, 11]);
  });

  it("customIndex appends the new IndexEntry and invalidates the assignment", async () => {
    qc.setQueryData(queryKeys.indices(5), []);
    qc.setQueryData(queryKeys.assignment(5), { exposure_id: 5, state: "indexed", members: [] });
    const inval = vi.spyOn(qc, "invalidateQueries");
    // Mock derived from routes_analysis.jl POST /custom-index response.
    // `peaks` is NON-EMPTY on purpose: insert_custom_index! claims the orders
    // that land, and the route reads them back. A fixture hardcoding `peaks: []`
    // would pin the pre-fix contract at this layer while layers 3/5 assert the
    // new one.
    mockFetchOnce({
      id: 77, exposure_id: 5, phase: "Pn3m", basis: 0.15, score: null, r_squared: null,
      lattice_d: 197, ngc: -1.5, status: "candidate", kind: "speculative", inputs_hash: "h",
      peaks: [{ peak_id: 3, ratio_position: 1, residual: 0.0, q_observed: 0.15 }],
      predicted_q: [0.15, 0.1837], bonnet: null,
      event_id: 30, view_row_id: 12,
    }, 200);
    await runMutator(qc, customIndexMutator, {
      kind: "custom_index_commit", clientOpId: "op-ci-1",
      exposureId: 5, username: "alice", clientId: "tab-1",
      phase: "Pn3m", basis: 0.15, ratios: [1, 1.2247448714],
      payload: { phase: "Pn3m", basis: 0.15, ratios: [1, 1.2247448714] },
    });
    type Claim = { peak_id: number; ratio_position: number; q_observed: number };
    const indices = qc.getQueryData<{ id: number; peaks: Claim[] }[]>(queryKeys.indices(5));
    const row = indices!.find((i) => i.id === 77)!;
    expect(row).toBeTruthy();
    // The field this PR exists to change must survive into the cache — an
    // IndexEntry whose `peaks` were dropped renders as claiming nothing.
    expect(row.peaks).toHaveLength(1);
    expect(row.peaks[0]!.ratio_position).toBe(1);
    expect(row.peaks[0]!.q_observed).toBe(0.15);
    expect(inval).toHaveBeenCalledWith({ queryKey: queryKeys.assignment(5) });
  });

  it("customIndex does NOT splice a phantom row on the SSE-wins synth (issue-#37)", () => {
    // On the SSE-wins own-tab race the deferred resolves off the FIRST frame
    // (speculative_created), whose mutator has no synthesizeFromSse — so the
    // response handed to customIndexMutator.onSuccess is the generic synth
    // {event_id, client_op_id, analysis_inputs_hash, index_id}, NOT a full
    // IndexEntry. The guard must invalidate instead of splicing a phantom
    // {id:undefined, phase:undefined} row.
    qc.setQueryData(queryKeys.indices(5), []);
    qc.setQueryData(queryKeys.assignment(5), { exposure_id: 5, state: "indexed", members: [] });
    const inval = vi.spyOn(qc, "invalidateQueries");
    const synth = {
      event_id: 30, client_op_id: "op-ci-1", analysis_inputs_hash: "h", index_id: 77,
    } as unknown as Parameters<typeof customIndexMutator.onSuccess>[1];
    const flat = {
      exposureId: 5, username: "alice", clientId: "tab-1",
      phase: "Pn3m", basis: 0.15,
    } as unknown as Parameters<typeof customIndexMutator.onSuccess>[0];
    customIndexMutator.onSuccess(flat, synth, qc);
    const indices = qc.getQueryData<{ id: number | undefined }[]>(queryKeys.indices(5));
    // No phantom row landed (cache stays empty; converges via the invalidate).
    expect(indices).toEqual([]);
    expect(indices!.some((i) => i.id === undefined)).toBe(false);
    expect(inval).toHaveBeenCalledWith({ queryKey: queryKeys.indices(5) });
    expect(inval).toHaveBeenCalledWith({ queryKey: queryKeys.assignment(5) });
  });

  it("createSpeculative writes an IndexEntry with exactly 13 keys", async () => {
    qc.setQueryData(queryKeys.indices(5), []);
    // Mock derived from routes_analysis.jl POST /speculative response (~line 374-397).
    mockFetchOnce({
      id: 99, exposure_id: 5, phase: "Pn3m", basis: 0.123,
      score: 0.85, r_squared: 0.99, lattice_d: 50.0, ngc: 0.5,
      status: "candidate", kind: "speculative", inputs_hash: "h-spec",
      peaks: [
        { peak_id: 7, ratio_position: 1, residual: 0.001, q_observed: 0.123 },
      ],
      predicted_q: [0.123, 0.174],
    }, 200);
    await runMutator(qc, createSpeculativeMutator, {
      kind: "speculative_created",
      clientOpId: "op-shape-8",
      exposureId: 5, username: "alice", clientId: "tab-1",
      phase: "Pn3m", anchor_peak_id: 7, anchor_ratio: 1, additional: [],
      payload: { phase: "Pn3m", anchor_peak_id: 7, anchor_ratio: 1, additional: [] },
    });
    const indices = qc.getQueryData<unknown[]>(queryKeys.indices(5));
    expect(indices).toHaveLength(1);
    assertKeys(indices![0], INDEX_ENTRY_KEYS, "createSpeculative cache row");
  });

  // -------------------------------------------------------------------------
  // Sample mutators
  // -------------------------------------------------------------------------

  it("updateSample merges only the patched fields (deep-scan Bug #2)", async () => {
    const initialSample = {
      id: 10, experiment_id: 1, display_name: "D1", name: "old", notes: "n",
    };
    qc.setQueryData(queryKeys.samples(1), [initialSample]);
    mockFetchOnce({
      id: 10, experiment_id: 1, display_name: "D1", name: "old", notes: "n",
      created_at: "2026-05-03",
    }, 200);
    await runMutator(qc, updateSampleMutator, {
      kind: "update_sample",
      clientOpId: "op-shape-9",
      sampleId: 10, experimentId: 1, username: "alice", clientId: "tab-1",
      display_name: "new",
      payload: { sampleId: 10 },
    });
    const row = qc.getQueryData<{ notes: string }[]>(queryKeys.samples(1))![0];
    assertKeys(row!, SAMPLE_KEYS, "updateSample list cache row");
    expect(row!.notes).toBe("n");
  });

  // -------------------------------------------------------------------------
  // Peak-id collision (manual usage test surfaced this — auto_peaks.id and
  // peak_curations.id share a namespace on the wire; cache filters/maps must
  // disambiguate by source.)
  // -------------------------------------------------------------------------

  it("peakRemove with id colliding between auto and manual only drops the manual peak", async () => {
    const colliding = [
      { id: 3, exposure_id: 5, q: 0.075, intensity: 1.2, prominence: 0.8,
        sharpness: 30, source: "auto", excluded: false },
      { id: 3, exposure_id: 5, q: 0.180, intensity: null, prominence: null,
        sharpness: null, source: "manual", excluded: false },
    ];
    qc.setQueryData(queryKeys.peaks(5), colliding);
    mockFetchOnce({ event_id: 9, view_row_id: null, analysis_inputs_hash: "h2" }, 200);
    await runMutator(qc, peakRemoveMutator, {
      kind: "peak_removed",
      clientOpId: "op-collide-rm",
      exposureId: 5, username: "alice", clientId: "tab-1",
      peakId: 3, payload: { peakId: 3 },
    });
    const after = qc.getQueryData<{ id: number; source: string }[]>(queryKeys.peaks(5))!;
    expect(after).toHaveLength(1);
    expect(after[0]!.source).toBe("auto");
    expect(after[0]!.id).toBe(3);
  });

  it("peakExclude with id colliding only flips the auto peak", async () => {
    const colliding = [
      { id: 3, exposure_id: 5, q: 0.075, intensity: 1.2, prominence: 0.8,
        sharpness: 30, source: "auto", excluded: false },
      { id: 3, exposure_id: 5, q: 0.180, intensity: null, prominence: null,
        sharpness: null, source: "manual", excluded: false },
    ];
    qc.setQueryData(queryKeys.peaks(5), colliding);
    mockFetchOnce({
      ...colliding[0], excluded: true,
      event_id: 10, view_row_id: 11, analysis_inputs_hash: "h2",
    }, 200);
    await runMutator(qc, peakExcludeMutator, {
      kind: "peak_excluded",
      clientOpId: "op-collide-ex",
      exposureId: 5, username: "alice", clientId: "tab-1",
      peakId: 3, q: 0.075, payload: { peakId: 3, q: 0.075 },
    });
    const after = qc.getQueryData<{ id: number; source: string; excluded: boolean }[]>(queryKeys.peaks(5))!;
    expect(after.find(p => p.source === "auto")!.excluded).toBe(true);
    expect(after.find(p => p.source === "manual")!.excluded).toBe(false);
  });

});
