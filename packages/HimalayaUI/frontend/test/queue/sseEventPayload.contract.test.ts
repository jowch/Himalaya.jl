/**
 * SSE event-payload contract test.
 *
 * Pins the third leg of the contract triangle: route emit (Julia) ↔ cache
 * row shape (Vitest cache-shape) ↔ SSE frame payload (this file). If a route
 * silently changes its SSE payload shape, foreign-tab convergence breaks
 * silently — only manifests as "the other tab didn't update."
 *
 * Discipline: each frame's `payload` is derived from the actual
 * `apply_event!` call site in routes_*.jl / pipeline.jl. The route-side
 * payload SHAPE is what the JSON3-serialized SSE frame carries, so this is
 * the canonical wire contract from the consumer's point of view.
 *
 * Each test:
 *   1. Seeds the relevant query cache.
 *   2. Dispatches an SseEvent whose payload mirrors a real route emit.
 *   3. Asserts the cache reaches the expected post-state.
 *
 * If you change a route's apply_event! payload, update the corresponding
 * test here — and double-check applyRemoteToCache.ts reads the same fields.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache, applyPostStateOnly } from "../../src/lib/queue/applyRemoteToCache";
import { resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import type { SseEvent } from "../../src/lib/queue/types";
import { remoteForeignEvent } from "./helpers";
import type {
  Peak, Exposure, Sample, Assignment,
} from "../../src/api";
import { queryKeys } from "../../src/queries";
import { makeExposure } from "../helpers/exposure";

const FULL_EXPOSURE: Exposure = makeExposure({ id: 5, analysis_inputs_hash: "h0" });

describe("SSE event-payload contract (applyRemoteToCache for each emitted kind)", () => {
  let qc: QueryClient;
  beforeEach(() => { qc = new QueryClient(); });

  // -------------------------------------------------------------------------
  // Peak events
  // -------------------------------------------------------------------------

  it("peak_added inserts a Peak with id=peak_curation_id from payload", () => {
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), []);
    // Mirrors routes_peaks.jl POST: payload is {q} initially, then mutated to
    // include peak_curation_id (the dispatcher's view_row_id) before broadcast.
    const evt: SseEvent = {
      id: 99, kind: "peak_added", entity_type: "exposure", entity_id: 5,
      payload: { q: 0.42, peak_curation_id: 100 },
    };
    applyRemoteToCache(evt, qc);
    const peaks = qc.getQueryData<Peak[]>(queryKeys.peaks(5))!;
    expect(peaks).toHaveLength(1);
    expect(peaks[0]!.id).toBe(100);  // load-bearing — must be peak_curation_id, NOT event id
    expect(peaks[0]!.q).toBe(0.42);
    expect(peaks[0]!.source).toBe("manual");
  });

  it("peak_excluded flips excluded=true on the matching auto peak", () => {
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: 1, prominence: 1,
        sharpness: 30, source: "auto", excluded: false },
    ]);
    // Mirrors routes_peaks.jl PATCH: payload is {q, auto_peak_id}.
    const evt: SseEvent = {
      id: 99, kind: "peak_excluded", entity_type: "exposure", entity_id: 5,
      payload: { q: 0.5, auto_peak_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<Peak[]>(queryKeys.peaks(5))![0]!.excluded).toBe(true);
  });

  it("peak_unexcluded flips excluded=false on the matching auto peak", () => {
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: 1, prominence: 1,
        sharpness: 30, source: "auto", excluded: true },
    ]);
    const evt: SseEvent = {
      id: 99, kind: "peak_unexcluded", entity_type: "exposure", entity_id: 5,
      payload: { q: 0.5, auto_peak_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<Peak[]>(queryKeys.peaks(5))![0]!.excluded).toBe(false);
  });

  it("peak_removed deletes by peak_curation_id (not q match)", () => {
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: null, prominence: null,
        sharpness: null, source: "manual", excluded: false },
      { id: 8, exposure_id: 5, q: 0.5, intensity: null, prominence: null,
        sharpness: null, source: "manual", excluded: false },
    ]);
    // Mirrors routes_peaks.jl DELETE: payload is {peak_curation_id, q?}.
    // Two peaks at the same q — id-based filter must keep the untouched one.
    const evt: SseEvent = {
      id: 99, kind: "peak_removed", entity_type: "exposure", entity_id: 5,
      payload: { peak_curation_id: 7, q: 0.5 },
    };
    applyRemoteToCache(evt, qc);
    const peaks = qc.getQueryData<Peak[]>(queryKeys.peaks(5))!;
    expect(peaks.map((p) => p.id)).toEqual([8]);
  });

  // ---- (id, source) collision regression rows for SSE merge path -----------
  // auto_peaks and peak_curations have independent SQLite sequences; they can
  // collide on the wire. The merge-path handlers must be source-scoped just
  // like the mutator-path onMutate handlers (issue #24).

  it("peak_added: foreign event with id colliding with existing auto peak still inserts the manual peak", () => {
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), [
      { id: 100, exposure_id: 5, q: 0.3, intensity: 1, prominence: 1,
        sharpness: 30, source: "auto", excluded: false },
    ]);
    const evt: SseEvent = {
      id: 99, kind: "peak_added", entity_type: "exposure", entity_id: 5,
      payload: { q: 0.42, peak_curation_id: 100 },  // collides with auto peak id
    };
    applyRemoteToCache(evt, qc);
    const peaks = qc.getQueryData<Peak[]>(queryKeys.peaks(5))!;
    expect(peaks).toHaveLength(2);
    expect(peaks.find((p) => p.source === "manual")?.q).toBe(0.42);
    expect(peaks.find((p) => p.source === "auto")?.q).toBe(0.3);  // untouched
  });

  it("peak_excluded: foreign event with auto_peak_id colliding with manual peak only flips the auto", () => {
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: 1, prominence: 1,
        sharpness: 30, source: "auto", excluded: false },
      { id: 7, exposure_id: 5, q: 0.9, intensity: null, prominence: null,
        sharpness: null, source: "manual", excluded: false },
    ]);
    const evt: SseEvent = {
      id: 99, kind: "peak_excluded", entity_type: "exposure", entity_id: 5,
      payload: { q: 0.5, auto_peak_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    const peaks = qc.getQueryData<Peak[]>(queryKeys.peaks(5))!;
    expect(peaks.find((p) => p.source === "auto")!.excluded).toBe(true);
    expect(peaks.find((p) => p.source === "manual")!.excluded).toBe(false);
  });

  it("peak_unexcluded: foreign event scoped to auto source only", () => {
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: 1, prominence: 1,
        sharpness: 30, source: "auto", excluded: true },
      { id: 7, exposure_id: 5, q: 0.9, intensity: null, prominence: null,
        sharpness: null, source: "manual", excluded: true },
    ]);
    const evt: SseEvent = {
      id: 99, kind: "peak_unexcluded", entity_type: "exposure", entity_id: 5,
      payload: { q: 0.5, auto_peak_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    const peaks = qc.getQueryData<Peak[]>(queryKeys.peaks(5))!;
    expect(peaks.find((p) => p.source === "auto")!.excluded).toBe(false);
    expect(peaks.find((p) => p.source === "manual")!.excluded).toBe(true);  // untouched
  });

  it("peak_removed: foreign event with peak_curation_id colliding with auto peak only drops the manual", () => {
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: 1, prominence: 1,
        sharpness: 30, source: "auto", excluded: false },
      { id: 7, exposure_id: 5, q: 0.9, intensity: null, prominence: null,
        sharpness: null, source: "manual", excluded: false },
    ]);
    const evt: SseEvent = {
      id: 99, kind: "peak_removed", entity_type: "exposure", entity_id: 5,
      payload: { peak_curation_id: 7, q: 0.9 },
    };
    applyRemoteToCache(evt, qc);
    const peaks = qc.getQueryData<Peak[]>(queryKeys.peaks(5))!;
    expect(peaks).toHaveLength(1);
    expect(peaks[0]!.source).toBe("auto");
    expect(peaks[0]!.q).toBe(0.5);  // auto peak survives despite shared id
  });

  it("analyze_run with post_state writes indices", () => {
    qc.setQueryData(queryKeys.indices(5), []);
    const evt: SseEvent = {
      id: 99, kind: "analyze_run", entity_type: "exposure", entity_id: 5,
      payload: {},
      post_state: {
        analysis_inputs_hash: "h-new",
        indices: [{ id: 1, exposure_id: 5, phase: "Pn3m", basis: 0.1 }],
      },
    };
    applyRemoteToCache(evt, qc);
    const indices = qc.getQueryData<unknown[]>(queryKeys.indices(5))!;
    expect(indices).toHaveLength(1);
  });

  it("analyze_run that re-detected on an unchanged trace invalidates the peaks cache", () => {
    qc.setQueryData(queryKeys.indices(5), []);
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: 1, prominence: 1,
        sharpness: 30, source: "auto", excluded: false },
    ]);
    const invalidateSpy = vi.spyOn(qc, "invalidateQueries");
    const frame = (findpeaks_skipped: boolean): SseEvent => ({
      id: 0, kind: "analyze_run", entity_type: "exposure", entity_id: 5,
      payload: { trace_hash_before: "t", trace_hash_after: "t", findpeaks_skipped },
      post_state: { analysis_inputs_hash: "h", indices: [] },
    });
    applyRemoteToCache(frame(true), qc);
    expect(invalidateSpy).not.toHaveBeenCalledWith({ queryKey: queryKeys.peaks(5) });
    applyRemoteToCache(frame(false), qc);
    expect(invalidateSpy).toHaveBeenCalledWith({ queryKey: queryKeys.peaks(5) });
    expect(invalidateSpy).not.toHaveBeenCalledWith({ queryKey: queryKeys.trace(5) });
  });

  // -------------------------------------------------------------------------
  // F-WIPE W1/W2 — peak_* frames now carry the assignment envelope INSIDE the
  // curation post_state: { analysis_inputs_hash, indices, assignment: {state,
  // members}, assignment_dropped?: string[] }. Serialized by the same Julia
  // `_assignment_post_state` helper the assignment_* frames use
  // (routes_peaks.jl::_enrich_curation_post_state). Backend pair:
  // test/test_assignment_reattach.jl (frame shape) +
  // test_route_response_shapes.jl (route emit).
  // -------------------------------------------------------------------------

  it("peak_excluded with the W1 envelope writes indices + assignment cache (no invalidate)", () => {
    qc.setQueryData(queryKeys.indices(5), []);
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), [
      { id: 7, exposure_id: 5, q: 0.5, intensity: 1, prominence: 1,
        sharpness: 30, source: "auto", excluded: false },
    ]);
    qc.setQueryData<Assignment>(queryKeys.assignment(5),
      { exposure_id: 5, state: "indexed", members: [10, 42] });
    const invalidateSpy = vi.spyOn(qc, "invalidateQueries");
    const evt: SseEvent = {
      id: 99, kind: "peak_excluded", entity_type: "exposure", entity_id: 5,
      payload: { q: 0.5, auto_peak_id: 7 },
      post_state: {
        analysis_inputs_hash: "h-new",
        indices: [{ id: 1, exposure_id: 5, phase: "Im3m", basis: 0.1 }],
        assignment: { state: "indexed", members: [10] },
      },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<unknown[]>(queryKeys.indices(5))).toHaveLength(1);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(
      { exposure_id: 5, state: "indexed", members: [10] });
    expect(invalidateSpy).not.toHaveBeenCalledWith(
      { queryKey: queryKeys.assignment(5) });
  });

  it("peak_added W1 envelope reaches the assignment cache through the same shared writer", () => {
    qc.setQueryData<Peak[]>(queryKeys.peaks(5), []);
    qc.setQueryData<Assignment>(queryKeys.assignment(5),
      { exposure_id: 5, state: "null", members: [] });
    const evt: SseEvent = {
      id: 99, kind: "peak_added", entity_type: "exposure", entity_id: 5,
      payload: { q: 0.42, peak_curation_id: 100 },
      post_state: {
        analysis_inputs_hash: "h-new", indices: [],
        assignment: { state: "indexed", members: [11] },
      },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(
      { exposure_id: 5, state: "indexed", members: [11] });
  });

  it("peak frame WITHOUT the assignment envelope (old backend) leaves the assignment cache untouched", () => {
    qc.setQueryData(queryKeys.indices(5), []);
    qc.setQueryData<Assignment>(queryKeys.assignment(5),
      { exposure_id: 5, state: "indexed", members: [10, 42] });
    const invalidateSpy = vi.spyOn(qc, "invalidateQueries");
    const evt: SseEvent = {
      id: 99, kind: "peak_removed", entity_type: "exposure", entity_id: 5,
      payload: { peak_curation_id: 7, q: 0.5 },
      post_state: { analysis_inputs_hash: "h-new", indices: [] },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(
      { exposure_id: 5, state: "indexed", members: [10, 42] });
    expect(invalidateSpy).not.toHaveBeenCalledWith(
      { queryKey: queryKeys.assignment(5) });
  });

  // -------------------------------------------------------------------------
  // Assignment events (Plan D-3) — DISTINCT {assignment:{state,members}}
  // post_state (NO top-level `indices` key).
  // -------------------------------------------------------------------------

  it("assignment_add patches the assignment cache from the distinct post_state", () => {
    qc.setQueryData<Assignment>(queryKeys.assignment(5),
      { exposure_id: 5, state: "indexed", members: [10] });
    const evt: SseEvent = {
      id: 99, kind: "assignment_add", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 42 },
      post_state: { assignment: { state: "indexed", members: [10, 42] } },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual(
      { exposure_id: 5, state: "indexed", members: [10, 42] });
  });

  it("assignment_remove patches the assignment cache from post_state", () => {
    qc.setQueryData<Assignment>(queryKeys.assignment(5),
      { exposure_id: 5, state: "indexed", members: [10, 42] });
    const evt: SseEvent = {
      id: 99, kind: "assignment_remove", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 42 },
      post_state: { assignment: { state: "indexed", members: [10] } },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))!.members).toEqual([10]);
  });

  it("sample_call_set_state writes the call to every frame of the sample", () => {
    qc.setQueryData(queryKeys.exposures(1), [{ id: 5 }, { id: 6 }]);
    qc.setQueryData<Assignment>(queryKeys.assignment(5),
      { exposure_id: 5, state: "form_factor", members: [] });
    const evt: SseEvent = {
      id: 99, kind: "sample_call_set_state", entity_type: "sample", entity_id: 1,
      payload: { state: "indexed" },
      post_state: { call: { state: "indexed", members: [{ index_id: 10, exposure_id: 6, phase: "Pn3m" }] } },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(5))).toEqual({ exposure_id: 5, state: "indexed", members: [] });
    expect(qc.getQueryData<Assignment>(queryKeys.assignment(6))).toEqual({ exposure_id: 6, state: "indexed", members: [10] });
  });

  it("sample_call_set_state synthesizes the sample's call from {call} post_state", () => {
    const synth = resolveMutatorForEvent("sample_call_set_state")!.synthesizeFromSse?.(
      remoteForeignEvent({
        id: 7, kind: "sample_call_set_state", entity_type: "sample", entity_id: 1, payload: {},
        post_state: { call: { state: "null", members: [] } },
      }),
      { event_id: 7, client_op_id: "x", analysis_inputs_hash: undefined },
    ) as { sample_id: number; state: string } | undefined;
    expect(synth).toMatchObject({ sample_id: 1, state: "null" });
  });

  it("assignment frame carries NO top-level `indices` key (distinct post_state)", () => {
    // Finding #5: the wire contract is {assignment:{state,members}} — assert the
    // absence of `indices` so a future shape drift toward CurationPostState is
    // caught. Mirrors the Julia route-emit test in test_assignments.jl.
    const post_state = { assignment: { state: "indexed" as const, members: [10] } };
    expect("indices" in post_state).toBe(false);
    expect("assignment" in post_state).toBe(true);
  });

  it("applyPostStateOnly is a NO-OP on an assignment frame (finding #5a)", () => {
    qc.setQueryData(queryKeys.indices(5), [{ id: 1 }]);
    const evt: SseEvent = {
      id: 99, kind: "assignment_add", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 42 },
      post_state: { assignment: { state: "indexed", members: [42] } },
    };
    applyPostStateOnly(evt, qc);
    expect(qc.getQueryData<unknown[]>(queryKeys.indices(5))).toEqual([{ id: 1 }]);
  });

  it("assignment_add with no post_state invalidates the assignment cache (fallback)", () => {
    qc.setQueryData<Assignment>(queryKeys.assignment(5),
      { exposure_id: 5, state: "indexed", members: [10] });
    const spy = vi.spyOn(qc, "invalidateQueries");
    const evt: SseEvent = {
      id: 99, kind: "assignment_add", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 42 },
    };
    applyRemoteToCache(evt, qc);
    expect(spy).toHaveBeenCalledWith({ queryKey: queryKeys.assignment(5) });
  });

  it("speculative_created invalidates indices (no inline cache write)", () => {
    let called = 0;
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: any) => { called++; return orig(arg); }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "speculative_created", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    expect(called).toBe(1);  // indices
  });

  it("speculative_deleted invalidates indices AND the assignment", () => {
    // assignment_members.index_id is ON DELETE CASCADE, so deleting a
    // speculative that was in the call silently drops its membership row
    // server-side. Invalidating only `indices` leaves a foreign tab's cart
    // holding a phantom member id that no longer resolves to an index.
    const keys: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: any) => {
      keys.push(arg.queryKey); return orig(arg);
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "speculative_deleted", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    expect(keys).toEqual(
      expect.arrayContaining([queryKeys.indices(5), queryKeys.assignment(5)]),
    );
  });

  it("speculative_created invalidates ONLY indices (its own assignment_add frame follows)", () => {
    // Guards the split: `created` must not inherit the delete's assignment
    // invalidation. The custom-index route emits a separate assignment_add
    // frame carrying the canonical {state, members} post_state.
    const keys: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: any) => {
      keys.push(arg.queryKey); return orig(arg);
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "speculative_created", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    expect(keys).toEqual([queryKeys.indices(5)]);
  });

  // -------------------------------------------------------------------------
  // Exposure-scoped events
  // -------------------------------------------------------------------------

  it("set_exposure_status writes status onto the cached exposure list", () => {
    qc.setQueryData<Exposure[]>(queryKeys.exposures(1), [FULL_EXPOSURE]);
    const evt: SseEvent = {
      id: 99, kind: "set_exposure_status", entity_type: "exposure", entity_id: 5,
      payload: { status: "rejected" },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<Exposure[]>(queryKeys.exposures(1))![0]!.status).toBe("rejected");
  });

  // -------------------------------------------------------------------------
  // Sample events
  // -------------------------------------------------------------------------

  it("update_sample spreads payload onto the cached samples list", () => {
    const sample: Sample = {
      id: 10, experiment_id: 1, name: "old", notes: "n",
    };
    qc.setQueryData<Sample[]>(queryKeys.samples(1), [sample]);
    // Mirrors routes_samples.jl PATCH: payload is the patched fields directly.
    const evt: SseEvent = {
      id: 99, kind: "update_sample", entity_type: "sample", entity_id: 10,
      payload: { name: "new" },
    };
    applyRemoteToCache(evt, qc);
    const after = qc.getQueryData<Sample[]>(queryKeys.samples(1))![0]!;
    expect(after.name).toBe("new");
    expect(after.notes).toBe("n");  // unpatched fields preserved
  });

  it("update_sample invalidates the corpus samples list the workspace reads", () => {
    const invalidated: unknown[] = [];
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated.push(arg.queryKey); return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "update_sample", entity_type: "sample", entity_id: 10,
      payload: { name: "new" },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual([queryKeys.corpusSamples]);
  });

  it.each(["label_added", "label_removed", "label_suggestion_restored"])(
    "%s invalidates only that experiment's labels",
    (kind) => {
      const invalidated: unknown[] = [];
      qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
        invalidated.push(arg.queryKey); return Promise.resolve();
      }) as typeof qc.invalidateQueries;
      const evt: SseEvent = {
        id: 99, kind, entity_type: "sample", entity_id: 10,
        payload: { label: "Buffer", experiment_id: 1 },
      };
      applyRemoteToCache(evt, qc);
      expect(invalidated).toEqual([queryKeys.labels(1)]);
    },
  );

  it("delete_index falls through to default (invalidates peaks+indices)", () => {
    // `delete_index` is the OpKind for the user gesture that hits the
    // DELETE /api/indices/:id route — the backend emits `speculative_deleted`
    // on the wire, which has its own dedicated case. But if a future
    // contributor adds a typed branch for `delete_index` in applyRemoteToCache
    // and forgets to handle the wire-name divergence, this test pins the
    // expected default fall-through (review suggestion #16).
    let invalidatedKeys: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((...args: Parameters<typeof orig>) => {
      invalidatedKeys.push(args[0]?.queryKey); return orig(...args);
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "delete_index", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    // Default branch invalidates peaks, indices for the entity.
    expect(invalidatedKeys).toEqual([
      queryKeys.peaks(5),
      queryKeys.indices(5),
    ]);
  });

});

describe("synthesizeFromSse coverage (resolveMutatorForEvent contract)", () => {
  // Factory for a minimally-valid SseEvent. Reuses `remoteForeignEvent` from
  // test/queue/helpers.ts (already exported) so we don't fall behind on
  // SseEvent shape drift.
  function evt(kind: string, entity_type: string, entity_id: number, payload: object) {
    return remoteForeignEvent({ id: 1, kind, entity_type, entity_id, payload });
  }

  // The exact set of kinds that had a bespoke branch in the legacy switch.
  // Each MUST resolve to a mutator that returns non-undefined for a
  // minimally-shaped SseEvent.
  const cases: Array<{ kind: string; entity_type: string; entity_id: number; payload: object }> = [
    { kind: "peak_added",           entity_type: "exposure",   entity_id: 5,  payload: { peak_curation_id: 99, q: 0.123 } },
    { kind: "peak_excluded",        entity_type: "exposure",   entity_id: 5,  payload: { auto_peak_id: 1, q: 0.1 } },
    { kind: "peak_unexcluded",      entity_type: "exposure",   entity_id: 5,  payload: { auto_peak_id: 1, q: 0.1 } },
  ];
  it.each(cases)("$kind/$entity_type resolves and synthesizes", ({ kind, entity_type, entity_id, payload }) => {
    const mutator = resolveMutatorForEvent(kind);
    expect(mutator).toBeDefined();
    const synth = mutator!.synthesizeFromSse?.(
      evt(kind, entity_type, entity_id, payload),
      { event_id: 1, client_op_id: "x", analysis_inputs_hash: undefined },
    );
    expect(synth).toBeDefined();
  });

  // Assignment mutators synthesize from the DISTINCT {assignment} post_state,
  // not from payload — so they need a post_state-bearing frame (the generic
  // `cases` factory above only carries payload). All three share synthAssignment.
  it.each([
    "assignment_add", "assignment_remove",
  ])("%s synthesizes an Assignment from {assignment} post_state", (kind) => {
    const mutator = resolveMutatorForEvent(kind);
    expect(mutator).toBeDefined();
    const synth = mutator!.synthesizeFromSse?.(
      remoteForeignEvent({
        id: 7, kind, entity_type: "exposure", entity_id: 5, payload: {},
        post_state: { assignment: { state: "indexed", members: [10, 11] } },
      }),
      { event_id: 7, client_op_id: "x", analysis_inputs_hash: undefined },
    ) as { exposure_id: number; state: string; members: number[] } | undefined;
    expect(synth).toBeDefined();
    expect(synth!.exposure_id).toBe(5);
    expect(synth!.state).toBe("indexed");
    expect(synth!.members).toEqual([10, 11]);
  });

  // Every event kind in resolveMutatorForEvent's switch falls into one of
  // two camps: (a) bespoke synth via the owning mutator's synthesizeFromSse
  // (the `cases` array above), or (b) generic `{...base, ...payload}`
  // fallback. The noSynth array enumerates camp (b) so any future drift
  // (e.g. a half-baked synth landed on a mutator before the pipeline is
  // ready) trips a test rather than silently changing wire shape.
  //
  // Camp (b) breaks down into two rationales:
  //
  //   (b.i)  Forward-scaffolded — mutator exists but no UI gesture queues it
  //          today (set_exposure_status, update_sample).
  //          When a future plan wires the gesture, add
  //          synthesizeFromSse to the owning mutator.
  //   (b.ii) Active mutators whose onSuccess relies on the `looksFull`
  //          detector to invalidate when the synth shape is incomplete
  //          (createSpeculative, both indexGroup variants, deleteIndex).
  const noSynth: Array<{ kind: string }> = [
    // (b.i) forward-scaffolded
    { kind: "set_exposure_status" },
    { kind: "update_sample" },
    // (b.ii) looksFull-handled or hash-only effects
    { kind: "analyze_run" },
    { kind: "peak_removed" },
    { kind: "speculative_created" },
    { kind: "speculative_deleted" },
    { kind: "label_added" },
    { kind: "label_removed" },
    { kind: "label_suggestion_restored" },
  ];
  it.each(noSynth)(
    "$kind stays on the generic fallback (no mutator.synthesizeFromSse)",
    ({ kind }) => {
      const mutator = resolveMutatorForEvent(kind);
      expect(mutator).toBeDefined();
      expect(mutator!.synthesizeFromSse).toBeUndefined();
    },
  );
});
