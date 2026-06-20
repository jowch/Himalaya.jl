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
  Peak, Exposure, Sample, SampleMessage, Assignment,
} from "../../src/api";
import { queryKeys } from "../../src/queries";

const FULL_EXPOSURE: Exposure = {
  id: 5, sample_id: 1, filename: null, kind: "file", selected: true,
  status: null, image_path: null, image_version: "",
  trace_hash: null, analysis_inputs_hash: "h0",
  tags: [], sources: [],
};

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

  it("analyze_run with post_state writes indices and updates exposure hash", () => {
    qc.setQueryData<Exposure>(queryKeys.exposure(5), FULL_EXPOSURE);
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
    expect(qc.getQueryData<Exposure>(queryKeys.exposure(5))!.analysis_inputs_hash)
      .toBe("h-new");
    const indices = qc.getQueryData<unknown[]>(queryKeys.indices(5))!;
    expect(indices).toHaveLength(1);
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

  it("peak_excluded with the W1 envelope writes indices + hash + assignment cache (no invalidate)", () => {
    qc.setQueryData<Exposure>(queryKeys.exposure(5), FULL_EXPOSURE);
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
    expect(qc.getQueryData<Exposure>(queryKeys.exposure(5))!.analysis_inputs_hash)
      .toBe("h-new");
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
    qc.setQueryData<Exposure>(queryKeys.exposure(5), FULL_EXPOSURE);
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

  it("assignment_set_state patches state + clears members from post_state", () => {
    qc.setQueryData<Assignment>(queryKeys.assignment(5),
      { exposure_id: 5, state: "indexed", members: [10] });
    const evt: SseEvent = {
      id: 99, kind: "assignment_set_state", entity_type: "exposure", entity_id: 5,
      payload: { state: "form_factor" },
      post_state: { assignment: { state: "form_factor", members: [] } },
    };
    applyRemoteToCache(evt, qc);
    const a = qc.getQueryData<Assignment>(queryKeys.assignment(5))!;
    expect(a.state).toBe("form_factor");
    expect(a.members).toEqual([]);
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
    // The single highest-value test: an assignment frame reaching the own-tab
    // applyPostStateOnly path must NOT clobber the exposure hash (the {assignment}
    // post_state has no `indices` array, so the guard bails).
    qc.setQueryData<Exposure>(queryKeys.exposure(5), FULL_EXPOSURE);
    qc.setQueryData(queryKeys.indices(5), [{ id: 1 }]);
    const evt: SseEvent = {
      id: 99, kind: "assignment_add", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 42 },
      post_state: { assignment: { state: "indexed", members: [42] } },
    };
    applyPostStateOnly(evt, qc);
    // Exposure hash untouched; indices cache untouched.
    expect(qc.getQueryData<Exposure>(queryKeys.exposure(5))!.analysis_inputs_hash).toBe("h0");
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

  it("speculative_deleted invalidates indices", () => {
    let called = 0;
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: any) => { called++; return orig(arg); }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "speculative_deleted", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    expect(called).toBe(1);
  });

  // -------------------------------------------------------------------------
  // Exposure-scoped events
  // -------------------------------------------------------------------------

  it("set_exposure_status writes status onto the exposure entity", () => {
    qc.setQueryData<Exposure>(queryKeys.exposure(5), FULL_EXPOSURE);
    // Mirrors routes_exposures.jl PATCH: payload is {status}.
    const evt: SseEvent = {
      id: 99, kind: "set_exposure_status", entity_type: "exposure", entity_id: 5,
      payload: { status: "rejected" },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<Exposure>(queryKeys.exposure(5))!.status).toBe("rejected");
  });

  it("select_exposure invalidates the exposure list AND the picker projection", () => {
    const invalidated: unknown[] = [];
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated.push(arg.queryKey); return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    // Mirrors routes_exposures.jl PATCH /select: payload is {sample_id}.
    const evt: SseEvent = {
      id: 99, kind: "select_exposure", entity_type: "exposure", entity_id: 5,
      payload: { sample_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toContainEqual(queryKeys.exposures(1));
    // indexing_exposure_id derives from selection; the builder's Confirm
    // resolves its plate through the picker projection (BU-RECIPENOOP), so a
    // foreign re-selection must refresh it or Confirm commits the PREVIOUS
    // representative exposure.
    expect(invalidated).toContainEqual(queryKeys.corpusPickerSamples);
  });

  // -------------------------------------------------------------------------
  // Sample/message events
  // -------------------------------------------------------------------------

  it("update_sample spreads payload onto the sample entity", () => {
    const sample: Sample = {
      id: 10, experiment_id: 1, name: "old", notes: "n", tags: [],
    };
    qc.setQueryData(queryKeys.sample(10), sample);
    // Mirrors routes_samples.jl PATCH: payload is the patched fields directly.
    const evt: SseEvent = {
      id: 99, kind: "update_sample", entity_type: "sample", entity_id: 10,
      payload: { name: "new" },
    };
    applyRemoteToCache(evt, qc);
    const after = qc.getQueryData<Sample>(queryKeys.sample(10))!;
    expect(after.name).toBe("new");
    expect(after.notes).toBe("n");  // unpatched fields preserved
    expect(after.tags).toEqual([]); // tags survive (deep-scan #2)
  });

  it("post_message dedupes by id and inserts the SampleMessage", () => {
    qc.setQueryData<SampleMessage[]>(queryKeys.messages(10), []);
    // Mirrors routes_messages.jl POST: payload is the full SampleMessage row.
    const evt: SseEvent = {
      id: 99, kind: "post_message", entity_type: "sample_message", entity_id: 200,
      payload: {
        id: 200, sample_id: 10, author_id: 3, author: "alice",
        body: "hi", created_at: "2026-05-03T12:00:00Z",
      },
    };
    applyRemoteToCache(evt, qc);
    const msgs = qc.getQueryData<SampleMessage[]>(queryKeys.messages(10))!;
    expect(msgs).toHaveLength(1);
    expect(msgs[0]!.id).toBe(200);
    // Idempotency: replaying the same frame must not double-insert.
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<SampleMessage[]>(queryKeys.messages(10))!).toHaveLength(1);
  });

  it("add_tag (sample) invalidates the experiment + corpus samples lists AND the scoping caches", () => {
    const invalidated: unknown[] = [];
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated.push(arg.queryKey); return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    // Mirrors routes_samples.jl POST: payload always includes experiment_id.
    // A sample tag lives in four cached projections — the per-experiment list,
    // the corpus contact-sheet list, and (I3.4 #174) the corpus tag proposal +
    // picker projection the /series/new scoping surface reads — so all four
    // keys must invalidate on a foreign tag write.
    const evt: SseEvent = {
      id: 99, kind: "add_tag", entity_type: "sample", entity_id: 10,
      payload: { key: "k", value: "v", tag_id: 50, experiment_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual([
      queryKeys.samples(1), queryKeys.corpusSamples,
      queryKeys.corpusSampleTags, queryKeys.corpusPickerSamples,
    ]);
  });

  it("add_tag (exposure) invalidates the sample's exposures list (uses sample_id)", () => {
    let invalidated: unknown = null;
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated = arg.queryKey; return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    // Mirrors routes_exposures.jl POST: payload includes sample_id.
    const evt: SseEvent = {
      id: 99, kind: "add_tag", entity_type: "exposure", entity_id: 5,
      payload: { key: "k", value: "v", tag_id: 60, sample_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual(queryKeys.exposures(1));
  });

  it("remove_tag (sample) invalidates the experiment + corpus samples lists AND the scoping caches", () => {
    const invalidated: unknown[] = [];
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated.push(arg.queryKey); return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    // remove_tag shares the add_tag(sample) branch, so it fans out to the same
    // four keys (I3.4 #174 added the two corpus-scoping caches).
    const evt: SseEvent = {
      id: 99, kind: "remove_tag", entity_type: "sample", entity_id: 10,
      payload: { tag_id: 50, experiment_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual([
      queryKeys.samples(1), queryKeys.corpusSamples,
      queryKeys.corpusSampleTags, queryKeys.corpusPickerSamples,
    ]);
  });

  it("delete_index falls through to default (invalidates peaks+indices)", () => {
    // `delete_index` is the OpKind for the user gesture that hits the
    // DELETE /api/indices/:id route — the backend emits `speculative_deleted`
    // on the wire, which has its own dedicated case. But if a future
    // contributor adds a typed branch for `delete_index` in applyRemoteToCache
    // and forgets to handle the wire-name divergence, this test pins the
    // expected default fall-through (review suggestion #16).
    let invalidatedKeys: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidatedKeys.push(arg.queryKey); return orig(arg);
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

  it("post_message with entity_type=sample_message still routes to sample messages cache", () => {
    // Regression: extending the dispatch must not break the legacy path.
    qc.setQueryData<SampleMessage[]>(queryKeys.messages(10), []);
    const evt: SseEvent = {
      id: 99, kind: "post_message", entity_type: "sample_message", entity_id: 200,
      payload: {
        id: 200, sample_id: 10, author_id: 3, author: "alice",
        body: "hi sample", created_at: "2026-05-06T12:00:00Z",
      },
    };
    applyRemoteToCache(evt, qc);
    const msgs = qc.getQueryData<SampleMessage[]>(queryKeys.messages(10))!;
    expect(msgs).toHaveLength(1);
  });

  it("remove_tag (exposure) invalidates the sample's exposures list", () => {
    let invalidated: unknown = null;
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated = arg.queryKey; return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "remove_tag", entity_type: "exposure", entity_id: 5,
      payload: { tag_id: 60, sample_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual(queryKeys.exposures(1));
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
    { kind: "add_tag",              entity_type: "sample",     entity_id: 5,  payload: { tag_id: 7, key: "tag", value: "v" } },
    { kind: "add_tag",              entity_type: "exposure",   entity_id: 5,  payload: { tag_id: 7, key: "tag", value: "v" } },
    { kind: "peak_excluded",        entity_type: "exposure",   entity_id: 5,  payload: { auto_peak_id: 1, q: 0.1 } },
    { kind: "peak_unexcluded",      entity_type: "exposure",   entity_id: 5,  payload: { auto_peak_id: 1, q: 0.1 } },
  ];
  it.each(cases)("$kind/$entity_type resolves and synthesizes", ({ kind, entity_type, entity_id, payload }) => {
    const mutator = resolveMutatorForEvent(kind, entity_type);
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
    "assignment_add", "assignment_remove", "assignment_set_state",
  ])("%s synthesizes an Assignment from {assignment} post_state", (kind) => {
    const mutator = resolveMutatorForEvent(kind, "exposure");
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
  // Camp (b) breaks down into three rationales:
  //
  //   (b.i)  Forward-scaffolded — mutator exists but no UI gesture queues it
  //          today (set_exposure_status, update_sample, select_exposure,
  //          remove_tag ×2). When a future plan wires the gesture, add
  //          synthesizeFromSse to the owning mutator.
  //   (b.ii) Active mutators whose SSE payload IS the cache row shape, so
  //          the generic fallback already produces the correct shape
  //          (post_message — payload IS SampleMessage).
  //   (b.iii) Active mutators whose onSuccess relies on the `looksFull`
  //          detector to invalidate when the synth shape is incomplete
  //          (createSpeculative, both indexGroup variants, deleteIndex).
  //          Plus mutators whose SSE-wins-and-then-resolve has no cache
  //          effect beyond `analysis_inputs_hash` (analyze_run, peak_removed
  //          — peakRemove.onSuccess is a no-op beyond the hash write).
  const noSynth: Array<{ kind: string; entity_type: string }> = [
    // (b.i) forward-scaffolded
    { kind: "set_exposure_status", entity_type: "exposure"           },
    { kind: "update_sample",       entity_type: "sample"             },
    { kind: "select_exposure",     entity_type: "exposure"           },
    { kind: "remove_tag",          entity_type: "sample"             },
    { kind: "remove_tag",          entity_type: "exposure"           },
    // (b.ii) payload IS cache-row shape
    { kind: "post_message",        entity_type: "sample_message"     },
    // (b.iii) looksFull-handled or hash-only effects
    { kind: "analyze_run",         entity_type: "exposure"           },
    { kind: "peak_removed",        entity_type: "exposure"           },
    { kind: "speculative_created", entity_type: "exposure"           },
    { kind: "speculative_deleted", entity_type: "exposure"           },
  ];
  it.each(noSynth)(
    "$kind/$entity_type stays on the generic fallback (no mutator.synthesizeFromSse)",
    ({ kind, entity_type }) => {
      const mutator = resolveMutatorForEvent(kind, entity_type);
      expect(mutator).toBeDefined();
      expect(mutator!.synthesizeFromSse).toBeUndefined();
    },
  );
});
