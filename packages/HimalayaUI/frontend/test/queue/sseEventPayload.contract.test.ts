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
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import type { SseEvent } from "../../src/lib/queue/types";
import type {
  Peak, Exposure, Sample, GroupEntry, SampleMessage,
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
  // Group/index events
  // -------------------------------------------------------------------------

  it("index_confirmed appends index_id to the matching group's members", () => {
    qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "custom", active: true, members: [10] },
    ]);
    // Mirrors routes_analysis.jl POST /groups/:id/members: payload is
    // {group_id, index_id}.
    const evt: SseEvent = {
      id: 99, kind: "index_confirmed", entity_type: "exposure", entity_id: 5,
      payload: { group_id: 1, index_id: 42 },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<GroupEntry[]>(queryKeys.groups(5))![0]!.members)
      .toEqual([10, 42]);
  });

  it("index_unconfirmed removes index_id from the matching group's members", () => {
    qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "custom", active: true, members: [10, 42] },
    ]);
    const evt: SseEvent = {
      id: 99, kind: "index_unconfirmed", entity_type: "exposure", entity_id: 5,
      payload: { group_id: 1, index_id: 42 },
    };
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<GroupEntry[]>(queryKeys.groups(5))![0]!.members)
      .toEqual([10]);
  });

  it("speculative_created invalidates indices+groups (no inline cache write)", () => {
    let called = 0;
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: any) => { called++; return orig(arg); }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "speculative_created", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    expect(called).toBe(2);  // indices + groups
  });

  it("speculative_deleted invalidates indices+groups", () => {
    let called = 0;
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: any) => { called++; return orig(arg); }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "speculative_deleted", entity_type: "exposure", entity_id: 5,
      payload: { index_id: 7 },
    };
    applyRemoteToCache(evt, qc);
    expect(called).toBe(2);
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

  it("select_exposure invalidates the parent sample's exposure list", () => {
    let invalidated: unknown = null;
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated = arg.queryKey; return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    // Mirrors routes_exposures.jl PATCH /select: payload is {sample_id}.
    const evt: SseEvent = {
      id: 99, kind: "select_exposure", entity_type: "exposure", entity_id: 5,
      payload: { sample_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual(queryKeys.exposures(1));
  });

  // -------------------------------------------------------------------------
  // Sample/message events
  // -------------------------------------------------------------------------

  it("update_sample spreads payload onto the sample entity", () => {
    const sample: Sample = {
      id: 10, experiment_id: 1, label: "D1", name: "old", notes: "n", tags: [],
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

  it("add_tag (sample) invalidates the experiment's samples list (uses experiment_id)", () => {
    let invalidated: unknown = null;
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated = arg.queryKey; return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    // Mirrors routes_samples.jl POST: payload includes experiment_id (the
    // parent invalidation key).
    const evt: SseEvent = {
      id: 99, kind: "add_tag", entity_type: "sample", entity_id: 10,
      payload: { key: "k", value: "v", tag_id: 50, experiment_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual(queryKeys.samples(1));
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

  it("remove_tag (sample) invalidates the experiment's samples list", () => {
    let invalidated: unknown = null;
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated = arg.queryKey; return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "remove_tag", entity_type: "sample", entity_id: 10,
      payload: { tag_id: 50, experiment_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual(queryKeys.samples(1));
  });

  it("delete_index falls through to default (invalidates peaks+indices+groups)", () => {
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
    // Default branch invalidates peaks, indices, groups for the entity.
    expect(invalidatedKeys).toEqual([
      queryKeys.peaks(5),
      queryKeys.indices(5),
      queryKeys.groups(5),
    ]);
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
