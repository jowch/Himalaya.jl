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
import { resolveMutatorForEvent } from "../../src/lib/queue/mutatorRegistry";
import type { SseEvent } from "../../src/lib/queue/types";
import { remoteForeignEvent } from "./helpers";
import type {
  Peak, Exposure, Sample, GroupEntry, SampleMessage,
  ComparisonMessage, ComparisonSummary,
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

  it("index_confirmed for an unknown group_id invalidates the groups list (issue #37 Bug 1c)", () => {
    // Foreign tab confirmed an index for the FIRST time on this exposure,
    // creating a fresh custom group on the backend. Other tabs only have the
    // auto group cached; the surgical update would silently miss the new
    // custom group, leaving the foreign confirmation invisible until refetch.
    qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "auto", active: true, members: [] },
    ]);
    let invalidated = false;
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = (filters: any) => {
      const k = filters?.queryKey ?? [];
      if (Array.isArray(k) && k[0] === "exposure" && k[1] === 5 && k[2] === "groups") {
        invalidated = true;
      }
      return orig(filters);
    };
    const evt: SseEvent = {
      id: 99, kind: "index_confirmed", entity_type: "exposure", entity_id: 5,
      payload: { group_id: 5, index_id: 42 },  // group_id=5 NOT in cache
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toBe(true);
  });

  it("index_unconfirmed for an unknown group_id invalidates the groups list", () => {
    qc.setQueryData<GroupEntry[]>(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "auto", active: true, members: [] },
    ]);
    let invalidated = false;
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = (filters: any) => {
      const k = filters?.queryKey ?? [];
      if (Array.isArray(k) && k[0] === "exposure" && k[1] === 5 && k[2] === "groups") {
        invalidated = true;
      }
      return orig(filters);
    };
    const evt: SseEvent = {
      id: 99, kind: "index_unconfirmed", entity_type: "exposure", entity_id: 5,
      payload: { group_id: 5, index_id: 42 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toBe(true);
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
      id: 10, experiment_id: 1, display_name: "D1", name: "old", notes: "n", tags: [],
    };
    qc.setQueryData(queryKeys.sample(10), sample);
    // Mirrors routes_samples.jl PATCH: payload is the patched fields directly.
    const evt: SseEvent = {
      id: 99, kind: "update_sample", entity_type: "sample", entity_id: 10,
      payload: { display_name: "new" },
    };
    applyRemoteToCache(evt, qc);
    const after = qc.getQueryData<Sample>(queryKeys.sample(10))!;
    expect(after.display_name).toBe("new");
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

  it("add_tag (sample) invalidates BOTH the experiment samples list and the corpus list", () => {
    const invalidated: unknown[] = [];
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated.push(arg.queryKey); return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    // Mirrors routes_samples.jl POST: payload always includes experiment_id.
    // A sample tag lives in two cached projections — the per-experiment list
    // and the corpus contact-sheet list — so both keys must invalidate.
    const evt: SseEvent = {
      id: 99, kind: "add_tag", entity_type: "sample", entity_id: 10,
      payload: { key: "k", value: "v", tag_id: 50, experiment_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual([queryKeys.samples(1), queryKeys.corpusSamples]);
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

  it("remove_tag (sample) invalidates BOTH the experiment samples list and the corpus list", () => {
    const invalidated: unknown[] = [];
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated.push(arg.queryKey); return Promise.resolve();
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "remove_tag", entity_type: "sample", entity_id: 10,
      payload: { tag_id: 50, experiment_id: 1 },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidated).toEqual([queryKeys.samples(1), queryKeys.corpusSamples]);
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

  // -------------------------------------------------------------------------
  // Compare page events (Phase 3). Three kinds:
  // - comparison_created  → invalidate comparison + members + listings
  // - comparison_submitted → invalidate same set (rebuilds full state)
  // - comparison_deleted  → REMOVE (not invalidate) entity + filter listings
  // - post_message (entity_type=comparison_message) → comparison thread cache
  // -------------------------------------------------------------------------

  it("comparison_created invalidates comparison(id), comparisonMembers(id), and listings", () => {
    const invalidatedKeys: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidatedKeys.push(arg.queryKey); return orig(arg);
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "comparison_created", entity_type: "comparison", entity_id: 42,
      payload: { title: "X", members: [] },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidatedKeys).toEqual([
      queryKeys.comparison(42),
      queryKeys.comparisonMembers(42),
      ["comparisons"],  // prefix invalidation hits all listing scopes
    ]);
  });

  it("comparison_submitted invalidates comparison(id), comparisonMembers(id), and listings", () => {
    const invalidatedKeys: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidatedKeys.push(arg.queryKey); return orig(arg);
    }) as typeof qc.invalidateQueries;
    const evt: SseEvent = {
      id: 99, kind: "comparison_submitted", entity_type: "comparison", entity_id: 42,
      payload: { title: "X edited", members: [] },
    };
    applyRemoteToCache(evt, qc);
    expect(invalidatedKeys).toEqual([
      queryKeys.comparison(42),
      queryKeys.comparisonMembers(42),
      ["comparisons"],
    ]);
  });

  it("comparison_deleted removes entity caches (NOT invalidates) and filters listings", () => {
    qc.setQueryData(queryKeys.comparison(42), { id: 42 });
    qc.setQueryData(queryKeys.comparisonMembers(42), []);
    qc.setQueryData(queryKeys.comparisonMessages(42), []);
    qc.setQueryData(queryKeys.comparisonForks(42), []);
    qc.setQueryData<ComparisonSummary[]>(queryKeys.comparisons("all"), [
      { id: 42, title: "doomed", description: null, content_hash: "h",
        created_by: 1, created_at: null, updated_at: null,
        forked_from_id: null, forked_at_hash: null,
        view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
        last_event_at: null, author_username: null, member_count: 0,
        member_phases: [], member_phase_count: 0, has_stale_members: false },
      { id: 99, title: "kept", description: null, content_hash: "h",
        created_by: 1, created_at: null, updated_at: null,
        forked_from_id: null, forked_at_hash: null,
        view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
        last_event_at: null, author_username: null, member_count: 0,
        member_phases: [], member_phase_count: 0, has_stale_members: false },
    ]);
    let invalidated = false;
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidated = true; return orig(arg);
    }) as typeof qc.invalidateQueries;

    const evt: SseEvent = {
      id: 99, kind: "comparison_deleted", entity_type: "comparison", entity_id: 42,
      payload: { id: 42 },
    };
    applyRemoteToCache(evt, qc);

    // Entity caches REMOVED, not invalidated
    expect(invalidated).toBe(false);
    expect(qc.getQueryState(queryKeys.comparison(42))).toBeUndefined();
    expect(qc.getQueryState(queryKeys.comparisonMembers(42))).toBeUndefined();
    expect(qc.getQueryState(queryKeys.comparisonMessages(42))).toBeUndefined();
    expect(qc.getQueryState(queryKeys.comparisonForks(42))).toBeUndefined();
    // Listing filtered (id pruned, others retained)
    const listing = qc.getQueryData<ComparisonSummary[]>(queryKeys.comparisons("all"))!;
    expect(listing.map((c) => c.id)).toEqual([99]);
  });

  it("post_message with entity_type=comparison_message routes to comparisonMessages cache", () => {
    qc.setQueryData<ComparisonMessage[]>(queryKeys.comparisonMessages(42), []);
    const evt: SseEvent = {
      id: 99, kind: "post_message", entity_type: "comparison_message", entity_id: 200,
      payload: {
        id: 200, comparison_id: 42, author_id: 3, author: "alice",
        body: "hi cmp", created_at: "2026-05-06T12:00:00Z",
      },
    };
    applyRemoteToCache(evt, qc);
    const msgs = qc.getQueryData<ComparisonMessage[]>(queryKeys.comparisonMessages(42))!;
    expect(msgs).toHaveLength(1);
    expect(msgs[0]!.id).toBe(200);
    expect(msgs[0]!.body).toBe("hi cmp");
    // Idempotency: replay must not double-insert
    applyRemoteToCache(evt, qc);
    expect(qc.getQueryData<ComparisonMessage[]>(queryKeys.comparisonMessages(42))!)
      .toHaveLength(1);
    // And the sample messages cache is untouched (no cross-thread pollution)
    expect(qc.getQueryData(queryKeys.messages(42))).toBeUndefined();
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
    { kind: "comparison_created",   entity_type: "comparison", entity_id: 11, payload: { title: "T" } },
    { kind: "comparison_submitted", entity_type: "comparison", entity_id: 11, payload: { title: "T" } },
    { kind: "comparison_deleted",   entity_type: "comparison", entity_id: 11, payload: {} },
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
  //          (post_message ×2 — payload IS SampleMessage / ComparisonMessage).
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
    { kind: "post_message",        entity_type: "comparison_message" },
    // (b.iii) looksFull-handled or hash-only effects
    { kind: "analyze_run",         entity_type: "exposure"           },
    { kind: "peak_removed",        entity_type: "exposure"           },
    { kind: "index_confirmed",     entity_type: "exposure"           },
    { kind: "index_unconfirmed",   entity_type: "exposure"           },
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
