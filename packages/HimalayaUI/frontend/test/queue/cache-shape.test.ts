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
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { peakAddMutator } from "../../src/lib/queue/mutators/peakAdd";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import {
  peakExcludeMutator, peakUnexcludeMutator,
} from "../../src/lib/queue/mutators/peakSetExcluded";
import {
  addIndexToGroupMutator, removeIndexFromGroupMutator,
} from "../../src/lib/queue/mutators/indexGroup";
import { createSpeculativeMutator } from "../../src/lib/queue/mutators/createSpeculative";
import { reanalyzeExposureMutator } from "../../src/lib/queue/mutators/reanalyzeExposure";
import {
  updateSampleMutator,
  addSampleTagMutator,
  addExposureTagMutator,
  postSampleMessageMutator,
} from "../../src/lib/queue/mutators/trivial";
import { saveComparisonMutator } from "../../src/lib/queue/mutators/saveComparison";
import { deleteComparisonMutator } from "../../src/lib/queue/mutators/deleteComparison";
import type { ComparisonSummary } from "../../src/api";
import { queryKeys } from "../../src/queries";
import { pendingDeferreds } from "../../src/lib/queue/deferred";

const PEAK_KEYS = new Set([
  "id", "exposure_id", "q", "intensity", "prominence", "sharpness",
  "source", "excluded",
]);
const GROUP_KEYS = new Set([
  "id", "exposure_id", "kind", "active", "members",
]);
const SAMPLE_KEYS = new Set([
  "id", "experiment_id", "name", "display_name", "notes", "tags",
]);
const EXPOSURE_KEYS = new Set([
  "id", "sample_id", "filename", "kind", "selected", "status",
  "image_path", "image_version", "tags", "sources", "trace_hash",
  "analysis_inputs_hash",
]);
const SAMPLE_TAG_KEYS = new Set(["id", "key", "value", "source"]);
const EXPOSURE_TAG_KEYS = new Set(["id", "key", "value", "source"]);
const SAMPLE_MESSAGE_KEYS = new Set([
  "id", "sample_id", "author_id", "author", "body", "created_at",
]);
const INDEX_ENTRY_KEYS = new Set([
  "id", "exposure_id", "phase", "basis", "score", "r_squared",
  "lattice_d", "ngc", "status", "kind", "inputs_hash",
  "peaks", "predicted_q",
]);
const COMPARISON_KEYS = new Set([
  "id", "title", "description", "content_hash",
  "created_by", "created_at", "updated_at",
  "forked_from_id", "forked_at_hash", "forked_from_title",
  "members",
]);
const COMPARISON_MEMBER_KEYS = new Set([
  "id", "comparison_id", "exposure_id", "display_order",
  "band_height", "y_offset", "normalization",
  "color_override", "label_override",
  "q_window_min", "q_window_max",
  "peak_display", "snapshot", "is_stale",
  "created_by", "created_at",
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

const FULL_EXPOSURE = {
  id: 5, sample_id: 1, filename: null, kind: "file" as const, selected: true,
  status: null, image_path: null, trace_hash: null,
  analysis_inputs_hash: "h0", tags: [], sources: [], image_version: "",
};

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
    qc.setQueryData(queryKeys.exposure(5), FULL_EXPOSURE);
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
    // Exposure cache hash was rewritten — make sure no extra fields leaked.
    assertKeys(qc.getQueryData(queryKeys.exposure(5)), EXPOSURE_KEYS,
      "peakAdd exposure cache");
  });

  it("peakExclude writes a Peak with exactly 8 keys (no event metadata leak)", async () => {
    const initial = {
      id: 7, exposure_id: 5, q: 0.5, intensity: 1.2, prominence: 0.8,
      sharpness: 30.0, source: "auto", excluded: false,
    };
    qc.setQueryData(queryKeys.peaks(5), [initial]);
    qc.setQueryData(queryKeys.exposure(5), FULL_EXPOSURE);
    // Mock derived from routes_peaks.jl PATCH /peaks/:id response.
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
    assertKeys(qc.getQueryData(queryKeys.exposure(5)), EXPOSURE_KEYS,
      "peakExclude exposure cache");
  });

  it("peakUnexclude writes a Peak with exactly 8 keys", async () => {
    const initial = {
      id: 7, exposure_id: 5, q: 0.5, intensity: 1.2, prominence: 0.8,
      sharpness: 30.0, source: "auto", excluded: true,
    };
    qc.setQueryData(queryKeys.peaks(5), [initial]);
    qc.setQueryData(queryKeys.exposure(5), FULL_EXPOSURE);
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

  it("peakRemove updates exposure cache shape with exactly 12 keys", async () => {
    // Only manual peaks are removable. Seed a manual peak so the optimistic
    // remove actually drops it.
    const initial = {
      id: 7, exposure_id: 5, q: 0.5, intensity: null, prominence: null,
      sharpness: null, source: "manual", excluded: false,
    };
    qc.setQueryData(queryKeys.peaks(5), [initial]);
    qc.setQueryData(queryKeys.exposure(5), FULL_EXPOSURE);
    // Mock derived from routes_peaks.jl DELETE /peaks/:id response.
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
    // The peaks list lost the row — but the survivors must still be Peak-shaped
    // and the exposure cache hash-rewrite must not have polluted the entity.
    const peaks = qc.getQueryData<unknown[]>(queryKeys.peaks(5));
    expect(peaks).toHaveLength(0);
    assertKeys(qc.getQueryData(queryKeys.exposure(5)), EXPOSURE_KEYS,
      "peakRemove exposure cache");
  });

  it("reanalyzeExposure updates exposure cache shape with exactly 12 keys", async () => {
    qc.setQueryData(queryKeys.exposure(5), FULL_EXPOSURE);
    // Mock derived from routes_exposures.jl POST /exposures/:id/analyze response.
    mockFetchOnce({
      id: 5, analyzed: true, analysis_inputs_hash: "h3",
    }, 200);
    await runMutator(qc, reanalyzeExposureMutator, {
      kind: "reanalyze_exposure",
      clientOpId: "op-shape-5",
      exposureId: 5, username: "alice", clientId: "tab-1",
      payload: {},
    });
    assertKeys(qc.getQueryData(queryKeys.exposure(5)), EXPOSURE_KEYS,
      "reanalyzeExposure exposure cache");
  });

  // -------------------------------------------------------------------------
  // Group / index mutators
  // -------------------------------------------------------------------------

  it("addIndexToGroup writes a GroupEntry with exactly 5 keys (issue #16)", async () => {
    qc.setQueryData(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "auto", active: true, members: [] },
    ]);
    mockFetchOnce({
      id: 1, exposure_id: 5, kind: "custom", active: true, members: [42],
      event_id: 11, view_row_id: 1,
    }, 200);
    await runMutator(qc, addIndexToGroupMutator, {
      kind: "index_confirmed",
      clientOpId: "op-shape-6",
      exposureId: 5, groupId: 1, username: "alice", clientId: "tab-1",
      indexId: 42, payload: { groupId: 1, indexId: 42 },
    });
    const groups = qc.getQueryData<unknown[]>(queryKeys.groups(5));
    assertKeys(groups![0], GROUP_KEYS, "addIndexToGroup cache row");
  });

  it("removeIndexFromGroup writes a GroupEntry with exactly 5 keys", async () => {
    qc.setQueryData(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "custom", active: true, members: [42] },
    ]);
    mockFetchOnce({
      id: 1, exposure_id: 5, kind: "custom", active: true, members: [],
      event_id: 12, view_row_id: 1,
    }, 200);
    await runMutator(qc, removeIndexFromGroupMutator, {
      kind: "index_unconfirmed",
      clientOpId: "op-shape-7",
      exposureId: 5, groupId: 1, username: "alice", clientId: "tab-1",
      indexId: 42, payload: { groupId: 1, indexId: 42 },
    });
    const groups = qc.getQueryData<unknown[]>(queryKeys.groups(5));
    assertKeys(groups![0], GROUP_KEYS, "removeIndexFromGroup cache row");
  });

  it("createSpeculative writes an IndexEntry with exactly 13 keys", async () => {
    qc.setQueryData(queryKeys.indices(5), []);
    qc.setQueryData(queryKeys.groups(5), []);
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
  // Sample / message / tag mutators
  // -------------------------------------------------------------------------

  it("updateSample preserves tags via field-merge (deep-scan Bug #2)", async () => {
    const initialSample = {
      id: 10, experiment_id: 1, display_name: "D1", name: "old", notes: "n",
      tags: [{ id: 1, key: "k", value: "v", source: "manual" }],
    };
    qc.setQueryData(queryKeys.sample(10), initialSample);
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
    const single = qc.getQueryData<{ tags: unknown[] }>(queryKeys.sample(10));
    assertKeys(single!, SAMPLE_KEYS, "updateSample single cache row");
    expect(single!.tags).toHaveLength(1);
  });

  it("addSampleTag inserts a SampleTag with exactly 4 keys (no sample_id leak)", async () => {
    // The route emits {id, sample_id, key, value, source} (routes_samples.jl:85-87)
    // but the SampleTag type has only {id, key, value, source}. If the mutator
    // spreads the response wholesale, the cached tag pollutes with sample_id.
    const initialSample = {
      id: 10, experiment_id: 1, display_name: "D1", name: "n", notes: null, tags: [],
    };
    qc.setQueryData(queryKeys.samples(1), [initialSample]);
    mockFetchOnce({
      id: 50, sample_id: 10, key: "color", value: "red", source: "manual",
    }, 201);
    await runMutator(qc, addSampleTagMutator, {
      kind: "add_tag",
      clientOpId: "op-shape-10",
      sampleId: 10, experimentId: 1, username: "alice", clientId: "tab-1",
      key: "color", value: "red",
      payload: { key: "color", value: "red" },
    });
    const list = qc.getQueryData<{ tags: unknown[] }[]>(queryKeys.samples(1));
    const tag = list![0]!.tags[0];
    assertKeys(tag, SAMPLE_TAG_KEYS, "addSampleTag cache row");
  });

  it("addExposureTag inserts an ExposureTag with exactly 4 keys (no exposure_id leak)", async () => {
    // The route emits {id, exposure_id, key, value, source}; ExposureTag is
    // only {id, key, value, source}.
    const initialExposure = {
      id: 5, sample_id: 1, filename: null, kind: "file", selected: true,
      status: null, image_path: null, trace_hash: null,
      analysis_inputs_hash: "h0", tags: [], sources: [], image_version: "",
    };
    qc.setQueryData(
      ["sample", 1, "exposures"] as const,
      [initialExposure],
    );
    mockFetchOnce({
      id: 60, exposure_id: 5, key: "noisy", value: "yes", source: "manual",
    }, 201);
    await runMutator(qc, addExposureTagMutator, {
      kind: "add_tag",
      clientOpId: "op-shape-11",
      sampleId: 1, exposureId: 5, username: "alice", clientId: "tab-1",
      key: "noisy", value: "yes",
      payload: { key: "noisy", value: "yes" },
    });
    const list = qc.getQueryData<{ tags: unknown[] }[]>(
      ["sample", 1, "exposures"] as const);
    const tag = list![0]!.tags[0];
    assertKeys(tag, EXPOSURE_TAG_KEYS, "addExposureTag cache row");
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
    qc.setQueryData(queryKeys.exposure(5), FULL_EXPOSURE);
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
    qc.setQueryData(queryKeys.exposure(5), FULL_EXPOSURE);
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

  // -------------------------------------------------------------------------
  // Compare page mutators (Phase 3). 3 event-shape rows: comparison_created
  // (saveComparison create), comparison_submitted (saveComparison update),
  // comparison_deleted (deleteComparison).
  // -------------------------------------------------------------------------

  function buildComparisonResponse(id: number, hash = "sha256:abc") {
    return {
      id, title: "X", description: null, content_hash: hash,
      created_by: 1, created_at: "2026-05-06",
      updated_at: "2026-05-06",
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      members: [
        {
          id: 999, comparison_id: id, exposure_id: 100, display_order: 0,
          band_height: 1, y_offset: 0, normalization: "none",
          color_override: null, label_override: null,
          q_window_min: null, q_window_max: null,
          peak_display: null,
          snapshot: { effective_peaks: [], confirmed_index: null,
                      analysis_inputs_hash: "sha256:zero" },
          is_stale: false, created_by: 1, created_at: "2026-05-06",
        },
      ],
    };
  }

  it("saveComparison (create → comparison_created) writes a Comparison with exactly 11 keys", async () => {
    mockFetchOnce(buildComparisonResponse(42), 201);
    await runMutator(qc, saveComparisonMutator, {
      kind: "comparison_save", clientOpId: "op-cmp-create",
      username: "alice", clientId: "tab-1",
      title: "X",
      members: [{ exposure_id: 100, display_order: 0,
                  snapshot: { effective_peaks: [], confirmed_index: null,
                              analysis_inputs_hash: "sha256:zero" } }],
      payload: {},
    });
    assertKeys(qc.getQueryData(queryKeys.comparison(42)), COMPARISON_KEYS,
      "saveComparison(create) cache row");
    const members = qc.getQueryData<unknown[]>(queryKeys.comparisonMembers(42))!;
    assertKeys(members[0], COMPARISON_MEMBER_KEYS,
      "saveComparison(create) member cache row");
  });

  it("saveComparison (submit → comparison_submitted) writes a Comparison with exactly 11 keys", async () => {
    mockFetchOnce(buildComparisonResponse(42, "sha256:new"), 200);
    await runMutator(qc, saveComparisonMutator, {
      kind: "comparison_save", clientOpId: "op-cmp-submit",
      username: "alice", clientId: "tab-1",
      id: 42, title: "X edited",
      members: [{ id: 999, exposure_id: 100, display_order: 0,
                  snapshot: { effective_peaks: [], confirmed_index: null,
                              analysis_inputs_hash: "sha256:zero" } }],
      expected_content_hash: "sha256:abc",
      payload: {},
    });
    assertKeys(qc.getQueryData(queryKeys.comparison(42)), COMPARISON_KEYS,
      "saveComparison(submit) cache row");
  });

  it("deleteComparison (→ comparison_deleted) removes entity caches and prunes listings", async () => {
    qc.setQueryData(queryKeys.comparison(42), { id: 42 });
    qc.setQueryData(queryKeys.comparisonMembers(42), []);
    qc.setQueryData(queryKeys.comparisons("all"), [
      { id: 42, title: "doomed", description: null, content_hash: "h",
        created_by: 1, created_at: null, updated_at: null,
        forked_from_id: null, forked_at_hash: null },
      { id: 99, title: "kept", description: null, content_hash: "h2",
        created_by: 1, created_at: null, updated_at: null,
        forked_from_id: null, forked_at_hash: null },
    ]);
    mockFetchOnce({ id: 42, deleted: true, event_id: 7 }, 200);
    await runMutator(qc, deleteComparisonMutator, {
      kind: "comparison_delete", clientOpId: "op-cmp-del",
      username: "alice", clientId: "tab-1",
      id: 42, payload: { id: 42 },
    });
    // Removed (no cache entry, no error state)
    expect(qc.getQueryState(queryKeys.comparison(42))).toBeUndefined();
    expect(qc.getQueryState(queryKeys.comparisonMembers(42))).toBeUndefined();
    // Listing filtered down
    const listing = qc.getQueryData<{ id: number }[]>(queryKeys.comparisons("all"))!;
    expect(listing.map((c) => c.id)).toEqual([99]);
  });

  it("postSampleMessage writes a SampleMessage with exactly 6 keys", async () => {
    qc.setQueryData(queryKeys.messages(10), []);
    // Mock derived from routes_messages.jl POST /samples/:id/messages response.
    mockFetchOnce({
      id: 200, sample_id: 10, author_id: 3, author: "alice",
      body: "hello", created_at: "2026-05-03T12:00:00Z",
    }, 201);
    await runMutator(qc, postSampleMessageMutator, {
      kind: "post_message",
      clientOpId: "op-shape-12",
      sampleId: 10, username: "alice", clientId: "tab-1",
      body: "hello",
      payload: { body: "hello" },
    });
    const list = qc.getQueryData<unknown[]>(queryKeys.messages(10));
    expect(list).toHaveLength(1);
    assertKeys(list![0], SAMPLE_MESSAGE_KEYS, "postSampleMessage cache row");
  });
});

// ---------------------------------------------------------------------------
// Compare UX A-8 — ComparisonSummary listing projection shape.
//
// `_comparison_listing_rows` (backend, #137) emits denormalised projection
// fields (author_username, member_count, member_phases, has_stale_members)
// plus the persisted view-choice columns (view_*). These tests pin the shape
// at the value layer — a compile-time `: ComparisonSummary` anchor only
// catches a missing field, not a rename that leaves the old name aliased.
// ---------------------------------------------------------------------------
describe("comparison listing cache shape — Compare UX A-8", () => {
  it("includes the new projection fields", () => {
    const row: ComparisonSummary = {
      id: 1,
      title: "x",
      description: null,
      content_hash: "h",
      created_by: 1,
      created_at: null,
      updated_at: null,
      forked_from_id: null,
      forked_at_hash: null,
      view_grouping_mode: null,
      view_show_peak_ticks: null,
      view_show_peak_labels: null,
      last_event_at: "2026-05-14T10:00:00Z",
      author_username: "alice",
      member_count: 3,
      member_phases: ["Pn3m", "Hex"],
      has_stale_members: false,
    };
    expect(row.member_count).toBe(3);
    expect(row.member_phases).toEqual(["Pn3m", "Hex"]);
    expect(row.author_username).toBe("alice");
    // Pin each new view_* field at the value layer — protects against a
    // future refactor that renames (e.g. viewGroupingMode) while leaving
    // the old name aliased; compile-time-only checks slip through.
    expect(row.view_grouping_mode).toBeNull();
    expect(row.view_show_peak_ticks).toBeNull();
    expect(row.view_show_peak_labels).toBeNull();
    expect(row.last_event_at).toBe("2026-05-14T10:00:00Z");
    expect(row.has_stale_members).toBe(false);
  });

  it("accepts populated view_* values too", () => {
    const row: ComparisonSummary = {
      id: 2, title: "y", description: null, content_hash: "h2",
      created_by: 1, created_at: null, updated_at: null,
      forked_from_id: null, forked_at_hash: null,
      view_grouping_mode: "byPhase",
      view_show_peak_ticks: true,
      view_show_peak_labels: false,
      last_event_at: null, author_username: null,
      member_count: 0, member_phases: [], has_stale_members: false,
    };
    expect(row.view_grouping_mode).toBe("byPhase");
    expect(row.view_show_peak_ticks).toBe(true);
    expect(row.view_show_peak_labels).toBe(false);
  });
});
