/**
 * Cache-shape integrity for mutators that spread a route response into a
 * typed query cache.
 *
 * The bug class this catches: a mutator's onSuccess writes the full
 * response into a cache typed as `Foo`, but the response has extra
 * (event_id, view_row_id, etc.) or missing fields. TypeScript can't
 * see this at runtime — issue #16 (group response pollution) was caught
 * only by reading the code; sixth-pass issue #17 (Peak missing
 * intensity/etc.) was caught only by re-running the route. This test
 * runs each cache-writing mutator end-to-end against a route-shaped
 * mock and asserts the cache row's keys match the type's exactly.
 *
 * Discipline: the mocked response shape MUST be derived from the actual
 * route's emit, not from the TypeScript interface (since the type might
 * be wrong). Each mock fixture below is annotated with the file:line of
 * the route handler it mirrors.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { peakAddMutator } from "../../src/lib/queue/mutators/peakAdd";
import {
  peakExcludeMutator, peakUnexcludeMutator,
} from "../../src/lib/queue/mutators/peakSetExcluded";
import {
  addIndexToGroupMutator, removeIndexFromGroupMutator,
} from "../../src/lib/queue/mutators/indexGroup";
import { updateSampleMutator } from "../../src/lib/queue/mutators/trivial";
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
  "id", "experiment_id", "label", "name", "notes", "tags",
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

async function runMutator<I, S, R>(
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
  // Both directions strict — extra keys mean cache pollution; missing
  // keys mean the route response omitted a type-required field.
  expect({ label, actual: [...actual].sort() })
    .toEqual({ label, actual: [...expected].sort() });
}

describe("Cache-shape integrity (mutator onSuccess writes type-shaped rows)", () => {
  let qc: QueryClient;
  beforeEach(() => {
    qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    pendingDeferreds.clear();
  });

  it("peakAdd writes a Peak with exactly the 8 declared keys (issue #17 + deep-scan #1)", async () => {
    qc.setQueryData(queryKeys.peaks(5), []);
    qc.setQueryData(queryKeys.exposure(5), {
      id: 5, sample_id: 1, filename: null, kind: "file", selected: true,
      status: null, image_path: null, trace_hash: null,
      analysis_inputs_hash: "h0", tags: [], sources: [], image_version: "",
    });
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
    const initialPeak = {
      id: 7, exposure_id: 5, q: 0.5, intensity: 1.2, prominence: 0.8,
      sharpness: 30.0, source: "auto", excluded: false,
    };
    qc.setQueryData(queryKeys.peaks(5), [initialPeak]);
    qc.setQueryData(queryKeys.exposure(5), {
      id: 5, sample_id: 1, filename: null, kind: "file", selected: true,
      status: null, image_path: null, trace_hash: null,
      analysis_inputs_hash: "h0", tags: [], sources: [], image_version: "",
    });
    // Mock derived from routes_peaks.jl PATCH /peaks/:id response (extends Peak).
    mockFetchOnce({
      ...initialPeak, excluded: true,
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
    const list = qc.getQueryData<unknown[]>(queryKeys.peaks(5));
    assertKeys((list as { id: number }[])!.find(p => p.id === 7)!,
      PEAK_KEYS, "peakExclude cache row");
  });

  it("addIndexToGroup writes a GroupEntry with exactly 5 keys (issue #16)", async () => {
    qc.setQueryData(queryKeys.groups(5), [
      { id: 1, exposure_id: 5, kind: "auto", active: true, members: [] },
    ]);
    // Mock derived from routes_analysis.jl _group_with_members + the
    // event_id/view_row_id metadata added in issue #13's fix.
    mockFetchOnce({
      id: 1, exposure_id: 5, kind: "custom", active: true, members: [42],
      event_id: 11, view_row_id: 1,
    }, 200);
    await runMutator(qc, addIndexToGroupMutator, {
      kind: "index_confirmed",
      clientOpId: "op-shape-3",
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
      clientOpId: "op-shape-4",
      exposureId: 5, groupId: 1, username: "alice", clientId: "tab-1",
      indexId: 42, payload: { groupId: 1, indexId: 42 },
    });
    const groups = qc.getQueryData<unknown[]>(queryKeys.groups(5));
    assertKeys(groups![0], GROUP_KEYS, "removeIndexFromGroup cache row");
  });

  it("updateSample preserves tags via field-merge (deep-scan Bug #2)", async () => {
    // Pre-populate the sample cache with a `tags` array that the PATCH
    // response will NOT include — the mutator must merge only the patched
    // fields, NOT spread the response wholesale (which would clobber
    // tags to undefined).
    const initialSample = {
      id: 10, experiment_id: 1, label: "D1", name: "old", notes: "n",
      tags: [{ id: 1, key: "k", value: "v", source: "manual" }],
    };
    qc.setQueryData(queryKeys.sample(10), initialSample);
    qc.setQueryData(queryKeys.samples(1), [initialSample]);
    // Mock derived from routes_samples.jl PATCH /samples/:id response —
    // bare samples row with NO tags. This is exactly the shape that
    // would clobber Sample.tags if the mutator spread the response.
    mockFetchOnce({
      id: 10, experiment_id: 1, label: "D1", name: "new", notes: "n",
      created_at: "2026-05-03",
    }, 200);
    await runMutator(qc, updateSampleMutator, {
      kind: "update_sample",
      clientOpId: "op-shape-5",
      sampleId: 10, experimentId: 1, username: "alice", clientId: "tab-1",
      name: "new",
      payload: { sampleId: 10 },
    });
    const single = qc.getQueryData<{ tags: unknown[] }>(queryKeys.sample(10));
    assertKeys(single!, SAMPLE_KEYS, "updateSample single cache row");
    // Tags survived — the bug was that a wholesale spread would set
    // tags = undefined on the cache.
    expect(single!.tags).toHaveLength(1);
  });
});
