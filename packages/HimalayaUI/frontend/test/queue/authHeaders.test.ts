/**
 * Auth-header propagation contract.
 *
 * Every mutator's HTTP request MUST carry all three audit/idempotency
 * headers: `X-Username`, `X-Client-Id`, `X-Client-Op-Id`.
 *
 * Why this matters:
 * - `X-Username`     → identifies the actor in `user_actions`. Missing →
 *   audit log says "unknown".
 * - `X-Client-Id`    → SSE self-echo filtering. Missing → the mutating tab
 *   re-applies its own SSE frame (double-apply / flicker).
 * - `X-Client-Op-Id` → idempotency token. Missing → retries re-execute
 *   instead of returning the cached response.
 *
 * Forgetting any one of these surfaces only as subtle multiplayer bugs
 * under concurrent load — code review has caught these so far, but a
 * mechanical test pins it as a contract.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { request } from "../../src/api";
import { peakAddMutator } from "../../src/lib/queue/mutators/peakAdd";
import { peakRemoveMutator } from "../../src/lib/queue/mutators/peakRemove";
import {
  peakExcludeMutator, peakUnexcludeMutator,
} from "../../src/lib/queue/mutators/peakSetExcluded";
import {
  addIndexToGroupMutator, removeIndexFromGroupMutator, deleteIndexMutator,
} from "../../src/lib/queue/mutators/indexGroup";
import { createSpeculativeMutator } from "../../src/lib/queue/mutators/createSpeculative";
import { reanalyzeExposureMutator } from "../../src/lib/queue/mutators/reanalyzeExposure";
import { saveComparisonMutator } from "../../src/lib/queue/mutators/saveComparison";
import { deleteComparisonMutator } from "../../src/lib/queue/mutators/deleteComparison";
import {
  updateSampleMutator,
  addSampleTagMutator, removeSampleTagMutator,
  addExposureTagMutator, removeExposureTagMutator,
  postSampleMessageMutator,
  setExposureStatusMutator, selectExposureMutator,
} from "../../src/lib/queue/mutators/trivial";

interface Captured {
  url: string;
  method: string;
  headers: Headers;
}

const REQUIRED_HEADERS = ["X-Username", "X-Client-Id", "X-Client-Op-Id"];

function captureFetch(): Captured[] {
  const captured: Captured[] = [];
  globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
    const url = typeof input === "string" ? input : String(input);
    const method = init?.method ?? "GET";
    const headers = new Headers(init?.headers ?? {});
    captured.push({ url, method, headers });
    // Generic body with ALL fields the various mutators read in onSuccess —
    // we only inspect outgoing headers, but onSuccess must not throw.
    return new Response(JSON.stringify({
      id: 100, exposure_id: 5, sample_id: 10, group_id: 1, index_id: 10,
      q: 0.5, intensity: null, prominence: null, sharpness: null,
      source: "manual", excluded: false, label: null,
      name: "n", notes: null, key: "k", value: "v", source_field: "manual",
      kind: "auto", active: true, members: [], peaks: [], predicted_q: [],
      basis: 0.1, score: 0.9, r_squared: 0.99, lattice_d: 50, ngc: 0.5,
      status: "candidate", inputs_hash: "h",
      author_id: 3, author: "alice", body: "hi", created_at: "2026-05-03",
      event_id: 9, view_row_id: 100, analysis_inputs_hash: "h1",
      analyzed: true, deleted: 1, selected: true,
    }), { status: 200, headers: { "Content-Type": "application/json" } });
  }) as typeof fetch;
  return captured;
}

interface Spec {
  name: string;
  run: (qc: QueryClient) => Promise<unknown>;
}

const FLAT_BASE = {
  clientOpId: "test-op-id",
  username: "alice",
  clientId: "tab-1",
} as const;

const SPECS: Spec[] = [
  {
    name: "peakAdd",
    run: (qc) => peakAddMutator.request(
      { ...FLAT_BASE, kind: "peak_added", payload: { q: 0.5 },
        exposureId: 5, q: 0.5 } as any,
      new AbortController().signal),
  },
  {
    name: "peakRemove",
    run: (qc) => peakRemoveMutator.request(
      { ...FLAT_BASE, kind: "peak_removed", payload: { peakId: 7 },
        exposureId: 5, peakId: 7 } as any,
      new AbortController().signal),
  },
  {
    name: "peakExclude",
    run: (qc) => peakExcludeMutator.request(
      { ...FLAT_BASE, kind: "peak_excluded",
        payload: { peakId: 7, q: 0.5 },
        exposureId: 5, peakId: 7, q: 0.5 } as any,
      new AbortController().signal),
  },
  {
    name: "peakUnexclude",
    run: (qc) => peakUnexcludeMutator.request(
      { ...FLAT_BASE, kind: "peak_unexcluded",
        payload: { peakId: 7, q: 0.5 },
        exposureId: 5, peakId: 7, q: 0.5 } as any,
      new AbortController().signal),
  },
  {
    name: "addIndexToGroup",
    run: (qc) => addIndexToGroupMutator.request(
      { ...FLAT_BASE, kind: "index_confirmed",
        payload: { groupId: 1, indexId: 10 },
        exposureId: 5, groupId: 1, indexId: 10 } as any,
      new AbortController().signal),
  },
  {
    name: "removeIndexFromGroup",
    run: (qc) => removeIndexFromGroupMutator.request(
      { ...FLAT_BASE, kind: "index_unconfirmed",
        payload: { groupId: 1, indexId: 10 },
        exposureId: 5, groupId: 1, indexId: 10 } as any,
      new AbortController().signal),
  },
  {
    name: "deleteIndex",
    run: (qc) => deleteIndexMutator.request(
      { ...FLAT_BASE, kind: "delete_index",
        payload: { indexId: 10 }, exposureId: 5, indexId: 10 } as any,
      new AbortController().signal),
  },
  {
    name: "createSpeculative",
    run: (qc) => createSpeculativeMutator.request(
      { ...FLAT_BASE, kind: "speculative_created",
        payload: { phase: "Pn3m", anchor_peak_id: 7, anchor_ratio: 1, additional: [] },
        exposureId: 5, phase: "Pn3m", anchor_peak_id: 7,
        anchor_ratio: 1, additional: [] } as any,
      new AbortController().signal),
  },
  {
    name: "reanalyzeExposure",
    run: (qc) => reanalyzeExposureMutator.request(
      { ...FLAT_BASE, kind: "reanalyze_exposure", payload: {},
        exposureId: 5 } as any,
      new AbortController().signal),
  },
  {
    name: "updateSample",
    run: (qc) => updateSampleMutator.request(
      { ...FLAT_BASE, kind: "update_sample",
        payload: { sampleId: 10 },
        sampleId: 10, experimentId: 1, display_name: "x" } as any,
      new AbortController().signal),
  },
  {
    name: "addSampleTag",
    run: (qc) => addSampleTagMutator.request(
      { ...FLAT_BASE, kind: "add_tag",
        payload: { key: "k", value: "v" },
        sampleId: 10, experimentId: 1, key: "k", value: "v" } as any,
      new AbortController().signal),
  },
  {
    name: "removeSampleTag",
    run: (qc) => removeSampleTagMutator.request(
      { ...FLAT_BASE, kind: "remove_tag",
        payload: { tagId: 1 },
        sampleId: 10, experimentId: 1, tagId: 1 } as any,
      new AbortController().signal),
  },
  {
    name: "addExposureTag",
    run: (qc) => addExposureTagMutator.request(
      { ...FLAT_BASE, kind: "add_tag",
        payload: { key: "k", value: "v" },
        sampleId: 1, exposureId: 5, key: "k", value: "v" } as any,
      new AbortController().signal),
  },
  {
    name: "removeExposureTag",
    run: (qc) => removeExposureTagMutator.request(
      { ...FLAT_BASE, kind: "remove_tag",
        payload: { tagId: 99 },
        sampleId: 1, exposureId: 5, tagId: 99 } as any,
      new AbortController().signal),
  },
  {
    name: "postSampleMessage",
    run: (qc) => postSampleMessageMutator.request(
      { ...FLAT_BASE, kind: "post_message", payload: { body: "hi" },
        sampleId: 10, body: "hi" } as any,
      new AbortController().signal),
  },
  {
    name: "setExposureStatus",
    run: (qc) => setExposureStatusMutator.request(
      { ...FLAT_BASE, kind: "set_exposure_status",
        payload: { exposureId: 5, status: "rejected" },
        sampleId: 1, exposureId: 5, status: "rejected" } as any,
      new AbortController().signal),
  },
  {
    name: "selectExposure",
    run: (qc) => selectExposureMutator.request(
      { ...FLAT_BASE, kind: "select_exposure",
        payload: { exposureId: 6 },
        sampleId: 1, exposureId: 6 } as any,
      new AbortController().signal),
  },
  // Compare page mutators (Phase 3). Both must carry the same audit +
  // idempotency headers as every other mutator.
  {
    name: "saveComparison",
    run: (qc) => saveComparisonMutator.request(
      { ...FLAT_BASE, kind: "comparison_save", payload: {},
        title: "X", description: null, members: [] } as any,
      new AbortController().signal),
  },
  {
    name: "deleteComparison",
    run: (qc) => deleteComparisonMutator.request(
      { ...FLAT_BASE, kind: "comparison_delete", payload: { id: 42 },
        id: 42 } as any,
      new AbortController().signal),
  },
];

describe("Auth-header propagation contract — every mutator carries audit + idempotency headers", () => {
  let originalFetch: typeof fetch;
  beforeEach(() => { originalFetch = globalThis.fetch; });
  for (const spec of SPECS) {
    it(`${spec.name} sends X-Username, X-Client-Id, X-Client-Op-Id`, async () => {
      const captured = captureFetch();
      try {
        const qc = new QueryClient();
        await spec.run(qc).catch(() => {/* response shape doesn't matter */});
        expect(captured.length).toBeGreaterThanOrEqual(1);
        const headers = captured[0]!.headers;
        for (const h of REQUIRED_HEADERS) {
          expect({ mutator: spec.name, header: h, value: headers.get(h) })
            .toEqual({ mutator: spec.name, header: h,
                       value: expect.stringMatching(/.+/) });
        }
        // Spot-check the actual values plumbed through.
        expect(headers.get("X-Username")).toBe("alice");
        expect(headers.get("X-Client-Id")).toBe("tab-1");
        expect(headers.get("X-Client-Op-Id")).toBe("test-op-id");
      } finally {
        globalThis.fetch = originalFetch;
      }
    });
  }
});

describe("Auth-header propagation contract — GET requests carry X-Username when opts provided", () => {
  let originalFetch: typeof fetch;
  beforeEach(() => { originalFetch = globalThis.fetch; });

  it("request() sends X-Username on GET when opts.username is set", async () => {
    const captured: { url: string; method: string; headers: Headers }[] = [];
    globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
      const url = typeof input === "string" ? input : String(input);
      const method = init?.method ?? "GET";
      const headers = new Headers(init?.headers ?? {});
      captured.push({ url, method, headers });
      return new Response(JSON.stringify({}), { status: 200, headers: { "Content-Type": "application/json" } });
    }) as typeof fetch;

    try {
      await request("GET", "/api/users/me/comparison-pins", undefined, {
        username: "alice",
        clientId: "tab-1",
        clientOpId: "test-op-id",
      }).catch(() => {/* response shape doesn't matter */});

      expect(captured.length).toBe(1);
      const headers = captured[0]!.headers;
      expect(headers.get("X-Username")).toBe("alice");
      expect(headers.get("X-Client-Id")).toBe("tab-1");
      expect(headers.get("X-Client-Op-Id")).toBe("test-op-id");
    } finally {
      globalThis.fetch = originalFetch;
    }
  });
});
