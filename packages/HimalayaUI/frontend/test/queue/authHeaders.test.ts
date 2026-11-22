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
import { deleteIndexMutator } from "../../src/lib/queue/mutators/indexGroup";
import { createSpeculativeMutator } from "../../src/lib/queue/mutators/createSpeculative";
import {
  addAssignmentPhaseMutator, removeAssignmentPhaseMutator, setSampleCallStateMutator,
} from "../../src/lib/queue/mutators/assignment";
import { customIndexMutator } from "../../src/lib/queue/mutators/customIndex";
import { reanalyzeExposureMutator } from "../../src/lib/queue/mutators/reanalyzeExposure";
import {
  updateSampleMutator,
  setExposureStatusMutator,
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
      basis: 0.1, score: 0.9, r_squared: 0.99, lattice_d: 50,
      status: "candidate", inputs_hash: "h",
      author_id: 3, author: "alice", body: "hi", created_at: "2026-05-03",
      event_id: 9, view_row_id: 100, analysis_inputs_hash: "h1",
      analyzed: true, deleted: 1,
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
    run: () => peakAddMutator.request(
      { ...FLAT_BASE, kind: "peak_added",
        exposureId: 5, q: 0.5 } as any,
      new AbortController().signal),
  },
  {
    name: "peakRemove",
    run: () => peakRemoveMutator.request(
      { ...FLAT_BASE, kind: "peak_removed",
        exposureId: 5, peakId: 7 } as any,
      new AbortController().signal),
  },
  {
    name: "peakExclude",
    run: () => peakExcludeMutator.request(
      { ...FLAT_BASE, kind: "peak_excluded",
        exposureId: 5, peakId: 7, q: 0.5 } as any,
      new AbortController().signal),
  },
  {
    name: "peakUnexclude",
    run: () => peakUnexcludeMutator.request(
      { ...FLAT_BASE, kind: "peak_unexcluded",
        exposureId: 5, peakId: 7, q: 0.5 } as any,
      new AbortController().signal),
  },
  {
    name: "deleteIndex",
    run: () => deleteIndexMutator.request(
      { ...FLAT_BASE, kind: "delete_index",
        exposureId: 5, indexId: 10 } as any,
      new AbortController().signal),
  },
  {
    name: "addAssignmentPhase",
    run: () => addAssignmentPhaseMutator.request(
      { ...FLAT_BASE, kind: "assignment_add",
        exposureId: 5, indexId: 10 } as any,
      new AbortController().signal),
  },
  {
    name: "removeAssignmentPhase",
    run: () => removeAssignmentPhaseMutator.request(
      { ...FLAT_BASE, kind: "assignment_remove",
        exposureId: 5, indexId: 10 } as any,
      new AbortController().signal),
  },
  {
    name: "setSampleCallState",
    run: () => setSampleCallStateMutator.request(
      { ...FLAT_BASE, kind: "sample_call_set_state",
        sampleId: 1, state: "form_factor" } as any,
      new AbortController().signal),
  },
  {
    name: "customIndex",
    run: () => customIndexMutator.request(
      { ...FLAT_BASE, kind: "custom_index_commit",
        exposureId: 5, phase: "Pn3m", basis: 0.15 } as any,
      new AbortController().signal),
  },
  {
    name: "createSpeculative",
    run: () => createSpeculativeMutator.request(
      { ...FLAT_BASE, kind: "speculative_created",
        exposureId: 5, phase: "Pn3m", anchor_peak_id: 7,
        anchor_ratio: 1, additional: [] } as any,
      new AbortController().signal),
  },
  {
    name: "reanalyzeExposure",
    run: () => reanalyzeExposureMutator.request(
      { ...FLAT_BASE, kind: "reanalyze_exposure",
        exposureId: 5 } as any,
      new AbortController().signal),
  },
  {
    name: "updateSample",
    run: () => updateSampleMutator.request(
      { ...FLAT_BASE, kind: "update_sample",
        sampleId: 10, experimentId: 1, display_name: "x" } as any,
      new AbortController().signal),
  },
  {
    name: "setExposureStatus",
    run: () => setExposureStatusMutator.request(
      { ...FLAT_BASE, kind: "set_exposure_status",
        sampleId: 1, exposureId: 5, status: "rejected" } as any,
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
