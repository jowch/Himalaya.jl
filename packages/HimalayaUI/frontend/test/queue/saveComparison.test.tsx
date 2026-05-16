/**
 * saveComparison mutator (Plan §Phase 3, Task 3.1).
 *
 * The OpKind-shaped contract test row. Pins:
 * - request routes to POST /api/comparisons (create) when id is absent
 * - request routes to POST /api/comparisons/:id/submit when id is present
 * - 409 from the server raises a typed ConflictError carrying current state
 * - onSuccess writes the canonical comparison + members + invalidates listings
 * - 403 surfaces as a non-retry validation error (no rollback to fix; no
 *   optimistic effect to roll back, but the contract still applies to error
 *   classification)
 * - onMutate is a no-op (no optimistic write per spec)
 */
import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { renderHook, act, waitFor } from "@testing-library/react";
import type { ReactNode } from "react";
import { saveComparisonMutator } from "../../src/lib/queue/mutators/saveComparison";
import { resolveMutator } from "../../src/lib/queue/mutatorRegistry";
import { useQueueMutation } from "../../src/lib/queue/useQueueMutation";
import { setToastImpl } from "../../src/lib/toast";
import { pendingDeferreds } from "../../src/lib/queue/deferred";
import { ConflictError } from "../../src/api";
import type { Comparison, ComparisonMemberInput } from "../../src/api";
import { queryKeys, useSaveComparison } from "../../src/queries";
import { useAppState } from "../../src/state";
import {
  attachConflictBridge, _resetConflictBridgeForTest,
} from "../../src/lib/queue/conflictBridge";
import { makeClient } from "../test-utils";

const SNAP = {
  effective_peaks: [],
  confirmed_index: null,
  analysis_inputs_hash: "sha256:zero",
};

const MEMBER_INPUT: ComparisonMemberInput = {
  exposure_id: 100,
  display_order: 0,
  snapshot: SNAP,
};

function buildComparison(id: number, hash = "sha256:abc"): Comparison {
  return {
    id, title: "X", description: null, content_hash: hash,
    created_by: 1, created_at: "2026-05-06T00:00:00Z",
    updated_at: "2026-05-06T00:00:00Z",
    forked_from_id: null, forked_at_hash: null, forked_from_title: null,
    members: [
      {
        id: 999, comparison_id: id, exposure_id: 100, display_order: 0,
        band_height: 1, y_offset: 0, normalization: "none",
        color_override: null, label_override: null,
        q_window_min: null, q_window_max: null,
        peak_display: null, snapshot: SNAP, is_stale: false,
        created_by: 1, created_at: "2026-05-06T00:00:00Z",
      },
    ],
  };
}

interface Captured {
  url: string;
  method: string;
  body: string | undefined;
}

function captureFetch(responseBody: unknown, status = 200): Captured[] {
  const captured: Captured[] = [];
  globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
    const url = typeof input === "string" ? input : String(input);
    const method = init?.method ?? "GET";
    const body = typeof init?.body === "string" ? init.body : undefined;
    captured.push({ url, method, body });
    return new Response(JSON.stringify(responseBody), {
      status, headers: { "Content-Type": "application/json" },
    });
  }) as typeof fetch;
  return captured;
}

describe("saveComparisonMutator (OpKind=comparison_save)", () => {
  let qc: QueryClient;
  const originalFetch = globalThis.fetch;
  beforeEach(() => { qc = new QueryClient(); });
  afterEach(() => { globalThis.fetch = originalFetch; vi.restoreAllMocks(); });

  it("kind is 'comparison_save'", () => {
    expect(saveComparisonMutator.kind).toBe("comparison_save");
  });

  it("onMutate is a no-op (returns a noop restore — submission shows a spinner)", () => {
    const ctx = saveComparisonMutator.onMutate(
      { kind: "comparison_save", clientOpId: "op-1",
        username: "alice", clientId: "tab",
        title: "T", members: [], payload: { title: "T", members: [] } } as any,
      qc,
    );
    expect(typeof ctx.restore).toBe("function");
    // restore() is a noop — must not throw, must not mutate cache
    qc.setQueryData(queryKeys.comparison(1), { id: 1 });
    ctx.restore();
    expect(qc.getQueryData(queryKeys.comparison(1))).toEqual({ id: 1 });
  });

  it("create flow: POSTs /api/comparisons when id is absent", async () => {
    const captured = captureFetch(buildComparison(42), 201);
    await saveComparisonMutator.request(
      { kind: "comparison_save", clientOpId: "op-create",
        username: "alice", clientId: "tab",
        title: "Hello", description: "world",
        members: [MEMBER_INPUT],
        payload: { title: "Hello", members: [MEMBER_INPUT] } } as any,
      new AbortController().signal,
    );
    expect(captured).toHaveLength(1);
    expect(captured[0]!.url).toBe("/api/comparisons");
    expect(captured[0]!.method).toBe("POST");
    const body = JSON.parse(captured[0]!.body!);
    expect(body.title).toBe("Hello");
    expect(body.description).toBe("world");
    expect(body.members).toEqual([MEMBER_INPUT]);
    // expected_content_hash MUST be absent on create — first submission has
    // nothing to expect against.
    expect("expected_content_hash" in body).toBe(false);
  });

  it("submit flow: POSTs /api/comparisons/:id/submit when id is present", async () => {
    const captured = captureFetch(buildComparison(42), 200);
    await saveComparisonMutator.request(
      { kind: "comparison_save", clientOpId: "op-submit",
        username: "alice", clientId: "tab",
        id: 42, title: "Edited", description: null,
        members: [MEMBER_INPUT],
        expected_content_hash: "sha256:base",
        payload: { id: 42, title: "Edited", members: [MEMBER_INPUT],
                   expected_content_hash: "sha256:base" } } as any,
      new AbortController().signal,
    );
    expect(captured).toHaveLength(1);
    expect(captured[0]!.url).toBe("/api/comparisons/42/submit");
    expect(captured[0]!.method).toBe("POST");
    const body = JSON.parse(captured[0]!.body!);
    expect(body.expected_content_hash).toBe("sha256:base");
  });

  it("409 response throws a ConflictError carrying current_hash + current_state", async () => {
    const conflictBody = {
      error: "conflict",
      current_hash: "sha256:server",
      current_state: buildComparison(42, "sha256:server"),
    };
    captureFetch(conflictBody, 409);
    await expect(saveComparisonMutator.request(
      { kind: "comparison_save", clientOpId: "op-conf",
        username: "alice", clientId: "tab",
        id: 42, title: "Edited", members: [MEMBER_INPUT],
        expected_content_hash: "sha256:stale",
        payload: { id: 42, title: "Edited", members: [MEMBER_INPUT],
                   expected_content_hash: "sha256:stale" } } as any,
      new AbortController().signal,
    )).rejects.toMatchObject({
      name: "ConflictError",
      status: 409,
      current_hash: "sha256:server",
    });
    // And the typed throw carries the parsed Comparison body.
    try {
      await saveComparisonMutator.request(
        { kind: "comparison_save", clientOpId: "op-conf2",
          username: "alice", clientId: "tab",
          id: 42, title: "Edited", members: [MEMBER_INPUT],
          expected_content_hash: "sha256:stale",
          payload: {} } as any,
        new AbortController().signal,
      );
    } catch (e) {
      expect(e).toBeInstanceOf(ConflictError);
      expect((e as ConflictError).current_state?.id).toBe(42);
    }
  });

  it("403 surfaces as ApiError with status=403 (validation-class)", async () => {
    captureFetch({ error: "only the author can submit" }, 403);
    await expect(saveComparisonMutator.request(
      { kind: "comparison_save", clientOpId: "op-403",
        username: "bob", clientId: "tab",
        id: 42, title: "Edited", members: [MEMBER_INPUT],
        expected_content_hash: "sha256:x",
        payload: {} } as any,
      new AbortController().signal,
    )).rejects.toMatchObject({ status: 403 });
  });

  // Fix 3 (spec auditor): 404 mirrors the 403 path — the existing comparison
  // was deleted between edit-mode entry and submit. Surfaces as a non-retry
  // ApiError (NOT a ConflictError, so the conflict modal stays closed and the
  // failure-class router classifies it as validation per `isValidationError`).
  it("404 on submit surfaces as ApiError with status=404 (non-retry, not ConflictError)", async () => {
    captureFetch({ error: "comparison not found" }, 404);
    let thrown: unknown;
    try {
      await saveComparisonMutator.request(
        { kind: "comparison_save", clientOpId: "op-404",
          username: "alice", clientId: "tab",
          id: 42, title: "Edited", members: [MEMBER_INPUT],
          expected_content_hash: "sha256:x",
          payload: {} } as any,
        new AbortController().signal,
      );
    } catch (e) {
      thrown = e;
    }
    expect(thrown).toBeDefined();
    expect((thrown as { status?: number }).status).toBe(404);
    // Critical: the 404 path must NOT throw a ConflictError. The conflict
    // modal is wired off `error instanceof ConflictError` — a stray Conflict
    // throw on the deleted-comparison path opens an empty/null conflict UI.
    expect(thrown).not.toBeInstanceOf(ConflictError);
  });

  it("onSuccess writes comparison + members into cache and invalidates listings", () => {
    const response = buildComparison(42);
    const invalidatedKeys: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: any) => {
      invalidatedKeys.push(arg?.queryKey);
      return orig(arg);
    }) as typeof qc.invalidateQueries;
    saveComparisonMutator.onSuccess(
      { kind: "comparison_save", clientOpId: "op-1",
        username: "alice", clientId: "tab",
        id: 42, title: "X", members: [MEMBER_INPUT],
        payload: {} } as any,
      response, qc,
    );
    expect(qc.getQueryData(queryKeys.comparison(42))).toEqual(response);
    expect(qc.getQueryData(queryKeys.comparisonMembers(42))).toEqual(response.members);
    // Prefix invalidation hits both the global "all" key and any per-experiment
    // key — pin via the prefix shape.
    expect(invalidatedKeys).toEqual([["comparisons"]]);
  });

  it("onSuccess works on create (no prior cache state to overwrite)", () => {
    const response = buildComparison(99);
    saveComparisonMutator.onSuccess(
      { kind: "comparison_save", clientOpId: "op-1",
        username: "alice", clientId: "tab",
        title: "X", members: [MEMBER_INPUT],
        payload: {} } as any,
      response, qc,
    );
    expect(qc.getQueryData(queryKeys.comparison(99))).toEqual(response);
  });
});

// Fix 1 (queue-reviewer): a 409 ConflictError is a validation-class error
// (status 400-499), but it has its own bespoke surface — the conflict modal,
// driven off `useMutation.error instanceof ConflictError`. Without the
// explicit instanceof gate in `useQueueMutation::onError`, the validation-
// toast branch fires AND the modal opens — duplicate "Couldn't save
// comparison" copy stacked behind the diff UI.
describe("saveComparisonMutator — useQueueMutation suppresses toast on 409 ConflictError", () => {
  let toastCalls: Array<{ msg: string; kind: string }> = [];

  beforeEach(() => {
    vi.restoreAllMocks();
    pendingDeferreds.clear();
    toastCalls = [];
    setToastImpl((msg, kind) => { toastCalls.push({ msg, kind }); });
  });
  afterEach(() => {
    setToastImpl(null);
  });

  function withClient() {
    const client = makeClient();
    const wrapper = ({ children }: { children: ReactNode }) => (
      <QueryClientProvider client={client}>{children}</QueryClientProvider>
    );
    return { client, wrapper };
  }

  function mockFetch409(): void {
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({
        error: "conflict",
        current_hash: "sha256:server",
        current_state: buildComparison(42, "sha256:server"),
      }), {
        status: 409, headers: { "Content-Type": "application/json" },
      }),
    );
  }

  it("does not surface a validation toast on 409 (modal handles it)", async () => {
    const { wrapper } = withClient();
    mockFetch409();
    const { result } = renderHook(
      () => useQueueMutation(saveComparisonMutator, {
        username: "alice", clientId: "tab",
      }),
      { wrapper },
    );
    act(() => {
      result.current.mutate({
        id: 42, title: "Edited", members: [MEMBER_INPUT],
        expected_content_hash: "sha256:stale",
      });
    });
    // Wait long enough for onError to have fired if it were going to.
    await waitFor(() => expect(result.current.isPending).toBe(false));
    expect(toastCalls).toEqual([]);
    // Sanity: the typed throw still landed on `error` for the modal to read.
    expect(result.current.error).toBeInstanceOf(ConflictError);
  });

  it("on 409, the App-level conflictBridge populates Zustand `pendingConflict` (Phase 12 bridge)", async () => {
    const { client, wrapper } = withClient();
    mockFetch409();
    // Reset slot + bridge state in case other tests left them set.
    useAppState.setState({ pendingConflict: null });
    _resetConflictBridgeForTest();
    // Mount the bridge for this test, mirroring App.tsx. Lifting it out of
    // useSaveComparison removed the per-hook race; the contract pinned here
    // is "a 409 on a comparison_save mutation lands in the slot via the
    // App-level subscriber" — independent of how many places mount the hook.
    const detachBridge = attachConflictBridge(
      client.getMutationCache(),
      useAppState.getState().setPendingConflict,
    );
    const { result } = renderHook(() => useSaveComparison(), { wrapper });
    act(() => {
      result.current.mutate({
        id: 42, title: "Edited", members: [MEMBER_INPUT],
        expected_content_hash: "sha256:stale",
      });
    });
    await waitFor(() => expect(result.current.isPending).toBe(false));
    // Bridge fired: the typed throw is now in Zustand for the modal to read.
    await waitFor(() => {
      const slot = useAppState.getState().pendingConflict;
      expect(slot).toBeInstanceOf(ConflictError);
    });
    const slot = useAppState.getState().pendingConflict;
    expect(slot?.current_hash).toBe("sha256:server");
    expect(slot?.current_state?.id).toBe(42);
    // Toast still suppressed (regression on Phase 3 follow-up).
    expect(toastCalls).toEqual([]);
    // Cleanup so subsequent tests start clean.
    detachBridge();
    useAppState.setState({ pendingConflict: null });
  });

  it("still surfaces a validation toast on a non-Conflict 4xx (e.g. 403)", async () => {
    const { wrapper } = withClient();
    vi.spyOn(global, "fetch").mockResolvedValue(
      new Response(JSON.stringify({ error: "only the author can submit" }), {
        status: 403, headers: { "Content-Type": "application/json" },
      }),
    );
    const { result } = renderHook(
      () => useQueueMutation(saveComparisonMutator, {
        username: "bob", clientId: "tab",
      }),
      { wrapper },
    );
    act(() => {
      result.current.mutate({
        id: 42, title: "Edited", members: [MEMBER_INPUT],
        expected_content_hash: "sha256:x",
      });
    });
    await waitFor(() => expect(result.current.isPending).toBe(false));
    // Pin: the toast did fire for the non-Conflict validation case so the
    // suppression is narrowly scoped to ConflictError, not blanket-disabling
    // toasts for the saveComparison mutator.
    expect(toastCalls).toHaveLength(1);
    expect(toastCalls[0]!.kind).toBe("error");
    expect(toastCalls[0]!.msg).toContain("Couldn't save comparison");
  });
});

// ---------------------------------------------------------------------------
// Compare UX A-10 — view_* fields ride the saveComparison mutator.
//
// Three layers of the six-layer contract (per docs/contract-testing.md):
//  - request body: view_* survive `buildBody` onto the HTTP payload.
//  - onSuccess (layer 6): the HTTP response writes view_* into BOTH cache
//    keys — pins the HTTP-response-wins-race path (A-9 pins the SSE path).
//  - registry: widening the input type does not break replay routing.
// ---------------------------------------------------------------------------
describe("saveComparison passes view_* into request body — Compare UX A-10", () => {
  const originalFetch = globalThis.fetch;
  afterEach(() => { globalThis.fetch = originalFetch; vi.restoreAllMocks(); });

  it("forwards view_grouping_mode and friends", async () => {
    const spy = vi.spyOn(globalThis, "fetch").mockResolvedValue(
      new Response(
        JSON.stringify({ id: 1, members: [], content_hash: "h" }),
        { status: 200, headers: { "content-type": "application/json" } },
      ));
    await saveComparisonMutator.request(
      { kind: "comparison_save", clientOpId: "op",
        username: "alice", clientId: "tab",
        id: 1, title: "t", members: [],
        expected_content_hash: "h0",
        view_grouping_mode: "byPhase",
        view_show_peak_ticks: true,
        view_show_peak_labels: false,
        payload: {} } as any,
      new AbortController().signal,
    );
    const init = spy.mock.calls[0]?.[1];
    const body = JSON.parse(String(init?.body ?? "{}"));
    expect(body.view_grouping_mode).toBe("byPhase");
    expect(body.view_show_peak_ticks).toBe(true);
    expect(body.view_show_peak_labels).toBe(false);
  });
});

describe("saveComparison onSuccess reconciliation — Compare UX A-10", () => {
  it("writes view_* fields AND members into both cache keys from the HTTP response", () => {
    const qc = new QueryClient();
    const fakeMembers = [{ exposure_id: 100, display_order: 0 }] as any;
    const fakeResponse = {
      id: 1, title: "t", description: null, content_hash: "h1",
      created_by: 1, created_at: null, updated_at: null,
      forked_from_id: null, forked_at_hash: null, forked_from_title: null,
      view_grouping_mode: "byPhase",
      view_show_peak_ticks: true,
      view_show_peak_labels: false,
      last_event_at: null, members: fakeMembers,
    };
    saveComparisonMutator.onSuccess(
      { kind: "comparison_save", clientOpId: "op",
        username: "alice", clientId: "tab",
        id: 1, title: "t", members: [],
        view_grouping_mode: "byPhase",
        view_show_peak_ticks: true,
        view_show_peak_labels: false,
        payload: {} } as any,
      fakeResponse as any,
      qc,
    );
    // Layer-6 contract: BOTH cache keys must reflect the response.
    const cachedComparison = qc.getQueryData(queryKeys.comparison(1)) as typeof fakeResponse;
    expect(cachedComparison.view_grouping_mode).toBe("byPhase");
    expect(cachedComparison.view_show_peak_ticks).toBe(true);
    expect(cachedComparison.view_show_peak_labels).toBe(false);
    expect(qc.getQueryData(queryKeys.comparisonMembers(1))).toEqual(fakeMembers);
  });
});

describe("mutatorRegistry resolves saveComparison with view_* — Compare UX A-10", () => {
  it("returns saveComparisonMutator for a comparison_save op carrying view_*", () => {
    // The registry keys on OpKind (`comparison_save`), not the input shape —
    // this pins that widening the input type with view_* does not break
    // replay-as-rerun routing.
    const op = {
      kind: "comparison_save" as const,
      payload: {
        id: 1, title: "t", members: [],
        view_grouping_mode: "byPhase",
        view_show_peak_ticks: true,
        view_show_peak_labels: false,
      },
    };
    expect(resolveMutator(op)).toBe(saveComparisonMutator);
  });
});
