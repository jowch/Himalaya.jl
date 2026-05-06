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
import { QueryClient } from "@tanstack/react-query";
import { saveComparisonMutator } from "../../src/lib/queue/mutators/saveComparison";
import { ConflictError } from "../../src/api";
import type { Comparison, ComparisonMemberInput } from "../../src/api";
import { queryKeys } from "../../src/queries";

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

  it("onSuccess writes comparison + members into cache and invalidates listings", () => {
    const response = buildComparison(42);
    let invalidatedKeys: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: { queryKey: unknown }) => {
      invalidatedKeys.push(arg.queryKey);
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
