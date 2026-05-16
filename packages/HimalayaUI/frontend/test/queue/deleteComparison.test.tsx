/**
 * deleteComparison mutator (Plan §Phase 3, Task 3.2).
 *
 * Pins:
 * - request issues DELETE /api/comparisons/:id
 * - onSuccess removes the entity caches AND filters the id out of every
 *   cached listing key (NOT invalidates — refetching a deleted resource
 *   404s and leaves stale error state)
 * - 404 is treated as success (idempotent remove, treats404AsSuccess flag)
 * - 403 surfaces as ApiError with status=403
 */
import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { deleteComparisonMutator } from "../../src/lib/queue/mutators/deleteComparison";
import { queryKeys } from "../../src/queries";
import type { ComparisonSummary } from "../../src/api";

interface Captured { url: string; method: string; }

function captureFetch(responseBody: unknown, status = 200): Captured[] {
  const captured: Captured[] = [];
  globalThis.fetch = (async (input: RequestInfo | URL, init?: RequestInit) => {
    const url = typeof input === "string" ? input : String(input);
    captured.push({ url, method: init?.method ?? "GET" });
    return new Response(JSON.stringify(responseBody), {
      status, headers: { "Content-Type": "application/json" },
    });
  }) as typeof fetch;
  return captured;
}

function summary(id: number, title: string): ComparisonSummary {
  return {
    id, title, description: null, content_hash: "sha256:" + id,
    created_by: 1, created_at: null, updated_at: null,
    forked_from_id: null, forked_at_hash: null,
    view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: null,
    author_username: null, member_count: 0,
    member_phases: [], member_phase_count: 0, has_stale_members: false,
  };
}

describe("deleteComparisonMutator (OpKind=comparison_delete)", () => {
  let qc: QueryClient;
  const originalFetch = globalThis.fetch;
  beforeEach(() => { qc = new QueryClient(); });
  afterEach(() => { globalThis.fetch = originalFetch; vi.restoreAllMocks(); });

  it("kind is 'comparison_delete'", () => {
    expect(deleteComparisonMutator.kind).toBe("comparison_delete");
  });

  it("treats404AsSuccess is true (idempotent remove)", () => {
    expect(deleteComparisonMutator.treats404AsSuccess).toBe(true);
  });

  it("onMutate is a no-op", () => {
    qc.setQueryData(queryKeys.comparison(42), { id: 42 });
    const ctx = deleteComparisonMutator.onMutate(
      { kind: "comparison_delete", clientOpId: "op-1",
        username: "alice", clientId: "tab",
        id: 42, payload: { id: 42 } } as any,
      qc,
    );
    ctx.restore();
    // restore() must not have mutated the cache
    expect(qc.getQueryData(queryKeys.comparison(42))).toEqual({ id: 42 });
  });

  it("request issues DELETE /api/comparisons/:id", async () => {
    const captured = captureFetch({ id: 42, deleted: true, event_id: 7 }, 200);
    await deleteComparisonMutator.request(
      { kind: "comparison_delete", clientOpId: "op-1",
        username: "alice", clientId: "tab",
        id: 42, payload: { id: 42 } } as any,
      new AbortController().signal,
    );
    expect(captured).toHaveLength(1);
    expect(captured[0]!.url).toBe("/api/comparisons/42");
    expect(captured[0]!.method).toBe("DELETE");
  });

  it("onSuccess removes entity caches (NOT invalidates) and filters listings", () => {
    qc.setQueryData(queryKeys.comparison(42), { id: 42 });
    qc.setQueryData(queryKeys.comparisonMembers(42), []);
    qc.setQueryData(queryKeys.comparisonMessages(42), []);
    qc.setQueryData(queryKeys.comparisonForks(42), []);
    qc.setQueryData(queryKeys.comparisons("all"), [
      summary(42, "doomed"),
      summary(99, "kept"),
    ]);
    qc.setQueryData(queryKeys.comparisons(7), [summary(42, "doomed")]);
    deleteComparisonMutator.onSuccess(
      { kind: "comparison_delete", clientOpId: "op-1",
        username: "alice", clientId: "tab",
        id: 42, payload: { id: 42 } } as any,
      { id: 42, deleted: true, event_id: 7 }, qc,
    );
    // Removed (NOT invalidated — refetch would 404)
    expect(qc.getQueryState(queryKeys.comparison(42))).toBeUndefined();
    expect(qc.getQueryState(queryKeys.comparisonMembers(42))).toBeUndefined();
    expect(qc.getQueryState(queryKeys.comparisonMessages(42))).toBeUndefined();
    expect(qc.getQueryState(queryKeys.comparisonForks(42))).toBeUndefined();
    // Listings filtered (id pruned, others retained)
    expect(qc.getQueryData<ComparisonSummary[]>(queryKeys.comparisons("all")))
      .toEqual([summary(99, "kept")]);
    expect(qc.getQueryData<ComparisonSummary[]>(queryKeys.comparisons(7)))
      .toEqual([]);
  });

  it("403 surfaces as ApiError with status=403 (validation-class)", async () => {
    captureFetch({ error: "only the author can delete" }, 403);
    await expect(deleteComparisonMutator.request(
      { kind: "comparison_delete", clientOpId: "op-1",
        username: "bob", clientId: "tab",
        id: 42, payload: { id: 42 } } as any,
      new AbortController().signal,
    )).rejects.toMatchObject({ status: 403 });
  });
});
