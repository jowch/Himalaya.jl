/**
 * applyRemoteToCache — comparison_submitted / comparison_created splice
 * (Compare UX A-9).
 *
 * As of A-5 Step 5b the SSE frame for a comparison submit/create carries
 * `post_state = fetch_comparison_with_members(db, id)` — the same projection
 * `GET /api/comparisons/:id` returns. The branch under test splices that
 * post_state straight into `queryKeys.comparison(id)` (+ members) instead of
 * invalidating, so the new view_* fields land without a refetch round-trip.
 *
 * Note on signature: `applyRemoteToCache(remote, qc)` — the SSE frame is the
 * first argument, the QueryClient the second.
 */
import { describe, it, expect } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { applyRemoteToCache } from "../../src/lib/queue/applyRemoteToCache";
import { queryKeys } from "../../src/queries";
import type { SseEvent } from "../../src/lib/queue/types";
import type { Comparison } from "../../src/api";

function baselineComparison(cid: number): Comparison {
  return {
    id: cid, title: "x", description: null, content_hash: "h0",
    created_by: 1, created_at: null, updated_at: null,
    forked_from_id: null, forked_at_hash: null, forked_from_title: null,
    view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, last_event_at: null, members: [],
  };
}

describe("applyRemoteToCache: comparison_submitted view_* — Compare UX A-9", () => {
  it("writes view_* fields into the cache", () => {
    const qc = new QueryClient();
    const cid = 7;
    const baseline = baselineComparison(cid);
    qc.setQueryData(queryKeys.comparison(cid), baseline);

    const remote: SseEvent = {
      id: 99, kind: "comparison_submitted",
      entity_type: "comparison", entity_id: cid,
      actor: "bob", client_id: "tab2", client_op_id: "op2",
      ts: "2026-05-14T10:00:00Z",
      payload: null,
      post_state: {
        ...baseline,
        title: "y",
        content_hash: "h1",
        view_grouping_mode: "byPhase",
        view_show_peak_ticks: true,
        view_show_peak_labels: false,
      },
    };
    applyRemoteToCache(remote, qc);

    const updated = qc.getQueryData<Comparison>(queryKeys.comparison(cid));
    expect(updated?.view_grouping_mode).toBe("byPhase");
    expect(updated?.view_show_peak_ticks).toBe(true);
    expect(updated?.view_show_peak_labels).toBe(false);
  });

  it("preserves forked_from_title through the splice", () => {
    const qc = new QueryClient();
    const cid = 7;
    // forked_from_title is server-computed — the contract the previous
    // invalidate-only design existed to protect. Splicing post_state (which
    // IS the server projection) must keep it intact.
    const baseline: Comparison = {
      ...baselineComparison(cid),
      forked_from_id: 3, forked_at_hash: "h_parent",
      forked_from_title: "Parent",
    };
    qc.setQueryData(queryKeys.comparison(cid), baseline);

    const remote: SseEvent = {
      id: 99, kind: "comparison_submitted",
      entity_type: "comparison", entity_id: cid,
      actor: "bob", client_id: "tab2", client_op_id: "op2",
      ts: "2026-05-14T10:00:00Z", payload: null,
      post_state: { ...baseline, title: "renamed", content_hash: "h1" },
    };
    applyRemoteToCache(remote, qc);

    const updated = qc.getQueryData<Comparison>(queryKeys.comparison(cid));
    expect(updated?.forked_from_title).toBe("Parent");
    expect(updated?.title).toBe("renamed");
  });

  it("splices the members array from post_state", () => {
    const qc = new QueryClient();
    const cid = 8;
    const baseline = baselineComparison(cid);
    qc.setQueryData(queryKeys.comparison(cid), baseline);

    const member = {
      id: 1, comparison_id: cid, exposure_id: 100, display_order: 0,
      band_height: 1, y_offset: 0, normalization: "none",
      color_override: null, label_override: null,
      q_window_min: null, q_window_max: null,
      peak_display: null, snapshot: null, is_stale: false,
      created_by: 1, created_at: null,
    };
    const remote: SseEvent = {
      id: 100, kind: "comparison_created",
      entity_type: "comparison", entity_id: cid,
      ts: "2026-05-14T10:00:00Z", payload: null,
      post_state: { ...baseline, content_hash: "h1", members: [member] },
    };
    applyRemoteToCache(remote, qc);

    expect(qc.getQueryData(queryKeys.comparisonMembers(cid))).toEqual([member]);
  });

  it("falls back to invalidate when post_state is absent (pre-A-5 frame)", () => {
    const qc = new QueryClient();
    const cid = 9;
    const invalidated: unknown[] = [];
    const orig = qc.invalidateQueries.bind(qc);
    qc.invalidateQueries = ((arg: { queryKey?: unknown }) => {
      invalidated.push(arg?.queryKey);
      return orig(arg);
    }) as typeof qc.invalidateQueries;

    const remote: SseEvent = {
      id: 101, kind: "comparison_submitted",
      entity_type: "comparison", entity_id: cid,
      ts: "2026-05-14T10:00:00Z", payload: null,
      // no post_state
    };
    applyRemoteToCache(remote, qc);

    expect(invalidated).toContainEqual(queryKeys.comparison(cid));
    expect(invalidated).toContainEqual(queryKeys.comparisonMembers(cid));
    expect(invalidated).toContainEqual(["comparisons"]);
  });
});

describe("applyRemoteToCache comparison whole-row overwrite — Compare UX A-9", () => {
  it("foreign-wins final state with an optimistic own-op write in the cache", () => {
    const qc = new QueryClient();
    const cid = 7;
    const baseline = baselineComparison(cid);
    qc.setQueryData(queryKeys.comparison(cid), baseline);

    // 1) Own-op onMutate writes optimistic title "x-mine" into cache.
    qc.setQueryData(queryKeys.comparison(cid), { ...baseline, title: "x-mine" });

    // 2) Foreign frame arrives.
    const remote: SseEvent = {
      id: 99, kind: "comparison_submitted",
      entity_type: "comparison", entity_id: cid,
      actor: "bob", client_id: "tab2", client_op_id: "op-foreign",
      ts: "2026-05-14T10:00:00Z", payload: null,
      post_state: { ...baseline, title: "x-foreign", content_hash: "h1" },
    };
    applyRemoteToCache(remote, qc);

    // 3) Final state — foreign-wins. NOTE: the own-op write at step (1) is
    //    NOT driven through useQueueMutation here; this test only pins the
    //    sequential-write contract for the new whole-row cache key. Genuine
    //    rollback / re-run-onMutate coverage against this cache shape lives
    //    in replayCoordinator.test.ts.
    const after = qc.getQueryData<Comparison>(queryKeys.comparison(cid));
    expect(after?.title).toBe("x-foreign");
    expect(after?.content_hash).toBe("h1");
  });
});
