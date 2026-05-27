/**
 * conflictBridge — series-commit widening (I3.5b).
 *
 * The bridge accepts BOTH `comparison_save` and `series_commit` mutationKeys.
 * Series recipe-save (`series_save`) never throws ConflictError (PATCH never
 * 409s), so a series_save error must NOT reach the slot. The comparison_save
 * path is regression-guarded here too.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient, MutationCache, MutationObserver } from "@tanstack/react-query";
import { ConflictError } from "../../src/api";
import type { Series, Comparison } from "../../src/api";
import { attachConflictBridge, _resetConflictBridgeForTest } from "../../src/lib/queue/conflictBridge";

function buildSeries(id: number, hash = "sha256:server"): Series {
  return {
    id, title: "S", description: null, content_hash: hash,
    created_by: 1, created_at: null, updated_at: null,
    forked_from_id: null, forked_at_hash: null, forked_from_title: null,
    view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, ordering_variable: null, order_rule: "manual",
    state: "committed", members: [], samples: [],
  };
}

function buildComparison(id: number): Comparison {
  return {
    id, title: "X", description: null, content_hash: "sha256:c",
    created_by: 1, created_at: null, updated_at: null,
    forked_from_id: null, forked_at_hash: null, forked_from_title: null,
    view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: null, members: [],
  };
}

async function runFailingMutation(
  qc: QueryClient,
  mutationKey: ReadonlyArray<unknown>,
  err: unknown,
): Promise<void> {
  const observer = new MutationObserver(qc, {
    mutationKey: mutationKey as unknown[],
    mutationFn: async () => { throw err; },
    retry: false,
  });
  await new Promise<void>((resolve) => {
    const unsub = observer.subscribe((res) => {
      if (res.status === "error" || res.status === "success") {
        unsub();
        resolve();
      }
    });
    observer.mutate(undefined as never).catch(() => {});
  });
}

describe("attachConflictBridge — series_commit", () => {
  let qc: QueryClient;
  let mc: MutationCache;
  let calls: Array<ConflictError | null>;
  const setPendingConflict = (c: ConflictError | null) => { calls.push(c); };

  beforeEach(() => {
    _resetConflictBridgeForTest();
    qc = new QueryClient({ defaultOptions: { mutations: { retry: false } } });
    mc = qc.getMutationCache();
    calls = [];
  });

  it("bridges a ConflictError on series_commit to the slot", async () => {
    const detach = attachConflictBridge(mc, setPendingConflict);
    const err = new ConflictError("sha256:server", buildSeries(5));
    await runFailingMutation(qc, ["series_commit"], err);
    const errorCalls = calls.filter((c) => c instanceof ConflictError);
    expect(errorCalls).toHaveLength(1);
    expect(errorCalls[0]).toBe(err);
    detach();
  });

  it("does NOT bridge a series_save error (recipe-save never conflicts)", async () => {
    const detach = attachConflictBridge(mc, setPendingConflict);
    // Even if a ConflictError somehow rode a series_save mutation, the bridge
    // must ignore the kind — series_save is not a conflict-bearing key.
    const err = new ConflictError("sha256:server", buildSeries(5));
    await runFailingMutation(qc, ["series_save"], err);
    expect(calls).toEqual([]);
    detach();
  });

  it("does NOT bridge comparison_save (Compare retired, #177)", async () => {
    // I3.6: the Compare page is gone, so the `comparison_save` arm was removed
    // from the bridge. A ConflictError on a (replay-only) comparison_save
    // mutation must NOT reach the slot — only `series_commit` does.
    const detach = attachConflictBridge(mc, setPendingConflict);
    const err = new ConflictError("sha256:server", buildComparison(42));
    await runFailingMutation(qc, ["comparison_save"], err);
    expect(calls.filter((c) => c instanceof ConflictError)).toHaveLength(0);
    detach();
  });
});
