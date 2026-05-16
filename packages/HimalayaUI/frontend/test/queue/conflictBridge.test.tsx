/**
 * conflictBridge — single MutationCache subscriber that bridges
 * ConflictError on `comparison_save` mutations to the Zustand
 * `pendingConflict` slot (Phase 12 follow-up; queue-reviewer Fix 1+2).
 *
 * The previous bridge lived inside `useSaveComparison`'s `useEffect`. With
 * two mount sites (`ComparePageEdit` + `ConflictModal`), the same error
 * routed through TWO bridges that both wrote to ONE Zustand slot — a race
 * that flipped `current_state` based on render ordering. This subscriber
 * lifts the bridge to a single source of truth: one mount at App startup,
 * one writer to the slot, idempotent on remount/re-subscribe.
 *
 * Asserts:
 * - Single ConflictError on comparison_save mutation → setPendingConflict
 *   fires exactly once with the error.
 * - Foreign mutations (different mutationKey) are ignored.
 * - Non-ConflictError errors are ignored.
 * - Resubscribe after detach does NOT re-fire on a stale error already
 *   bridged on a previous attach (last-seen mutationId tracking).
 * - A second mutation that lands a fresh ConflictError fires the bridge
 *   again with the new error (re-callability).
 */
import { describe, it, expect, beforeEach } from "vitest";
import { QueryClient, MutationCache, MutationObserver } from "@tanstack/react-query";
import { ConflictError } from "../../src/api";
import type { Comparison } from "../../src/api";
import { attachConflictBridge, _resetConflictBridgeForTest } from "../../src/lib/queue/conflictBridge";

function buildComparison(id: number, hash = "sha256:server"): Comparison {
  return {
    id, title: "X", description: null, content_hash: hash,
    created_by: 1, created_at: "2026-05-06T00:00:00Z",
    updated_at: "2026-05-06T00:00:00Z",
    forked_from_id: null, forked_at_hash: null, forked_from_title: null,
    view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: null,
    members: [],
  };
}

/**
 * Run a mutation that throws `err` to completion against the given client,
 * with the given mutationKey. Returns when the mutation has settled.
 */
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
  // Subscribe and wait for terminal state.
  await new Promise<void>((resolve) => {
    const unsub = observer.subscribe((res) => {
      if (res.status === "error" || res.status === "success") {
        unsub();
        resolve();
      }
    });
    // Swallow the rejected promise — we only care about the cache state.
    observer.mutate(undefined as never).catch(() => {});
  });
}

describe("attachConflictBridge", () => {
  let qc: QueryClient;
  let mc: MutationCache;
  let calls: Array<ConflictError | null>;
  const setPendingConflict = (c: ConflictError | null) => { calls.push(c); };

  beforeEach(() => {
    _resetConflictBridgeForTest();
    qc = new QueryClient({
      defaultOptions: { mutations: { retry: false } },
    });
    mc = qc.getMutationCache();
    calls = [];
  });

  it("fires setPendingConflict exactly once for a ConflictError on comparison_save", async () => {
    const detach = attachConflictBridge(mc, setPendingConflict);
    const err = new ConflictError("sha256:server", buildComparison(42));
    await runFailingMutation(qc, ["comparison_save"], err);
    // Filter to non-null calls (the bridge does not write null on errors).
    const errorCalls = calls.filter((c) => c instanceof ConflictError);
    expect(errorCalls).toHaveLength(1);
    expect(errorCalls[0]).toBe(err);
    detach();
  });

  it("ignores foreign mutations (different mutationKey)", async () => {
    const detach = attachConflictBridge(mc, setPendingConflict);
    const err = new ConflictError("sha256:server", buildComparison(42));
    // Wrong key — bridge must NOT fire even though the error type is right.
    await runFailingMutation(qc, ["peak_add"], err);
    expect(calls).toEqual([]);
    detach();
  });

  it("ignores non-ConflictError errors on comparison_save", async () => {
    const detach = attachConflictBridge(mc, setPendingConflict);
    await runFailingMutation(qc, ["comparison_save"], new Error("plain"));
    expect(calls).toEqual([]);
    detach();
  });

  it("two simultaneous subscribers don't double-write the slot", async () => {
    // Simulates Fix 1's race: previously two `useSaveComparison()` mount sites
    // each wired their own bridge effect; both fired on the same mutation.
    // With the lifted subscriber there's a single subscription, but we go
    // further and pin: even if `attachConflictBridge` is called twice, the
    // shared module-state guard keeps the slot writes deduped. (Defense in
    // depth — App.tsx still mounts it once.)
    const detach1 = attachConflictBridge(mc, setPendingConflict);
    const detach2 = attachConflictBridge(mc, setPendingConflict);
    const err = new ConflictError("sha256:server", buildComparison(42));
    await runFailingMutation(qc, ["comparison_save"], err);
    const errorCalls = calls.filter((c) => c instanceof ConflictError);
    // Two subscribers may each see the error event, but the shared
    // last-seen-mutationId guard collapses them to one slot write.
    expect(errorCalls).toHaveLength(1);
    detach1();
    detach2();
  });

  it("resubscribe after detach does not re-pop the modal on a stale error (Fix 2)", async () => {
    // First attach + mutation: bridges the error.
    let detach = attachConflictBridge(mc, setPendingConflict);
    const err = new ConflictError("sha256:server", buildComparison(42));
    await runFailingMutation(qc, ["comparison_save"], err);
    expect(calls.filter((c) => c instanceof ConflictError)).toHaveLength(1);
    // Detach (e.g. App unmounts under React StrictMode or HMR).
    detach();
    // Reset call log for the second attach; the mutation in the cache is
    // still there with its terminal error state.
    calls.length = 0;
    // Reattach. The subscriber must NOT re-fire on the stale error.
    detach = attachConflictBridge(mc, setPendingConflict);
    // Spin a microtask to flush any synchronous initial-scan emissions.
    await Promise.resolve();
    expect(calls).toEqual([]);
    detach();
  });

  it("a NEW conflict (different mutationId) re-fires the bridge", async () => {
    const detach = attachConflictBridge(mc, setPendingConflict);
    const err1 = new ConflictError("sha256:v1", buildComparison(42, "sha256:v1"));
    await runFailingMutation(qc, ["comparison_save"], err1);
    const err2 = new ConflictError("sha256:v2", buildComparison(42, "sha256:v2"));
    await runFailingMutation(qc, ["comparison_save"], err2);
    const errorCalls = calls.filter((c) => c instanceof ConflictError);
    // Each mutation has its own mutationId — both fire the bridge.
    expect(errorCalls).toHaveLength(2);
    expect(errorCalls[0]).toBe(err1);
    expect(errorCalls[1]).toBe(err2);
    detach();
  });
});
