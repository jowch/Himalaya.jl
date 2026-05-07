import type { MutationCache } from "@tanstack/react-query";
import { ConflictError } from "../../api";

/**
 * Single-source-of-truth bridge from `comparison_save` mutation errors to
 * the Zustand `pendingConflict` slot (Phase 12 follow-up; queue-reviewer
 * Fix 1 + Fix 2).
 *
 * Why a MutationCache subscriber, not a hook
 * ------------------------------------------
 * `useSaveComparison()` is mounted at TWO sites (`ComparePageEdit` for the
 * edit-mode Save button; `ConflictModal` for the Overwrite re-submit).
 * Each call instantiates its own `useMutation` registration with its own
 * `result.error`. If the bridge lives inside the hook's `useEffect`, both
 * effects race to write the same Zustand slot — the second-409 race could
 * flip-flop `current_state` based on render ordering.
 *
 * Lifting the bridge to a single MutationCache subscriber gives us:
 *   - one writer to `pendingConflict`, regardless of how many components
 *     call `useSaveComparison`,
 *   - last-seen `mutationId` tracking so a remount/re-subscribe (StrictMode
 *     double-invocation, App-level remount, HMR) does NOT re-pop the modal
 *     on a stale terminal-error mutation that's still in the cache.
 *
 * The dedup state lives in module scope (not in a closure on each `attach`
 * call) so multiple `attach` calls share a single guard, and detach +
 * re-attach preserves the "already bridged" set.
 *
 * The MutationCache is the same one used by `replayCoordinator` and
 * `attachPersistence` — its `subscribe` callback fires on every cache event
 * (added, updated, removed, observer*). We filter by `type === "updated"` +
 * terminal `error` state + `mutationKey === ["comparison_save"]` +
 * `error instanceof ConflictError`.
 *
 * Usage: mount once at App startup
 * --------------------------------
 *   useEffect(() => attachConflictBridge(mc, setPendingConflict), [mc]);
 */

const COMPARISON_SAVE_KIND = "comparison_save";

/** mutationIds we have already bridged to the slot. Module-scoped so detach
 *  + re-attach (StrictMode, HMR) does not re-fire on a still-cached error. */
const bridged: Set<number> = new Set();

/**
 * Test-only reset. Production code should never call this — clearing the
 * `bridged` set in flight would re-pop the modal on the next subscribe.
 */
export function _resetConflictBridgeForTest(): void {
  bridged.clear();
}

/**
 * Subscribe to the MutationCache, bridging ConflictError on
 * `comparison_save` mutations to `setPendingConflict`. Returns an
 * unsubscribe function.
 *
 * Idempotent re-attach: the module-scoped `bridged` set persists across
 * detach/re-attach so a stale error in the cache does not re-pop the modal.
 */
export function attachConflictBridge(
  mc: MutationCache,
  setPendingConflict: (err: ConflictError | null) => void,
): () => void {
  return mc.subscribe((event) => {
    // We only care about state transitions that yield a terminal error.
    // The subscribe callback fires for many event types; pin the one
    // carrying error-state changes.
    if (event.type !== "updated") return;
    const mutation = event.mutation;
    if (mutation.state.status !== "error") return;
    const key = mutation.options.mutationKey;
    // mutationKey is `[mutator.kind]` — see useQueueMutation. The
    // saveComparison mutator's kind is "comparison_save".
    if (!Array.isArray(key) || key[0] !== COMPARISON_SAVE_KIND) return;
    const err = mutation.state.error;
    if (!(err instanceof ConflictError)) return;
    // Dedupe: each `Mutation` instance has a stable `mutationId`. We bridge
    // exactly once per mutation — a re-subscribe sees the same id and skips.
    // A NEW mutation (e.g. the user clicks Overwrite, second 409 arrives)
    // gets its own id and re-fires the bridge.
    if (bridged.has(mutation.mutationId)) return;
    bridged.add(mutation.mutationId);
    setPendingConflict(err);
  });
}
