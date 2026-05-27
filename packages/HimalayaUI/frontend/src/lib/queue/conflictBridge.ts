import type { MutationCache } from "@tanstack/react-query";
import { ConflictError } from "../../api";

/**
 * Single-source-of-truth bridge from conflict-bearing mutation errors
 * (`comparison_save` submit, `series_commit` plate-commit) to the Zustand
 * `pendingConflict` slot (Phase 12 follow-up; queue-reviewer Fix 1 + Fix 2;
 * series-commit widening I3.5b).
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
// I3.5b — series plate-commit is the only series mutation that can 409
// (recipe-save `series_save` never reads `expected_content_hash`). The bridge
// accepts both conflict-bearing kinds; `series_save` is deliberately excluded.
const SERIES_COMMIT_KIND = "series_commit";
const CONFLICT_KINDS: ReadonlySet<string> = new Set([
  COMPARISON_SAVE_KIND,
  SERIES_COMMIT_KIND,
]);

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
    // mutationKey is `[mutator.kind]` — see useQueueMutation. Conflict-bearing
    // kinds: "comparison_save" (submit) and "series_commit" (plate commit).
    // "series_save" (recipe-save PATCH) is excluded — it never 409s.
    if (!Array.isArray(key) || typeof key[0] !== "string" || !CONFLICT_KINDS.has(key[0])) return;
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
