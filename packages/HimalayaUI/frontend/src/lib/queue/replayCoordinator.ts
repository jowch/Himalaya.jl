import type { QueryClient, MutationCache } from "@tanstack/react-query";
import type { SseEvent } from "./types";
import { getDeferred, clearDeferred } from "./deferred";
import { applyRemoteToCache } from "./applyRemoteToCache";

/**
 * Process an SSE frame against the local cache and pending mutation queue.
 *
 * Two cases:
 *   1. The frame's client_op_id matches a registered deferred: this is the
 *      server confirming our own pending op. Resolve the deferred with a
 *      synthesized response so the originating useQueueMutation can settle.
 *      Aborts any in-flight HTTP request via the deferred's controller — the
 *      HTTP response would be redundant with the SSE-synthesized one.
 *   2. No matching deferred: the event came from another tab or user.
 *      Rollback pending optimistic effects in reverse queue order, apply
 *      the remote effect to the cache, then re-run pending optimistic
 *      effects in insertion order. This is replay-as-rerun.
 *
 * MutationCache.getAll() preserves insertion order — a regression test
 * asserts this against TanStack version drift.
 *
 * Forward-scaffolded branches (post_message, set_exposure_status, add_tag,
 * remove_tag, update_sample, select_exposure) are unreachable today: their
 * routes still use log_action! and don't emit SSE. They will become live
 * once those routes migrate to apply_event! in M2.1. The default branch's
 * invalidate fallback would also handle them in the meantime.
 */
export function handleRemoteEvent(
  remote: SseEvent,
  qc: QueryClient,
  mc: MutationCache,
): void {
  // Case 1: SSE confirms a pending op of ours.
  if (remote.client_op_id) {
    const deferred = getDeferred(remote.client_op_id);
    if (deferred) {
      // Resolve first, THEN abort — Promises only settle once, so the
      // resolution wins and the abort-handler's reject becomes a no-op.
      // Order matters: aborting first triggers onAbort → reject, which
      // would settle the deferred before resolve() lands.
      deferred.resolve(synthesizeResponseFromSse(remote));
      // Abort the in-flight HTTP request — its response is redundant with
      // the SSE-synthesized one. abort() is a no-op if already aborted or
      // never wired up.
      deferred.controller?.abort();
      clearDeferred(remote.client_op_id);
      return;
    }
  }

  // Case 2: foreign event. Replay-as-rerun.
  const pending = mc.getAll().filter(
    (m) => m.state.status === "pending"
  );

  // Rollback in reverse queue order — last-in first-out so dependent ops
  // unwind before the ops they depend on.
  for (const m of [...pending].reverse()) {
    const ctx = m.state.context as { restore?: () => void } | undefined;
    ctx?.restore?.();
  }

  applyRemoteToCache(remote, qc);

  // Re-apply optimistic effects in insertion order so the cache reflects
  // the remote effect plus our pending pipeline on top.
  for (const m of pending) {
    const onMutate = m.options.onMutate as ((vars: unknown) => void) | undefined;
    onMutate?.(m.state.variables);
  }
}

/**
 * Build a synthetic response object mirroring what the HTTP route would
 * have returned. The deferred is awaited inside useQueueMutation's
 * mutationFn (M1.5); the resolution lets the mutation transition to
 * "success" without the HTTP call having returned. event_id is lifted
 * from the SSE frame; analysis_inputs_hash from post_state if present;
 * remaining fields are the route-specific payload (q, group_id, etc.)
 * that the original route handler would have echoed.
 */
function synthesizeResponseFromSse(remote: SseEvent): unknown {
  const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
  return {
    event_id: remote.id,
    client_op_id: remote.client_op_id,
    analysis_inputs_hash: remote.post_state?.analysis_inputs_hash,
    ...payload,
  };
}
