import type { QueryClient, MutationCache } from "@tanstack/react-query";
import type { SseEvent } from "./types";
import { getDeferred, clearDeferred } from "./deferred";
import { applyRemoteToCache, applyPostStateOnly } from "./applyRemoteToCache";
import { getClientId } from "../clientId";

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
  const isOwnTab = remote.client_id !== undefined &&
                   remote.client_id !== null &&
                   remote.client_id === getClientId();

  // Case 1: SSE confirms a pending op of ours.
  if (remote.client_op_id) {
    const deferred = getDeferred(remote.client_op_id);
    if (deferred) {
      // Apply post_state (indices + exposure hash) BEFORE resolving the
      // deferred so that by the time the mutator's onSuccess runs the
      // indices cache is already fresh. Without this, the indices stay
      // frozen at the pre-mutation `inputs_hash` and StaleIndicesBanner
      // sticks until a hard refetch.
      applyPostStateOnly(remote, qc);
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

  // Self-echo guard (issue #8 from PR review). No deferred matches but the
  // frame originated in this same tab — i.e. HTTP-first won and the deferred
  // was already cleared in `mutationFn`'s finally. The mutator's onSuccess
  // already wrote the canonical peak row + exposure hash; falling into Case
  // 2 would double-apply the per-kind body (e.g. duplicate peak row for
  // peak_added). But we still need to propagate `post_state.indices` —
  // mutator onSuccess paths don't write the indices cache, so without this
  // the StaleIndicesBanner sticks.
  if (isOwnTab) {
    applyPostStateOnly(remote, qc);
    return;
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
  // the remote effect plus our pending pipeline on top. Capture each fresh
  // restore closure and swap it onto `m.state.context` so a later `onError`
  // rolls back to the post-foreign-event snapshot rather than the original
  // pre-rollback state (issue #3 from PR review). TanStack v5 sets
  // `state.context` once from the original `onMutate` return value and
  // doesn't replace it on its own — the assignment below overrides that.
  for (const m of pending) {
    const onMutate = m.options.onMutate as
      | ((vars: unknown) => unknown)
      | undefined;
    const fresh = onMutate?.(m.state.variables);
    if (fresh !== undefined) {
      // Cast through `unknown` because TanStack types `state.context` as
      // readonly at the public API surface; we're surgically updating it
      // here to keep the per-mutation rollback closure consistent.
      (m.state as unknown as { context: unknown }).context = fresh;
    }
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
 *
 * Kind-aware for `peak_added`: the SSE payload only carries `q` and
 * `peak_curation_id`, but `peakAddMutator.onSuccess` reads `id`, `source`,
 * and `exposure_id` off the response. Without this branch, the SSE-wins
 * path lands a peak row with `id === undefined` in the cache, which then
 * matches `hoveredPeakId` (also undefined) and renders a permanent halo.
 */
function synthesizeResponseFromSse(remote: SseEvent): unknown {
  const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
  const base = {
    event_id: remote.id,
    client_op_id: remote.client_op_id,
    analysis_inputs_hash: remote.post_state?.analysis_inputs_hash,
  };
  if (remote.kind === "peak_added") {
    const peakId = payload.peak_curation_id as number | undefined;
    return {
      ...base,
      id: peakId,
      exposure_id: remote.entity_id,
      q: payload.q as number,
      intensity: null,
      prominence: null,
      sharpness: null,
      source: "manual",
      excluded: false,
      view_row_id: peakId,
    };
  }
  if (remote.kind === "add_tag") {
    // SSE payload uses `tag_id`; HTTP response uses `id`. Mutators
    // (addSampleTagMutator / addExposureTagMutator) read response.id and
    // response.source. Without this branch, SSE-wins lands a malformed tag
    // with id=undefined and source=undefined — undeletable, mis-classified.
    return {
      ...base,
      id: payload.tag_id,
      key: payload.key,
      value: payload.value,
      source: "manual",
    };
  }
  return { ...base, ...payload };
}
