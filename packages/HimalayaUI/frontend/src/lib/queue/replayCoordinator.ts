import type { QueryClient, MutationCache } from "@tanstack/react-query";
import type { SseEvent } from "./types";
import type { Peak, GroupEntry, Exposure, Sample, SampleMessage } from "../../api";
import { getDeferred, clearDeferred } from "./deferred";
import { queryKeys } from "../../queries";

/**
 * Process an SSE frame against the local cache and pending mutation queue.
 *
 * Two cases:
 *   1. The frame's client_op_id matches a registered deferred: this is the
 *      server confirming our own pending op. Resolve the deferred with a
 *      synthesized response so the originating useQueueMutation can settle.
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
      deferred.resolve(synthesizeResponseFromSse(remote));
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

/**
 * Apply a remote SSE event to the local query cache. Per-kind logic mirrors
 * the spec's "replay-without-refetch where post_state covers it; refetch
 * fallback where the event payload is insufficient or update is rare."
 */
function applyRemoteToCache(remote: SseEvent, qc: QueryClient): void {
  const id = remote.entity_id;
  const payload = remote.payload as Record<string, unknown> | undefined;

  const applyPostState = (): void => {
    if (!remote.post_state) return;
    qc.setQueryData(queryKeys.indices(id), remote.post_state.indices);
    qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
      old
        ? { ...old, analysis_inputs_hash: remote.post_state!.analysis_inputs_hash }
        : old);
  };

  switch (remote.kind) {
    case "peak_added": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) => [
        ...old,
        // Negative id placeholder per the optimistic-id invariant in types.ts.
        // `kind: "add"` is the cache-shape used by tests/legacy code; the
        // canonical Peak shape uses `source: "manual"` — set both so existing
        // consumers and the test contract both pass.
        {
          id: -remote.id,
          exposure_id: id,
          q: payload?.q as number,
          intensity: null,
          prominence: null,
          sharpness: null,
          source: "manual",
          excluded: false,
          kind: "add",
        } as unknown as Peak,
      ]);
      applyPostState();
      break;
    }
    case "peak_excluded": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.map((p) =>
          Math.abs(p.q - (payload?.q as number)) < 1e-6
            ? { ...p, excluded: true }
            : p));
      applyPostState();
      break;
    }
    case "peak_unexcluded": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.map((p) =>
          Math.abs(p.q - (payload?.q as number)) < 1e-6
            ? { ...p, excluded: false }
            : p));
      applyPostState();
      break;
    }
    case "peak_removed": {
      qc.setQueryData<Peak[]>(queryKeys.peaks(id), (old = []) =>
        old.filter((p) => Math.abs(p.q - (payload?.q as number)) >= 1e-6));
      applyPostState();
      break;
    }
    case "index_confirmed": {
      const groupId = payload?.group_id as number;
      const indexId = payload?.index_id as number;
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(id), (old = []) =>
        old.map((g) =>
          g.id === groupId ? { ...g, members: [...g.members, indexId] } : g));
      break;
    }
    case "index_unconfirmed": {
      const groupId = payload?.group_id as number;
      const indexId = payload?.index_id as number;
      qc.setQueryData<GroupEntry[]>(queryKeys.groups(id), (old = []) =>
        old.map((g) =>
          g.id === groupId
            ? { ...g, members: g.members.filter((m) => m !== indexId) }
            : g));
      break;
    }
    case "speculative_created":
    case "speculative_deleted": {
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
      break;
    }
    case "set_exposure_status": {
      // Forward-scaffolded: route still uses log_action!.
      qc.setQueryData(queryKeys.exposure(id), (old: Exposure | undefined) =>
        old ? { ...old, status: payload?.status as Exposure["status"] } : old);
      break;
    }
    case "post_message": {
      // Forward-scaffolded: route still uses log_action!.
      const sampleId = payload?.sample_id as number;
      qc.setQueryData<SampleMessage[]>(queryKeys.messages(sampleId), (old = []) => [
        ...old,
        payload as unknown as SampleMessage,
      ]);
      break;
    }
    case "add_tag":
    case "remove_tag": {
      // Forward-scaffolded: route still uses log_action!.
      const parentKey = remote.entity_type === "sample"
        ? queryKeys.samples(payload?.experiment_id as number)
        : queryKeys.exposures(payload?.sample_id as number);
      qc.invalidateQueries({ queryKey: parentKey });
      break;
    }
    case "update_sample": {
      // Forward-scaffolded: route still uses log_action!.
      qc.setQueryData(queryKeys.sample(id), (old: Sample | undefined) =>
        old ? { ...old, ...(payload ?? {}) } : old);
      break;
    }
    case "select_exposure": {
      // Forward-scaffolded: route still uses log_action!.
      const sampleId = payload?.sample_id as number;
      qc.invalidateQueries({ queryKey: queryKeys.exposures(sampleId) });
      break;
    }
    case "analyze_run": {
      applyPostState();
      break;
    }
    default: {
      qc.invalidateQueries({ queryKey: queryKeys.peaks(id) });
      qc.invalidateQueries({ queryKey: queryKeys.indices(id) });
      qc.invalidateQueries({ queryKey: queryKeys.groups(id) });
    }
  }
}
