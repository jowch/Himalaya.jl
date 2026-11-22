import type { QueryClient, MutationCache } from "@tanstack/react-query";
import type { SseEvent, QueueResponseMeta, CurationPostState } from "./types";
import { getDeferred, clearDeferred } from "./deferred";
import { applyRemoteToCache, applyPostStateOnly } from "./applyRemoteToCache";
import { getClientId } from "../clientId";
import { resolveMutatorForEvent } from "./mutatorRegistry";

export function handleRemoteEvent(
  remote: SseEvent,
  qc: QueryClient,
  mc: MutationCache,
): void {
  const isOwnTab = remote.client_id !== undefined &&
                   remote.client_id !== null &&
                   remote.client_id === getClientId();
  if (remote.client_op_id) {
    const deferred = getDeferred(remote.client_op_id);
    if (deferred) {
      applyPostStateOnly(remote, qc);
      deferred.resolve(synthesizeResponseFromSse(remote));
      deferred.controller?.abort();
      clearDeferred(remote.client_op_id);
      return;
    }
  }
  if (isOwnTab) {
    applyPostStateOnly(remote, qc);
    return;
  }
  const pending = mc.getAll().filter(
    (m) => m.state.status === "pending"
  );
  for (const m of [...pending].reverse()) {
    const ctx = m.state.context as { restore?: () => void } | undefined;
    ctx?.restore?.();
  }
  applyRemoteToCache(remote, qc);
  for (const m of pending) {
    const onMutate = m.options.onMutate as
      | ((vars: unknown) => unknown)
      | undefined;
    try {
      const fresh = onMutate?.(m.state.variables);
      if (fresh !== undefined) {
        (m.state as unknown as { context: unknown }).context = fresh;
      }
    } catch (e) {
      console.error("[handleRemoteEvent] onMutate threw during re-run:", e);
    }
  }
}
function synthesizeResponseFromSse(remote: SseEvent): unknown {
  const base: QueueResponseMeta = {
    event_id: remote.id,
    client_op_id: remote.client_op_id,
    analysis_inputs_hash: (remote.post_state as CurationPostState | undefined)
      ?.analysis_inputs_hash,
  };
  const mutator = resolveMutatorForEvent(remote.kind);
  const synth = mutator?.synthesizeFromSse?.(remote, base);
  if (synth !== undefined) return synth;
  const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
  return { ...base, ...payload };
}