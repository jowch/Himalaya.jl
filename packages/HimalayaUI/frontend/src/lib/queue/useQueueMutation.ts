import { useMutation, useQueryClient } from "@tanstack/react-query";
import { newClientOpId } from "../clientOpId";
import { showToast } from "../toast";
import { makeDeferred, clearDeferred } from "./deferred";
import {
  isValidationError,
  isInfrastructureError,
  buildValidationMessage,
} from "./errors";
import type { FlatPayload, Mutator, RollbackContext } from "./types";

export interface UseQueueMutationResult<TInput> {
  mutate: (input: TInput) => void;
  isPending: boolean;
  isSuccess: boolean;
  error: unknown;
  reset: () => void;
}

const MAX_RETRIES = 5;
const MAX_BACKOFF_MS = 30_000;

/**
 * Wires a Mutator into TanStack Query's useMutation with:
 * - per-call client_op_id mint (newClientOpId() on every mutate())
 * - optimistic effect via mutator.onMutate, restored on error or replay-rollback
 * - deferred-promise resolution: HTTP and SSE race; first to land wins
 * - failure-class routing: 4xx → Validation toast, 5xx/network → infrastructure
 *   (banner is mounted by M2 at App.tsx and reads from useMutationState directly)
 * - retry policy: validation never retries; infrastructure retries up to 5x
 *   with exponential backoff capped at 30s
 *
 * `scope` is merged into every payload before mutator callbacks run. Use it
 * for invariants the mutator needs but that aren't part of the per-call input
 * (exposureId, sampleId, username, clientId).
 *
 * **Important:** the per-call clientOpId is minted INSIDE `mutate(input)`, not
 * at hook construction. Each mutate() call produces a unique idempotency key.
 */
export function useQueueMutation<TInput, TScope, TResponse>(
  mutator: Mutator<TInput, TScope, TResponse>,
  scope: TScope,
): UseQueueMutationResult<TInput> {
  const qc = useQueryClient();

  type Payload = FlatPayload<TInput, TScope>;

  const mutation = useMutation<TResponse, unknown, Payload, RollbackContext>(
    {
      mutationKey: [mutator.kind],
      mutationFn: async (payload) => {
        // TanStack v5 doesn't surface a per-call AbortSignal on the
        // MutationFunctionContext, so we mint our own. mutator.request still
        // honours signal.aborted; if the mutation is reset/cancelled future
        // wiring can call controller.abort().
        const controller = new AbortController();
        const { signal } = controller;
        const deferred = makeDeferred<TResponse>(payload.clientOpId);
        // Stash the controller on the deferred so the SSE-resolution path in
        // replayCoordinator can abort the HTTP request when SSE wins the race.
        deferred.controller = controller;
        // Wire HTTP into the deferred. If SSE arrives first, the HTTP
        // resolution becomes a no-op (deferred already cleared).
        mutator
          .request(payload, signal)
          .then((response) => deferred.resolve(response))
          .catch((err) => deferred.reject(err));
        // AbortSignal cleanup: clear the registry entry and reject deferred
        // so a leftover entry doesn't accumulate after early abort.
        const onAbort = () => {
          clearDeferred(payload.clientOpId);
          deferred.reject(new DOMException("aborted", "AbortError"));
        };
        signal.addEventListener("abort", onAbort);
        try {
          return await deferred.promise;
        } finally {
          signal.removeEventListener("abort", onAbort);
          clearDeferred(payload.clientOpId);
        }
      },
      onMutate: (payload) => mutator.onMutate(payload, qc),
      onSuccess: (response, payload) =>
        mutator.onSuccess(payload, response, qc),
      onError: (err, _payload, context) => {
        context?.restore?.();
        if (isValidationError(err)) {
          showToast(buildValidationMessage(mutator.kind, err), "error");
        }
        // Infrastructure errors: handled by retry; the banner reads from
        // useMutationState — nothing per-mutation to do here.
      },
      retry: (failureCount, err) => {
        if (isValidationError(err)) return false;
        if (isInfrastructureError(err)) return failureCount < MAX_RETRIES;
        return false;
      },
      retryDelay: (attempt) => Math.min(1000 * 2 ** attempt, MAX_BACKOFF_MS),
    },
  );

  const mutate = (input: TInput): void => {
    // The single cast at the framework layer: useMutation flat-spreads
    // {kind, clientOpId, ...scope, ...input} into the variables, but TS
    // can't statically prove that the resulting object satisfies the
    // intersection. Mutator callbacks receive this flat shape directly,
    // so consumer code never has to cast.
    const payload = {
      kind: mutator.kind,
      clientOpId: newClientOpId(),
      ...scope,
      ...(input as object),
    } as unknown as Payload;
    mutation.mutate(payload);
  };

  return {
    mutate,
    isPending: mutation.isPending,
    isSuccess: mutation.isSuccess,
    error: mutation.error,
    reset: mutation.reset,
  };
}
