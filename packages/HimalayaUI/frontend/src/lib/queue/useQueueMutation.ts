import { useMutation, useQueryClient } from "@tanstack/react-query";
import { newClientOpId } from "../clientOpId";
import { showToast } from "../toast";
import { ConflictError } from "../../api";
import { makeDeferred, clearDeferred } from "./deferred";
import {
  isValidationError,
  isInfrastructureError,
  is404Error,
  buildValidationMessage,
} from "./errors";
import type { FlatPayload, Mutator, RollbackContext } from "./types";

/**
 * Per-call lifecycle callbacks for `mutate`. They fire AFTER the queue's own
 * global handlers (optimistic rollback, validation toast, retry policy):
 * - `onSuccess` fires on CONFIRMATION — HTTP response or own-op SSE frame,
 *   whichever wins the race. Consumers use it for honest success toasts
 *   (never toast optimistically at mutate() time).
 * - `onError` fires on TERMINAL failure — for infrastructure errors that is
 *   after the retry policy exhausts (the InfrastructureBanner only covers the
 *   pending window; once the mutation settles as error it disappears, so a
 *   consumer that claimed nothing would leave a silent rollback).
 * TanStack semantics apply: when mutate() is called again before the previous
 * call settles, only the LATEST call's callbacks fire (no toast pile-up).
  * NOTE: onError also fires for 404s a `treats404AsSuccess` mutator treats
 * as success (the mutation still settles as error; the global handler merely
 * skips the rollback) - guard with isValidationError/is404Error if your
 * mutator opts in.
 */
export interface MutateCallbacks<TResponse = unknown> {
  onSuccess?: (response: TResponse) => void;
  onError?: (error: unknown) => void;
}

export interface UseQueueMutationResult<TInput, TResponse = unknown> {
  mutate: (input: TInput, callbacks?: MutateCallbacks<TResponse>) => void;
  isPending: boolean;
  isSuccess: boolean;
  /** Last successful response, or undefined if no mutation has succeeded yet. */
  data: TResponse | undefined;
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
): UseQueueMutationResult<TInput, TResponse> {
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
        // Mutators that opt into `treats404AsSuccess` (idempotent removes)
        // swallow a 404 entirely: the optimistic effect already reflects the
        // desired end state, and re-running rollback would re-insert a row
        // the server has already deleted. Skip both restore and toast.
        if (mutator.treats404AsSuccess && is404Error(err)) return;
        context?.restore?.();
        if (isValidationError(err)) {
          // ConflictError (409, content_hash drift) is caught here solely to
          // suppress the generic error toast — there is no conflict modal
          // (multiplayer is last-write-wins; conflict UI was cancelled).
          // The same pattern applies to any future typed-throw that has its
          // own bespoke error surface; gate via instanceof.
          if (err instanceof ConflictError) return;
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

  const mutate = (input: TInput, callbacks?: MutateCallbacks<TResponse>): void => {
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
    // Conditional spread (exactOptionalPropertyTypes): never pass an explicit
    // `undefined` handler into TanStack's MutateOptions.
    mutation.mutate(payload, {
      ...(callbacks?.onSuccess
        ? { onSuccess: (response: TResponse) => callbacks.onSuccess!(response) }
        : {}),
      ...(callbacks?.onError
        ? { onError: (err: unknown) => callbacks.onError!(err) }
        : {}),
    });
  };

  return {
    mutate,
    isPending: mutation.isPending,
    isSuccess: mutation.isSuccess,
    data: mutation.data ?? undefined,
    error: mutation.error,
    reset: mutation.reset,
  };
}
