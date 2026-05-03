import type { PendingDeferred } from "./types";

/**
 * Process-global registry of deferred promises keyed by client_op_id. The
 * replay coordinator looks up by SSE frame's `client_op_id`; the queue
 * mutator framework registers on dispatch and clears on confirmation.
 *
 * Exported (vs. encapsulated behind getter/setter) so test code can clear
 * it between cases and observe registration directly. Production code
 * should use makeDeferred / getDeferred / clearDeferred — direct mutation
 * is reserved for tests and rehydration paths.
 */
export const pendingDeferreds = new Map<string, PendingDeferred<unknown>>();

/**
 * Mint a deferred and register it under `clientOpId`. The returned object
 * exposes the Promise plus its resolve/reject hooks so callers can resolve
 * out-of-band (e.g. from an SSE handler).
 *
 * Re-registering under an existing key overwrites silently — caller is
 * responsible for minting unique clientOpIds (newClientOpId() guarantees this
 * via crypto.randomUUID()).
 */
export function makeDeferred<T>(clientOpId: string): PendingDeferred<T> {
  let resolve!: (value: T) => void;
  let reject!: (reason: unknown) => void;
  const promise = new Promise<T>((res, rej) => {
    resolve = res;
    reject = rej;
  });
  const d: PendingDeferred<T> = { promise, resolve, reject };
  pendingDeferreds.set(clientOpId, d as PendingDeferred<unknown>);
  return d;
}

/** Look up a deferred by clientOpId. Returns undefined if not registered. */
export function getDeferred(clientOpId: string): PendingDeferred<unknown> | undefined {
  return pendingDeferreds.get(clientOpId);
}

/** Remove a deferred from the registry. Idempotent. */
export function clearDeferred(clientOpId: string): void {
  pendingDeferreds.delete(clientOpId);
}
