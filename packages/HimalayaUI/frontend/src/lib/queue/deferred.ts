import type { PendingDeferred } from "./types";

export const pendingDeferreds = new Map<string, PendingDeferred<unknown>>();
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
export function getDeferred(clientOpId: string): PendingDeferred<unknown> | undefined {
  return pendingDeferreds.get(clientOpId);
}
export function clearDeferred(clientOpId: string): void {
  pendingDeferreds.delete(clientOpId);
}