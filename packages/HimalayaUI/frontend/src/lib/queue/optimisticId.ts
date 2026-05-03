/**
 * Shared monotonic source for optimistic placeholder ids. All mutators that
 * mint negative-id placeholders must call this so that two concurrent
 * mutators in the same Date.now() ms can't collide on the same id.
 *
 * Encoding: `-(Date.now() * 1000 + counter)`. The 1000-multiplier reserves
 * three decimal digits for the per-ms counter — enough headroom for the
 * worst case (a burst of optimistic writes from one user gesture). Negative
 * sign keeps these distinct from server-assigned positive ids; consumers
 * should treat any `id < 0` as a placeholder until reconciliation.
 */
let optimisticCounter = 0;
export function nextOptimisticId(): number {
  optimisticCounter += 1;
  return -(Date.now() * 1000 + optimisticCounter);
}
