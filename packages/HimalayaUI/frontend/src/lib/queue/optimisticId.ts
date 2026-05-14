import { getClientId } from "../clientId";

/**
 * Shared monotonic source for optimistic placeholder ids. All mutators that
 * mint negative-id placeholders must call this so that two concurrent
 * mutators in the same tab can't collide on the same id, AND so that two
 * tabs of the same user can't collide either.
 *
 * Encoding: `-(tabFingerprint * 2^28 + counter)`.
 *   - `tabFingerprint`: 24-bit FNV-1a hash of the per-tab `clientId` from
 *     `lib/clientId.ts`. Stable for the lifetime of the tab; ~1-in-16M
 *     collision probability between any two tabs (effective space is
 *     2^24 − 1 — value 0 is remapped to 1 to keep `-counter`-only ids out
 *     of the encoding space).
 *   - `counter`: 28-bit monotonic per-tab counter (~268M ids per tab, far
 *     beyond any realistic session burst).
 *   - Total magnitude: 52 bits — within `Number.MAX_SAFE_INTEGER` (53 bits).
 *
 * The previous encoding used `Date.now() * 1000 + counter`, which collided
 * across tabs whose `Date.now()` ms happened to match at first call. This
 * encoding drops the wall-clock prefix because it never solved a real
 * problem — Date.now() was decorative, not load-bearing for uniqueness.
 *
 * Negative sign keeps placeholders distinct from server-assigned positive
 * ids; consumers should treat any `id < 0` as a placeholder until
 * reconciliation.
 */
let cachedTabFingerprint: number | null = null;
function tabFingerprint(): number {
  if (cachedTabFingerprint !== null) return cachedTabFingerprint;
  const id = getClientId();
  // FNV-1a 32-bit, folded to 24 bits.
  let h = 0x811c9dc5;
  for (let i = 0; i < id.length; i++) {
    h ^= id.charCodeAt(i);
    h = Math.imul(h, 0x01000193);
  }
  // Force into 24-bit unsigned range; remap a literal-zero fingerprint to 1
  // (zero would degenerate into bare -counter, indistinguishable across
  // tabs). `masked === 0 ? 1 : masked` preserves the full 24-bit cardinality
  // except for the single value 0 — an earlier version OR-ed with 0x1, which
  // halved the effective space by collapsing every even number onto its
  // odd-1 neighbour.
  const masked = (h >>> 8) & 0xffffff;
  cachedTabFingerprint = masked === 0 ? 1 : masked;
  return cachedTabFingerprint;
}

let optimisticCounter = 0;
export function nextOptimisticId(): number {
  optimisticCounter += 1;
  return -(tabFingerprint() * 0x10000000 + optimisticCounter);
}

// Test-only: reset the module-level counter and cached fingerprint so each
// test starts from a known state. Not part of the public surface.
export function __resetOptimisticIdForTest(): void {
  optimisticCounter = 0;
  cachedTabFingerprint = null;
}
