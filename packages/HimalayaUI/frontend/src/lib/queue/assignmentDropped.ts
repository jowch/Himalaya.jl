/**
 * F-WIPE W2 — user-visible announcement when a reanalysis drops assignment
 * members from the phase call.
 *
 * Backend context (W1): reanalysis re-attaches assignment members
 * semantically; a member only drops when its index genuinely no longer
 * exists in the new candidate set (or merged into a sibling). The peak_*
 * SSE frame then carries `post_state.assignment_dropped: string[]` — one
 * phase name PER MEMBER, possibly repeating. This module aggregates that
 * list into one message and surfaces it.
 *
 * Surface choice: a visible toast (kind "warning") — losing part of the
 * phase call is a consequential state change, same tier as the queue's
 * validation-failure toasts. No separate `announce()` call: the toast
 * surface itself is an aria-live region (Toast.tsx renders warning toasts
 * with role="alert" / aria-live="assertive" and the severity word), so a
 * parallel LiveRegion announcement would double-speak to screen readers.
 */
import { showToast } from "../toast";

/**
 * Aggregate the per-member dropped-phase list into one or two short
 * sentences. Counts repeats per phase (first-seen order):
 *
 *   ["Pn3m"]                     → "Pn3m index dropped from the call. …"
 *   ["Pn3m", "Pn3m"]             → "2 Pn3m indices dropped from the call. …"
 *   ["Pn3m", "Lamellar"]         → "Pn3m and Lamellar indices dropped …"
 *   ["Pn3m", "Im3m", "Lamellar"] → "Pn3m, Im3m and Lamellar indices …"
 */
export function buildAssignmentDroppedMessage(dropped: string[]): string {
  const counts = new Map<string, number>();
  for (const phase of dropped) counts.set(phase, (counts.get(phase) ?? 0) + 1);
  const groups = [...counts.entries()].map(
    ([phase, n]) => (n > 1 ? `${n} ${phase}` : phase),
  );
  const head = groups.length === 1
    ? groups[0]!
    : `${groups.slice(0, -1).join(", ")} and ${groups[groups.length - 1]!}`;
  const total = dropped.length;
  const noun = total === 1 ? "index" : "indices";
  const tail = total === 1
    ? "Its index no longer fits the peaks."
    : "They no longer fit the peaks.";
  return `${head} ${noun} dropped from the call. ${tail}`;
}

/**
 * Fire the dropped-members toast. No-op on absent/empty input so callers
 * can pass `post_state.assignment_dropped` straight through.
 */
export function announceAssignmentDropped(dropped: string[] | undefined): void {
  if (!dropped || dropped.length === 0) return;
  showToast(buildAssignmentDroppedMessage(dropped), "warning");
}
