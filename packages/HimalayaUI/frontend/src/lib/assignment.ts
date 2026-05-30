import type { Assignment, AssignmentState, IndexEntry } from "../api";

/**
 * The active phase set, derived from the durable assignment cart. Replaces the
 * legacy `groups.find(g => g.active)?.members` read.
 *
 * Only an `indexed` assignment has active indices; `form_factor` and `null`
 * declarations always yield [] (state is explicit, never inferred from
 * members.length — see Plan D §"The 3-state model on the wire"). Member ids
 * not present in `indices` (e.g. a just-deleted candidate) are dropped.
 */
export function deriveActiveIndices(
  a: Assignment | undefined,
  indices: IndexEntry[],
): IndexEntry[] {
  if (!a || a.state !== "indexed") return [];
  const ids = new Set(a.members);
  return indices.filter((ix) => ids.has(ix.id));
}

/**
 * The assignment's 3-state, defaulting to "indexed" when no assignment row
 * exists yet (mirrors the backend's neutral default in _assignment_body).
 */
export function assignmentState(a: Assignment | undefined): AssignmentState {
  return a?.state ?? "indexed";
}
