/**
 * Client-side `computeMemberSnapshot` (Plan §Phase 3, Task 3.3).
 *
 * Reads from the TanStack cache to derive the per-member snapshot at submit
 * time:
 * - `effective_peaks` from `queryKeys.peaks(exposureId)` (already mirrors
 *   the server's `auto − exclude ∪ add` union shape — see
 *   `routes_peaks.jl::GET /api/exposures/:id/peaks`). Excluded auto peaks
 *   are filtered out client-side.
 * - `confirmed_index` from `queryKeys.indices(exposureId)` filtered by
 *   R² ≥ `R2_GATE` and selected as the highest-scored member of the active
 *   custom group from `queryKeys.groups(exposureId)`. Mirrors the server-
 *   side `compute_member_snapshot` selection rule (score DESC NULLS LAST,
 *   id ASC).
 * - `analysis_inputs_hash` from the exposure entity cache.
 *
 * The MemberSnapshot type is defined in `api.ts` (single source of truth
 * for HTTP parser AND SSE handler). The submit-time output is sent verbatim
 * to the server — the dispatcher writes it to `comparison_members.snapshot`
 * without ever re-deriving from DB state.
 */
import type { QueryClient } from "@tanstack/react-query";
import type {
  Peak, IndexEntry, GroupEntry, Exposure,
  MemberSnapshot, MemberSnapshotPeak, MemberSnapshotConfirmedIndex,
} from "../../api";
import { queryKeys } from "../../queries";

/**
 * R²-gate for `confirmed_index`. Mirrors `CONFIRMED_INDEX_R2_GATE` in
 * `packages/HimalayaUI/src/comparisons.jl` and the UI hide-low-R² filter
 * in `PhasePanel`. Bumping must happen in both places in lockstep.
 */
export const CONFIRMED_INDEX_R2_GATE = 0.98;

export function computeMemberSnapshot(
  exposureId: number,
  qc: QueryClient,
): MemberSnapshot {
  const peaks = qc.getQueryData<Peak[]>(queryKeys.peaks(exposureId)) ?? [];
  const indices = qc.getQueryData<IndexEntry[]>(queryKeys.indices(exposureId)) ?? [];
  const groups = qc.getQueryData<GroupEntry[]>(queryKeys.groups(exposureId)) ?? [];
  const exposure = qc.getQueryData<Exposure>(queryKeys.exposure(exposureId));

  const effective_peaks: MemberSnapshotPeak[] = peaks
    .filter((p) => !(p.source === "auto" && p.excluded))
    // Sort by q to mirror server's `ORDER BY q` (test fixture pin).
    .slice()
    .sort((a, b) => a.q - b.q)
    .map((p) => ({
      id: p.id,
      q: p.q,
      // Manual peaks have intensity:null on the API; auto peaks carry it.
      intensity: p.source === "manual" ? null : p.intensity,
      // MemberSnapshotPeak.sharpness is `number` (not nullable). Coerce
      // null/undefined to 0 — defensive: cache rows can carry null
      // sharpness for legacy rows. The type contract holds; the worst
      // case is a member peak with `sharpness: 0`, which downstream
      // metadata rendering can detect (sharpness > 0 floor).
      sharpness: p.sharpness ?? 0,
      source: p.source,
    }));

  const confirmed_index = pickConfirmedIndex(indices, groups);
  const analysis_inputs_hash = exposure?.analysis_inputs_hash ?? "";

  return { effective_peaks, confirmed_index, analysis_inputs_hash };
}

/**
 * Pick the highest-scored R²-passing index from the active custom group.
 * Mirrors the SQL in `compute_member_snapshot`:
 *   ORDER BY i.score DESC NULLS LAST, i.id ASC LIMIT 1
 *
 * Returns null when:
 * - no `kind: "custom", active: true` group exists for the exposure
 * - no member of that group has `r_squared >= CONFIRMED_INDEX_R2_GATE`
 */
function pickConfirmedIndex(
  indices: IndexEntry[],
  groups: GroupEntry[],
): MemberSnapshotConfirmedIndex | null {
  const customGroup = groups.find((g) => g.kind === "custom" && g.active);
  if (!customGroup) return null;
  const memberSet = new Set<number>(customGroup.members);
  const candidates = indices.filter((ix) =>
    memberSet.has(ix.id)
    && ix.r_squared !== null
    && ix.r_squared >= CONFIRMED_INDEX_R2_GATE,
  );
  if (candidates.length === 0) return null;
  // Sort: score DESC NULLS LAST, then id ASC. JS sort is stable since ES2019.
  candidates.sort((a, b) => {
    // NULLS LAST: null/undefined → +Infinity sentinel for sort key, but we
    // want them to land *after* finite scores — flip via comparator.
    const aHasScore = typeof a.score === "number";
    const bHasScore = typeof b.score === "number";
    if (aHasScore && !bHasScore) return -1;
    if (!aHasScore && bHasScore) return 1;
    if (aHasScore && bHasScore) {
      // score DESC
      if (b.score! !== a.score!) return b.score! - a.score!;
    }
    // tie-break: id ASC
    return a.id - b.id;
  });
  const winner = candidates[0]!;
  // The IndexEntry.peaks list orders by ratio_position (server-side). The
  // server's `compute_member_snapshot` queries `index_peaks ORDER BY
  // ratio_position`. We mirror by reading `winner.peaks` (already in that
  // order from the API response) and projecting peak_id.
  const peak_ids = winner.peaks.map((p) => p.peak_id);
  return {
    id: winner.id,
    phase: winner.phase,
    // Spec types these as non-null `number`. The server-side compute helper
    // may produce null when the underlying `indices` row had a NULL
    // lattice_d; if that happens we still satisfy the type by coercing to 0
    // (and a confirmed index with lattice_d=0 surfaces as "unknown" in the
    // metadata gutter — same as a missing index).
    lattice_d: winner.lattice_d ?? 0,
    r_squared: winner.r_squared ?? 0,
    ngc: winner.ngc ?? 0,
    peak_ids,
  };
}
