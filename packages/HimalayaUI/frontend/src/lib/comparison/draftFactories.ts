/**
 * Snapshot-deriving factories for the Compare draft slot
 * (Plan §Phase 4, Task 4.3).
 *
 * Lives in its own file so `state.ts` (which `draft.ts` imports) doesn't
 * pull `queries.ts` (and through it the queue/mutator graph) at module
 * init time. `state.ts` calls these factories *inside* action bodies,
 * which evaluate after all modules have finished loading.
 *
 * Recovery rule (load-bearing): `memberFromSaved` recomputes the snapshot
 * against the **current** TanStack cache, not the saved server snapshot.
 * Resubmitting an existing comparison therefore "refreshes to current
 * truth" — see Plan §Task 4.3 and the spec's stale-comparison flow.
 */
import type { QueryClient } from "@tanstack/react-query";
import type { Comparison, ComparisonMember, MemberSnapshot } from "../../api";
import { computeMemberSnapshot } from "./snapshot";
import type { ActiveDraft, DraftMember, DraftMemberNormalization } from "./draft";

export function memberFromSaved(
  saved: ComparisonMember,
  qc: QueryClient,
): DraftMember {
  const fresh: MemberSnapshot | undefined =
    saved.exposure_id !== null
      ? computeMemberSnapshot(saved.exposure_id, qc)
      : (saved.snapshot ?? undefined);
  return {
    id: saved.id,
    exposure_id: saved.exposure_id,
    display_order: saved.display_order,
    band_height: saved.band_height,
    y_offset: saved.y_offset,
    normalization: (saved.normalization as DraftMemberNormalization) ?? "none",
    color_override: saved.color_override ?? undefined,
    label_override: saved.label_override ?? undefined,
    q_window_min: saved.q_window_min ?? undefined,
    q_window_max: saved.q_window_max ?? undefined,
    peak_display:
      saved.peak_display
      && typeof saved.peak_display === "object"
      && "hidden" in (saved.peak_display as object)
      && "labeled" in (saved.peak_display as object)
        ? (saved.peak_display as { hidden: number[]; labeled: number[] })
        : undefined,
    snapshot: fresh,
  };
}

export function memberFromNewExposure(
  exposureId: number,
  displayOrder: number,
  qc: QueryClient,
): DraftMember {
  return {
    id: undefined,
    exposure_id: exposureId,
    display_order: displayOrder,
    band_height: 1.0,
    y_offset: 0.0,
    normalization: "none",
    color_override: undefined,
    label_override: undefined,
    q_window_min: undefined,
    q_window_max: undefined,
    peak_display: undefined,
    snapshot: computeMemberSnapshot(exposureId, qc),
  };
}

export function fromComparison(c: Comparison, qc: QueryClient): ActiveDraft {
  return {
    id: c.id,
    baseHash: c.content_hash,
    title: c.title ?? "",
    description: c.description ?? "",
    members: c.members.map((m) => memberFromSaved(m, qc)),
  };
}
