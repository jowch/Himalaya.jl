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
    normalization: (saved.normalization as DraftMemberNormalization) ?? "max",
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
    normalization: "max",
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
    forkedFromId: undefined,
    forkedAtHash: undefined,
    viewGroupingMode:  (c.view_grouping_mode  as ActiveDraft["viewGroupingMode"])  ?? undefined,
    viewShowPeakTicks:  c.view_show_peak_ticks  ?? undefined,
    viewShowPeakLabels: c.view_show_peak_labels ?? undefined,
  };
}

/**
 * Fork-draft factory (Plan §Phase 11, Task 11.2). Creates a brand-new
 * draft (no `id`, no `baseHash`) pre-populated from a parent comparison's
 * members and lineage. The submit will route to `POST /api/comparisons`
 * (create), with `forked_from_id` + `forked_at_hash` riding the body.
 *
 * Member id is deliberately dropped — each member becomes a fresh INSERT
 * on the new comparison. Snapshot recomputes against the current cache so
 * the fork starts at "current truth" rather than the parent's frozen
 * snapshot (matches `loadDraftFromComparison` recovery semantics).
 *
 * Title defaults to "Fork of <parent title>" so the user has something
 * meaningful in the title field; they can rename before submit.
 */
export function fromComparisonAsFork(c: Comparison, qc: QueryClient): ActiveDraft {
  return {
    id: undefined,
    baseHash: undefined,
    title: `Fork of ${c.title ?? "comparison"}`,
    description: c.description ?? "",
    members: c.members.map((m) => {
      const dm = memberFromSaved(m, qc);
      // Drop server id — each member becomes an INSERT on the new
      // comparison. memberFromSaved already recomputed the snapshot.
      return { ...dm, id: undefined };
    }),
    forkedFromId: c.id,
    forkedAtHash: c.content_hash,
    // Carry view choices from the parent so the fork opens with the same
    // visual defaults; the author can change them before saving.
    viewGroupingMode:  (c.view_grouping_mode  as ActiveDraft["viewGroupingMode"])  ?? undefined,
    viewShowPeakTicks:  c.view_show_peak_ticks  ?? undefined,
    viewShowPeakLabels: c.view_show_peak_labels ?? undefined,
  };
}
