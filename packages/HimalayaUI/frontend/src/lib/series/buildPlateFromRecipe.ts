/**
 * Pure recipe→plate resolver for the builder's Confirm chain (P0
 * BU-RECIPENOOP). `PATCH /api/series/:id` persists the RECIPE
 * (`series_samples`) but does NOT rebuild the plate (`series_members`), and
 * `POST /api/series/:id/commit` takes its members verbatim from the request
 * body — so plate resolution is the FRONTEND's job. Committing the cached
 * `series.members` re-posts the old plate byte-for-byte: reorders never
 * surface and an added sample never becomes a member.
 *
 * Resolution mirrors the backend's create-path resolver
 * (`_resolve_series_plate!`, events.jl) field for field:
 *   - order = recipe `position` (tiebreak recipe-row id), excluded rows skipped
 *     (a deliberate exclusion, not a failure — same `excluded = 0` filter);
 *   - each sample resolves to its indexing/representative exposure (the
 *     picker's `indexing_exposure_id`, which encodes the SAME highest-id
 *     selected-else-highest-id rule — comparisons.jl
 *     `_picker_samples_projection`, mirrored by events.jl
 *     `_resolve_series_plate!` on the create path);
 *   - `display_order` is sequential over the RESOLVED members (no gaps);
 *   - a sample with no resolvable exposure is left out and REPORTED in
 *     `unresolvedSampleIds` — the caller owes the user an honest toast
 *     (silently dropping it would lie about what was published).
 *
 * Display props (band_height, y_offset, normalization, overrides, q-window,
 * peak_display, snapshot) CARRY OVER from an old member with the same
 * exposure_id (consume-once, so a duplicated exposure can't double-spend one
 * old member). Brand-new members send only `{exposure_id, display_order}` —
 * the commit route's `_series_member_payload` fills the same defaults the
 * create path uses (band_height 1.0, y_offset 0, normalization "none", null
 * overrides/windows) and computes the snapshot server-side.
 */
import type { SeriesMember, SeriesMemberInput, SeriesSample } from "../../api";
import { memberToInput } from "./buildSeriesCommitBody";

export interface PlateFromRecipe {
  /** The resolved plate, in recipe order, display_order = index. */
  members: SeriesMemberInput[];
  /** Recipe samples (sample_ids, recipe order) that could not be resolved to
   *  an exposure and were therefore left off the plate. */
  unresolvedSampleIds: number[];
}

export function buildPlateFromRecipe(
  samples: SeriesSample[],
  oldMembers: SeriesMember[],
  exposureIdForSample: (sampleId: number) => number | null | undefined,
): PlateFromRecipe {
  // Consume-once pool of old members keyed by exposure_id: a carried member
  // donates its display props to exactly one new slot.
  const pool = new Map<number, SeriesMember[]>();
  for (const m of oldMembers) {
    if (m.exposure_id === null) continue;
    const bucket = pool.get(m.exposure_id);
    if (bucket) bucket.push(m);
    else pool.set(m.exposure_id, [m]);
  }

  const ordered = [...samples]
    .filter((s) => !s.excluded)
    .sort((a, b) => a.position - b.position || a.id - b.id);

  const members: SeriesMemberInput[] = [];
  const unresolvedSampleIds: number[] = [];
  for (const s of ordered) {
    const eid = exposureIdForSample(s.sample_id);
    if (eid === null || eid === undefined) {
      unresolvedSampleIds.push(s.sample_id);
      continue;
    }
    const carried = pool.get(eid)?.shift();
    members.push(
      carried
        ? { ...memberToInput(carried), display_order: members.length }
        : { exposure_id: eid, display_order: members.length },
    );
  }
  return { members, unresolvedSampleIds };
}
