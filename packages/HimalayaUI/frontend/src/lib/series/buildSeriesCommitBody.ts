/**
 * Pure `SeriesMember[] → CommitSeriesPlateBody` builder for the plate-commit
 * round-trip (I3.5b). Plate members are POSITIONAL and id-less on the wire —
 * `series_plate_committed` mints the `series_members` ids dispatcher-side
 * (master plan §5.2), so the local member ids are stripped.
 *
 * The body NO LONGER carries `expected_content_hash` (Plan 6a — the commit
 * route is last-write-wins; the backend stopped 409ing on a stale hash, so the
 * frontend stopped sending it).
 *
 * NOTE: the builder page now commits via `buildPlateFromRecipe` (BU-RECIPENOOP),
 * so this whole-members builder has no production caller. It stays, deliberately,
 * as the documented home of `memberToInput` (the ONE SeriesMember -> wire
 * projection) and the wire-contract pin its test suite exercises.
 *
 * `memberToInput` is the ONE `SeriesMember → SeriesMemberInput` wire
 * projection. `buildPlateFromRecipe` (the builder page's recipe→plate
 * resolver, BU-RECIPENOOP) reuses it for carried-over members — do not
 * duplicate the field list anywhere else.
 */
import type {
  CommitSeriesPlateBody, SeriesMember, SeriesMemberInput,
} from "../../api";

/** Project one resolved member onto the id-less wire shape. `snapshot` is
 *  sent when present, else omitted (server-filled by the commit route). */
export function memberToInput(m: SeriesMember): SeriesMemberInput {
  const out: SeriesMemberInput = {
    exposure_id: m.exposure_id,
    display_order: m.display_order,
    band_height: m.band_height,
    y_offset: m.y_offset,
    normalization: m.normalization,
    color_override: m.color_override,
    label_override: m.label_override,
    q_window_min: m.q_window_min,
    q_window_max: m.q_window_max,
    peak_display: m.peak_display,
  };
  if (m.snapshot !== null) out.snapshot = m.snapshot;
  return out;
}

export function buildSeriesCommitBody(
  members: SeriesMember[],
): CommitSeriesPlateBody {
  const plate: SeriesMemberInput[] = [...members]
    .sort((a, b) => a.display_order - b.display_order)
    .map(memberToInput);
  return { members: plate };
}
