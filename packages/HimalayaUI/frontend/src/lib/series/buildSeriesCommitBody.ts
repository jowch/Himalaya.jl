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
 * The plate source is the loaded series' current `members` (the recipe-resolved
 * exposures). `snapshot` is sent when present, else omitted (server-filled).
 */
import type {
  CommitSeriesPlateBody, SeriesMember, SeriesMemberInput,
} from "../../api";

export function buildSeriesCommitBody(
  members: SeriesMember[],
): CommitSeriesPlateBody {
  const plate: SeriesMemberInput[] = [...members]
    .sort((a, b) => a.display_order - b.display_order)
    .map((m) => {
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
    });
  return { members: plate };
}
