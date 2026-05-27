/**
 * Pure `SeriesDraft + SeriesMember[] → CommitSeriesPlateBody` builder for the
 * plate-commit round-trip (I3.5b). Plate members are POSITIONAL and id-less on
 * the wire — `series_plate_committed` mints the `series_members` ids
 * dispatcher-side (master plan §5.2), so the local member ids are stripped.
 *
 * Carries `expected_content_hash = draft.baseHash` (the only series body that
 * carries the hash; commit is the only conflict-bearing route). Omits it when
 * `baseHash` is undefined (a never-committed draft).
 *
 * The plate source is the loaded series' current `members` (the recipe-resolved
 * exposures). `snapshot` is sent when present, else omitted (server-filled).
 */
import type {
  CommitSeriesPlateBody, SeriesMember, SeriesMemberInput,
} from "../../api";
import type { SeriesDraft } from "./seriesDraft";

export function buildSeriesCommitBody(
  draft: SeriesDraft,
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
  const body: CommitSeriesPlateBody = { members: plate };
  if (draft.baseHash !== undefined) body.expected_content_hash = draft.baseHash;
  return body;
}
