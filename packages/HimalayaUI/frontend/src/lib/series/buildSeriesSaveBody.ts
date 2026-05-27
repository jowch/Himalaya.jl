/**
 * Pure `SeriesDraft → SaveSeriesBody` builder for the recipe-save (PATCH)
 * round-trip (I3.5b). Membership is POSITIONAL and id-less on the wire: the
 * recipe rows' local ids (positive or negative placeholders) are stripped.
 *
 * Carries NO `expected_content_hash` — `SaveSeriesBody` has no such field and
 * `PATCH /api/series/:id` never reads it (routes_series.jl:170-214). Recipe
 * save therefore cannot 409 / conflict. (Commit is the only carrier — see
 * `buildSeriesCommitBody`.)
 */
import type { SaveSeriesBody, SeriesSampleInput } from "../../api";
import type { SeriesDraft } from "./seriesDraft";

export function buildSeriesSaveBody(draft: SeriesDraft): SaveSeriesBody {
  const samples: SeriesSampleInput[] = [...draft.recipe]
    .sort((a, b) => a.position - b.position)
    .map((r) => ({
      sample_id: r.sample_id,
      position: r.position,
      pinned: r.pinned,
      excluded: r.excluded,
    }));
  const body: SaveSeriesBody = {
    title: draft.title,
    samples,
    order_rule: draft.orderRule,
    // `ordering_variable` is `string | null` on the wire; null clears it.
    ordering_variable: draft.orderingVariable,
  };
  // Description: empty string means "no description" — send null to clear.
  body.description = draft.description === "" ? null : draft.description;
  return body;
}
