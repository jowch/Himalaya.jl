/**
 * Series draft factories + pure recipe transforms (I3.5b).
 *
 * `fromSeries` seeds a `SeriesDraft` from a loaded `Series` with REAL positive
 * recipe ids. `addSampleToRecipe` appends a row with a NEGATIVE placeholder id
 * minted by `nextOptimisticId()` (master plan §11) — a local render-key only,
 * stripped before the wire. `reorderRecipe` moves a row and renumbers positions
 * densely. All transforms are pure (return a fresh draft).
 *
 * Imports `nextOptimisticId` from the queue's shared source so a placeholder id
 * can never collide with a concurrent peak/comparison placeholder in the same
 * tab (or across tabs).
 */
import type { OrderRule, Series } from "../../api";
import { nextOptimisticId } from "../queue/optimisticId";
import type { SeriesDraft, SeriesRecipeRow } from "./seriesDraft";

const ORDER_RULES: readonly OrderRule[] = ["ascending", "descending", "manual"];

function coerceOrderRule(value: string): OrderRule {
  return ORDER_RULES.includes(value as OrderRule) ? (value as OrderRule) : "manual";
}

/** Seed an editable draft from a loaded series. Recipe rows carry the real
 *  `series_samples.id` (positive). */
export function fromSeries(s: Series): SeriesDraft {
  const recipe: SeriesRecipeRow[] = [...s.samples]
    .sort((a, b) => a.position - b.position)
    .map((ss) => ({
      id: ss.id,
      sample_id: ss.sample_id,
      position: ss.position,
      pinned: ss.pinned,
      excluded: ss.excluded,
    }));
  return {
    id: s.id,
    baseHash: s.content_hash === "" ? undefined : s.content_hash,
    title: s.title,
    description: s.description ?? "",
    orderingVariable: s.ordering_variable,
    orderRule: coerceOrderRule(s.order_rule),
    recipe,
  };
}

/** Append a sample to the recipe with a negative placeholder id. Pure. */
export function addSampleToRecipe(draft: SeriesDraft, sampleId: number): SeriesDraft {
  const row: SeriesRecipeRow = {
    id: nextOptimisticId(),
    sample_id: sampleId,
    position: draft.recipe.length,
    pinned: false,
    excluded: false,
  };
  return { ...draft, recipe: [...draft.recipe, row] };
}

/** Remove a recipe row by its (possibly-negative) local id; renumber. Pure. */
export function removeRecipeRow(draft: SeriesDraft, rowId: number): SeriesDraft {
  const recipe = draft.recipe
    .filter((r) => r.id !== rowId)
    .map((r, i) => ({ ...r, position: i }));
  return { ...draft, recipe };
}

/** Move a recipe row from `from` to `to` (array indices); renumber. Pure. */
export function reorderRecipe(draft: SeriesDraft, from: number, to: number): SeriesDraft {
  if (from === to || from < 0 || to < 0 || from >= draft.recipe.length || to >= draft.recipe.length) {
    return draft;
  }
  const next = [...draft.recipe];
  const [moved] = next.splice(from, 1);
  next.splice(to, 0, moved);
  return { ...draft, recipe: next.map((r, i) => ({ ...r, position: i })) };
}
