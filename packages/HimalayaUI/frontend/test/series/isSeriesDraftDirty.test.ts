import { describe, it, expect, beforeEach } from "vitest";
import { __resetOptimisticIdForTest } from "../../src/lib/queue/optimisticId";
import { isSeriesDraftDirty } from "../../src/lib/series/isSeriesDraftDirty";
import {
  fromSeries,
  addSampleToRecipe,
  removeRecipeRow,
  reorderRecipe,
} from "../../src/lib/series/seriesDraftFactories";
import type { Series, SeriesSample } from "../../src/api";

function sample(over: Partial<SeriesSample> = {}): SeriesSample {
  return { id: 1, series_id: 5, sample_id: 10, position: 0, pinned: false, excluded: false, ...over };
}

function series(over: Partial<Series> = {}): Series {
  return {
    id: 5, title: "Titration", description: "desc", content_hash: "sha256:base",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, forked_from_title: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    ordering_variable: "LL37 : lipid", order_rule: "ascending",
    state: "committed", members: [],
    samples: [
      sample({ id: 11, sample_id: 10, position: 0 }),
      sample({ id: 12, sample_id: 20, position: 1 }),
    ],
    ...over,
  };
}

describe("isSeriesDraftDirty", () => {
  beforeEach(() => __resetOptimisticIdForTest());

  it("a pristine fork (Edit opened, nothing changed) is NOT dirty", () => {
    const s = series();
    expect(isSeriesDraftDirty(fromSeries(s), s)).toBe(false);
  });

  it("a title edit is dirty", () => {
    const s = series();
    expect(isSeriesDraftDirty({ ...fromSeries(s), title: "Renamed" }, s)).toBe(true);
  });

  it("adding a sample is dirty", () => {
    const s = series();
    expect(isSeriesDraftDirty(addSampleToRecipe(fromSeries(s), 30), s)).toBe(true);
  });

  it("removing a sample is dirty", () => {
    const s = series();
    const d = fromSeries(s);
    expect(isSeriesDraftDirty(removeRecipeRow(d, d.recipe[0]!.id), s)).toBe(true);
  });

  it("a manual reorder is dirty", () => {
    const s = series();
    expect(isSeriesDraftDirty(reorderRecipe(fromSeries(s), 0, 1), s)).toBe(true);
  });

  it("toggling a member's excluded flag is dirty", () => {
    const s = series();
    const d = fromSeries(s);
    const flipped = { ...d, recipe: d.recipe.map((r, i) => (i === 0 ? { ...r, excluded: true } : r)) };
    expect(isSeriesDraftDirty(flipped, s)).toBe(true);
  });

  it("a description change (empty string vs null both clear) is canonicalized — '' equals a null-described series", () => {
    // buildSeriesSaveBody maps "" → null, so a blank draft over a null-description
    // series is NOT a spurious dirty.
    const s = series({ description: null });
    const d = fromSeries(s); // fromSeries seeds description "" when null
    expect(isSeriesDraftDirty(d, s)).toBe(false);
  });
});
