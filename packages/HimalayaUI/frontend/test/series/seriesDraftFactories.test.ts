import { describe, it, expect, beforeEach } from "vitest";
import { __resetOptimisticIdForTest } from "../../src/lib/queue/optimisticId";
import {
  fromSeries,
  addSampleToRecipe,
  reorderRecipe,
} from "../../src/lib/series/seriesDraftFactories";
import { buildSeriesSaveBody } from "../../src/lib/series/buildSeriesSaveBody";
import { buildSeriesCommitBody } from "../../src/lib/series/buildSeriesCommitBody";
import type { Series, SeriesSample, SeriesMember } from "../../src/api";

function sample(over: Partial<SeriesSample> = {}): SeriesSample {
  return { id: 1, series_id: 5, sample_id: 10, position: 0, pinned: false, excluded: false, ...over };
}

function member(over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id: 1, series_id: 5, exposure_id: 101, display_order: 0,
    band_height: 1, y_offset: 0, normalization: "max",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: null, is_stale: false, created_by: 1, created_at: null, ...over,
  };
}

function series(over: Partial<Series> = {}): Series {
  return {
    id: 5, title: "Titration", description: "desc", content_hash: "sha256:base",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, forked_from_title: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    ordering_variable: "LL37 : lipid", order_rule: "ascending",
    state: "committed", members: [], samples: [], ...over,
  };
}

describe("fromSeries", () => {
  it("seeds the draft from a Series with real positive recipe ids", () => {
    const s = series({
      samples: [sample({ id: 11, sample_id: 10, position: 0 }), sample({ id: 12, sample_id: 20, position: 1 })],
    });
    const d = fromSeries(s);
    expect(d.id).toBe(5);
    expect(d.title).toBe("Titration");
    expect(d.description).toBe("desc");
    expect(d.orderingVariable).toBe("LL37 : lipid");
    expect(d.orderRule).toBe("ascending");
    expect(d.recipe.map((r) => [r.id, r.sample_id, r.position])).toEqual([
      [11, 10, 0], [12, 20, 1],
    ]);
  });

  it("coerces an out-of-union order_rule to manual", () => {
    const d = fromSeries(series({ order_rule: "weird" }));
    expect(d.orderRule).toBe("manual");
  });

  it("null description maps to an empty string", () => {
    const d = fromSeries(series({ description: null }));
    expect(d.description).toBe("");
  });
});

describe("addSampleToRecipe", () => {
  beforeEach(() => __resetOptimisticIdForTest());

  it("appends a row with a NEGATIVE placeholder id from nextOptimisticId()", () => {
    const d0 = fromSeries(series({ samples: [sample({ id: 11, sample_id: 10, position: 0 })] }));
    const d1 = addSampleToRecipe(d0, 20);
    expect(d1.recipe).toHaveLength(2);
    const added = d1.recipe[1];
    expect(added.sample_id).toBe(20);
    expect(added.position).toBe(1); // appended at the tail
    expect(added.id).toBeLessThan(0); // placeholder
    expect(added.pinned).toBe(false);
    expect(added.excluded).toBe(false);
  });

  it("does not mutate the input draft", () => {
    const d0 = fromSeries(series({ samples: [] }));
    addSampleToRecipe(d0, 20);
    expect(d0.recipe).toHaveLength(0);
  });
});

describe("reorderRecipe", () => {
  it("moves a row and renumbers positions densely", () => {
    let d = fromSeries(series({
      samples: [
        sample({ id: 11, sample_id: 10, position: 0 }),
        sample({ id: 12, sample_id: 20, position: 1 }),
        sample({ id: 13, sample_id: 30, position: 2 }),
      ],
    }));
    d = reorderRecipe(d, 2, 0); // move 3rd → front
    expect(d.recipe.map((r) => r.sample_id)).toEqual([30, 10, 20]);
    expect(d.recipe.map((r) => r.position)).toEqual([0, 1, 2]);
  });
});

describe("buildSeriesSaveBody", () => {
  it("emits positional id-less SeriesSampleInput[] and carries ordering meta", () => {
    let d = fromSeries(series({ samples: [sample({ id: 11, sample_id: 10, position: 0 })] }));
    d = addSampleToRecipe(d, 20); // placeholder negative id
    const body = buildSeriesSaveBody(d);
    expect(body.title).toBe("Titration");
    expect(body.samples).toEqual([
      { sample_id: 10, position: 0, pinned: false, excluded: false },
      { sample_id: 20, position: 1, pinned: false, excluded: false },
    ]);
    // No id leaks onto the wire (membership is positional).
    for (const s of body.samples) expect("id" in s).toBe(false);
    expect(body.ordering_variable).toBe("LL37 : lipid");
    expect(body.order_rule).toBe("ascending");
  });

  it("NEVER carries expected_content_hash (recipe-save cannot conflict)", () => {
    const d = fromSeries(series());
    const body = buildSeriesSaveBody(d);
    expect("expected_content_hash" in body).toBe(false);
  });
});

describe("buildSeriesCommitBody", () => {
  it("emits id-less positional members and NEVER carries expected_content_hash (LWW relax, Plan 6a)", () => {
    const s = series({
      members: [member({ id: 1, exposure_id: 101, display_order: 0 }),
                member({ id: 2, exposure_id: 102, display_order: 1 })],
    });
    const body = buildSeriesCommitBody(s.members);
    expect(body).not.toHaveProperty("expected_content_hash");
    expect(body.members.map((m) => [m.exposure_id, m.display_order])).toEqual([
      [101, 0], [102, 1],
    ]);
    for (const m of body.members) expect("id" in m).toBe(false);
  });
});
