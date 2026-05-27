import { describe, it, expect, beforeEach } from "vitest";
import { useAppState } from "../src/state";
import { __resetOptimisticIdForTest } from "../src/lib/queue/optimisticId";
import type { Series, SeriesSample } from "../src/api";

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
    samples: [sample({ id: 11, sample_id: 10, position: 0 })], ...over,
  };
}

describe("series-draft Zustand actions", () => {
  beforeEach(() => {
    __resetOptimisticIdForTest();
    sessionStorage.clear();
    useAppState.getState().discardSeriesDraft();
  });

  it("startSeriesDraftFromSeries seeds the slot from the series", () => {
    useAppState.getState().startSeriesDraftFromSeries(series());
    const d = useAppState.getState().seriesDraft;
    expect(d?.id).toBe(5);
    expect(d?.title).toBe("Titration");
    expect(d?.recipe).toHaveLength(1);
    expect(d?.baseHash).toBe("sha256:base");
  });

  it("is idempotent for the same series id (preserves in-progress edits)", () => {
    const st = useAppState.getState();
    st.startSeriesDraftFromSeries(series());
    st.setSeriesDraftTitle("Edited");
    // Re-seeding the same id must NOT clobber the edit.
    st.startSeriesDraftFromSeries(series());
    expect(useAppState.getState().seriesDraft?.title).toBe("Edited");
  });

  it("addSeriesSample appends a negative-id placeholder row", () => {
    const st = useAppState.getState();
    st.startSeriesDraftFromSeries(series());
    st.addSeriesSample(20);
    const recipe = useAppState.getState().seriesDraft!.recipe;
    expect(recipe).toHaveLength(2);
    expect(recipe[1].sample_id).toBe(20);
    expect(recipe[1].id).toBeLessThan(0);
  });

  it("removeSeriesSample drops a row by local id and renumbers", () => {
    const st = useAppState.getState();
    st.startSeriesDraftFromSeries(series());
    st.addSeriesSample(20);
    const placeholderId = useAppState.getState().seriesDraft!.recipe[1].id;
    st.removeSeriesSample(placeholderId);
    const recipe = useAppState.getState().seriesDraft!.recipe;
    expect(recipe.map((r) => r.sample_id)).toEqual([10]);
    expect(recipe.map((r) => r.position)).toEqual([0]);
  });

  it("reorderSeriesSample moves a row and renumbers densely", () => {
    const st = useAppState.getState();
    st.startSeriesDraftFromSeries(series({
      samples: [
        sample({ id: 11, sample_id: 10, position: 0 }),
        sample({ id: 12, sample_id: 20, position: 1 }),
      ],
    }));
    st.reorderSeriesSample(1, 0);
    const recipe = useAppState.getState().seriesDraft!.recipe;
    expect(recipe.map((r) => r.sample_id)).toEqual([20, 10]);
    expect(recipe.map((r) => r.position)).toEqual([0, 1]);
  });

  it("ordering-variable + order-rule setters update the draft", () => {
    const st = useAppState.getState();
    st.startSeriesDraftFromSeries(series());
    st.setSeriesOrderingVariable("new var");
    st.setSeriesOrderRule("manual");
    const d = useAppState.getState().seriesDraft!;
    expect(d.orderingVariable).toBe("new var");
    expect(d.orderRule).toBe("manual");
  });

  it("mirrors to sessionStorage and clears it on discard", () => {
    const st = useAppState.getState();
    st.startSeriesDraftFromSeries(series());
    expect(sessionStorage.getItem("himalaya-ui:series-draft")).not.toBeNull();
    st.discardSeriesDraft();
    expect(sessionStorage.getItem("himalaya-ui:series-draft")).toBeNull();
    expect(useAppState.getState().seriesDraft).toBeNull();
  });

  it("recipe-mutating actions are no-ops with no active draft", () => {
    const st = useAppState.getState();
    st.discardSeriesDraft();
    st.addSeriesSample(99);
    st.setSeriesDraftTitle("nope");
    expect(useAppState.getState().seriesDraft).toBeNull();
  });

  it("seriesDraft is NOT in the persisted (localStorage) partition", () => {
    const st = useAppState.getState();
    st.startSeriesDraftFromSeries(series());
    const ls = localStorage.getItem("himalaya-ui:state");
    if (ls) expect(ls).not.toContain("seriesDraft");
  });
});
