import { describe, it, expect, beforeEach } from "vitest";
import {
  emptySeriesDraft,
  persistSeriesDraftToSession,
  loadSeriesDraftFromSession,
  SERIES_DRAFT_KEY,
  type SeriesDraft,
} from "../../src/lib/series/seriesDraft";

function draft(over: Partial<SeriesDraft> = {}): SeriesDraft {
  return {
    ...emptySeriesDraft(5),
    ...over,
  };
}

describe("seriesDraft", () => {
  beforeEach(() => sessionStorage.clear());

  it("emptySeriesDraft seeds id, empty recipe, default order_rule", () => {
    const d = emptySeriesDraft(7);
    expect(d.id).toBe(7);
    expect(d.recipe).toEqual([]);
    expect(d.orderRule).toBe("manual");
    expect(d.orderingVariable).toBeNull();
  });

  it("round-trips through sessionStorage", () => {
    const d = draft({
      title: "Titration",
      recipe: [{ id: 1, sample_id: 10, position: 0, pinned: false, excluded: false }],
    });
    persistSeriesDraftToSession(d);
    expect(loadSeriesDraftFromSession()).toEqual(d);
  });

  it("persisting null clears the key", () => {
    persistSeriesDraftToSession(draft());
    persistSeriesDraftToSession(null);
    expect(sessionStorage.getItem(SERIES_DRAFT_KEY)).toBeNull();
    expect(loadSeriesDraftFromSession()).toBeNull();
  });

  it("loadSeriesDraftFromSession returns null on a version mismatch", () => {
    sessionStorage.setItem(
      SERIES_DRAFT_KEY,
      JSON.stringify({ version: 999, draft: draft() }),
    );
    expect(loadSeriesDraftFromSession()).toBeNull();
  });

  it("loadSeriesDraftFromSession returns null on garbage", () => {
    sessionStorage.setItem(SERIES_DRAFT_KEY, "{not json");
    expect(loadSeriesDraftFromSession()).toBeNull();
  });
});
