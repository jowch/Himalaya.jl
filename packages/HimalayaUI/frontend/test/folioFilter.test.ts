import { describe, it, expect } from "vitest";
import { filterSort } from "../src/lib/series/folioFilter";
import type { SeriesSummary } from "../src/api";

function summary(over: Partial<SeriesSummary> = {}): SeriesSummary {
  return {
    id: 1, title: "Alpha", description: null, content_hash: "h",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, last_event_at: "2026-05-02 10:00:00",
    author_username: "jc", member_count: 3, member_phases: ["Pn3m"],
    member_phase_count: 1, has_stale_members: false, ...over,
  };
}

describe("filterSort", () => {
  it("returns everything for the default controls", () => {
    const rows = [summary({ id: 1 }), summary({ id: 2 })];
    expect(filterSort(rows, { search: "", filter: "all", sort: "recent" })).toHaveLength(2);
  });

  it("filters by case-insensitive title search", () => {
    const rows = [summary({ id: 1, title: "Alpha series" }), summary({ id: 2, title: "Beta series" })];
    const out = filterSort(rows, { search: "BETA", filter: "all", sort: "recent" });
    expect(out.map((r) => r.id)).toEqual([2]);
  });

  it("'transition' filter keeps only multi-phase series", () => {
    const rows = [
      summary({ id: 1, member_phase_count: 1 }),
      summary({ id: 2, member_phase_count: 3 }),
    ];
    const out = filterSort(rows, { search: "", filter: "transition", sort: "recent" });
    expect(out.map((r) => r.id)).toEqual([2]);
  });

  it("'cross' filter yields nothing (no cross-experiment data on the listing)", () => {
    const rows = [summary({ id: 1 }), summary({ id: 2 })];
    expect(filterSort(rows, { search: "", filter: "cross", sort: "recent" })).toHaveLength(0);
  });

  it("'size' sort orders by member_count descending", () => {
    const rows = [
      summary({ id: 1, member_count: 2 }),
      summary({ id: 2, member_count: 9 }),
      summary({ id: 3, member_count: 5 }),
    ];
    const out = filterSort(rows, { search: "", filter: "all", sort: "size" });
    expect(out.map((r) => r.id)).toEqual([2, 3, 1]);
  });

  it("'recent' sort preserves the backend input order", () => {
    const rows = [summary({ id: 3 }), summary({ id: 1 }), summary({ id: 2 })];
    const out = filterSort(rows, { search: "", filter: "all", sort: "recent" });
    expect(out.map((r) => r.id)).toEqual([3, 1, 2]);
  });

  it("'variable' sort orders by title (the only stable listing text key)", () => {
    const rows = [
      summary({ id: 1, title: "Charlie" }),
      summary({ id: 2, title: "Alpha" }),
      summary({ id: 3, title: "Bravo" }),
    ];
    const out = filterSort(rows, { search: "", filter: "all", sort: "variable" });
    expect(out.map((r) => r.id)).toEqual([2, 3, 1]);
  });

  it("does not mutate the input array", () => {
    const rows = [summary({ id: 1, member_count: 2 }), summary({ id: 2, member_count: 9 })];
    const snapshot = rows.map((r) => r.id);
    filterSort(rows, { search: "", filter: "all", sort: "size" });
    expect(rows.map((r) => r.id)).toEqual(snapshot);
  });
});
