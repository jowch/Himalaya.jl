import { describe, it, expect } from "vitest";
import {
  filterSort,
  folioControlsFromParams,
  writeFolioControlsToParams,
} from "../src/lib/series/folioFilter";
import type { SeriesSummary } from "../src/api";

function summary(over: Partial<SeriesSummary> = {}): SeriesSummary {
  return {
    id: 1, title: "Alpha", description: null, content_hash: "h",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, last_event_at: "2026-05-02 10:00:00",
    author_username: "jc", member_count: 3, member_phases: ["Pn3m"],
    member_phase_count: 1, has_stale_members: false,
    ordering_variable: null, spans_experiments: false, ...over,
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

  it("'cross' filter keeps only series that span more than one experiment", () => {
    const rows = [
      summary({ id: 1, spans_experiments: false }),
      summary({ id: 2, spans_experiments: true }),
      summary({ id: 3, spans_experiments: true }),
    ];
    const out = filterSort(rows, { search: "", filter: "cross", sort: "recent" });
    expect(out.map((r) => r.id)).toEqual([2, 3]);
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

  it("'variable' sort orders by the recipe ordering variable, nulls last", () => {
    const rows = [
      summary({ id: 1, ordering_variable: "temperature", title: "Z" }),
      summary({ id: 2, ordering_variable: null, title: "A" }),
      summary({ id: 3, ordering_variable: "dose", title: "M" }),
    ];
    const out = filterSort(rows, { search: "", filter: "all", sort: "variable" });
    // dose < temperature (by variable); the null-variable series sorts last.
    expect(out.map((r) => r.id)).toEqual([3, 1, 2]);
  });

  it("'variable' sort tiebreaks by title when the variable matches", () => {
    const rows = [
      summary({ id: 1, ordering_variable: "dose", title: "Charlie" }),
      summary({ id: 2, ordering_variable: "dose", title: "Alpha" }),
      summary({ id: 3, ordering_variable: "dose", title: "Bravo" }),
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

describe("folio URL codec (FOL P2-3)", () => {
  it("no params parse to the defaults", () => {
    expect(folioControlsFromParams(new URLSearchParams(""))).toEqual({
      search: "", filter: "all", sort: "recent",
    });
  });

  it("?q=&filter=&sort= round-trip into controls", () => {
    expect(
      folioControlsFromParams(new URLSearchParams("q=LL37&filter=transition&sort=size")),
    ).toEqual({ search: "LL37", filter: "transition", sort: "size" });
    expect(
      folioControlsFromParams(new URLSearchParams("filter=cross&sort=variable")),
    ).toEqual({ search: "", filter: "cross", sort: "variable" });
  });

  it("unknown filter/sort tokens fall back to the defaults (stale shared link still renders)", () => {
    expect(
      folioControlsFromParams(new URLSearchParams("filter=bogus&sort=nope")),
    ).toEqual({ search: "", filter: "all", sort: "recent" });
  });

  it("default controls write NO params (an empty search never leaves ?q=)", () => {
    const params = writeFolioControlsToParams(
      new URLSearchParams("q=old&filter=transition&sort=size"),
      { search: "", filter: "all", sort: "recent" },
    );
    expect(params.toString()).toBe("");
  });

  it("non-default controls write their params", () => {
    const params = writeFolioControlsToParams(new URLSearchParams(""), {
      search: "LL37", filter: "cross", sort: "variable",
    });
    expect(params.get("q")).toBe("LL37");
    expect(params.get("filter")).toBe("cross");
    expect(params.get("sort")).toBe("variable");
  });

  it("leaves foreign params untouched", () => {
    const params = writeFolioControlsToParams(new URLSearchParams("other=1&q=x"), {
      search: "", filter: "all", sort: "recent",
    });
    expect(params.get("other")).toBe("1");
    expect(params.get("q")).toBe(null);
  });

  it("controls round-trip through the params unchanged", () => {
    const controls = { search: "lipid 1-2", filter: "transition", sort: "size" } as const;
    const params = writeFolioControlsToParams(new URLSearchParams(""), controls);
    expect(folioControlsFromParams(params)).toEqual(controls);
  });
});
