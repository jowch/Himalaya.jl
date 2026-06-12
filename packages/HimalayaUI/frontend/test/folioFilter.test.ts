import { describe, it, expect } from "vitest";
import {
  filterSort,
  folioControlsFromParams,
  writeFolioControlsToParams,
  defaultDirForSort,
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
    ordering_variable: null, spans_experiments: false, experiment_name: null, ...over,
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

  describe("broadened search haystack (F1)", () => {
    it("matches a member-phase token (e.g. 'Im3m')", () => {
      const rows = [
        summary({ id: 1, title: "Cubic run", member_phases: ["Pn3m", "Im3m"] }),
        summary({ id: 2, title: "Lamellar run", member_phases: ["Lamellar"] }),
      ];
      const out = filterSort(rows, { search: "im3m", filter: "all", sort: "recent" });
      expect(out.map((r) => r.id)).toEqual([1]);
    });

    it("matches an ordering-variable token (e.g. part of 'lipid ratio')", () => {
      const rows = [
        summary({ id: 1, title: "Run A", ordering_variable: "LL37 : lipid ratio" }),
        summary({ id: 2, title: "Run B", ordering_variable: "temperature" }),
      ];
      const out = filterSort(rows, { search: "lipid", filter: "all", sort: "recent" });
      expect(out.map((r) => r.id)).toEqual([1]);
    });

    it("matches an author-username token", () => {
      const rows = [
        summary({ id: 1, title: "Run A", author_username: "marie" }),
        summary({ id: 2, title: "Run B", author_username: "jc" }),
      ];
      const out = filterSort(rows, { search: "marie", filter: "all", sort: "recent" });
      expect(out.map((r) => r.id)).toEqual([1]);
    });

    it("matches a description word", () => {
      const rows = [
        summary({ id: 1, title: "Run A", description: "the swelling titration baseline" }),
        summary({ id: 2, title: "Run B", description: null }),
      ];
      const out = filterSort(rows, { search: "swelling", filter: "all", sort: "recent" });
      expect(out.map((r) => r.id)).toEqual([1]);
    });

    it("a token matching NONE of the searched fields returns empty", () => {
      const rows = [
        summary({ id: 1, title: "Run A", description: "baseline", author_username: "jc",
          member_phases: ["Pn3m"], ordering_variable: "dose" }),
      ];
      const out = filterSort(rows, { search: "zzznomatch", filter: "all", sort: "recent" });
      expect(out).toHaveLength(0);
    });
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

describe("filterSort sort direction (F2)", () => {
  it("defaultDirForSort: recent/size default desc, variable defaults asc", () => {
    expect(defaultDirForSort("recent")).toBe("desc");
    expect(defaultDirForSort("size")).toBe("desc");
    expect(defaultDirForSort("variable")).toBe("asc");
  });

  it("'size' asc reverses the default (smallest first)", () => {
    const rows = [
      summary({ id: 1, member_count: 2 }),
      summary({ id: 2, member_count: 9 }),
      summary({ id: 3, member_count: 5 }),
    ];
    const out = filterSort(rows, { search: "", filter: "all", sort: "size", dir: "asc" });
    expect(out.map((r) => r.id)).toEqual([1, 3, 2]);
  });

  it("'size' desc (the default) is unchanged when dir is given explicitly", () => {
    const rows = [
      summary({ id: 1, member_count: 2 }),
      summary({ id: 2, member_count: 9 }),
      summary({ id: 3, member_count: 5 }),
    ];
    const out = filterSort(rows, { search: "", filter: "all", sort: "size", dir: "desc" });
    expect(out.map((r) => r.id)).toEqual([2, 3, 1]);
  });

  it("'variable' desc reverses the A→Z order (Z→A, nulls FIRST)", () => {
    const rows = [
      summary({ id: 1, ordering_variable: "temperature", title: "Z" }),
      summary({ id: 2, ordering_variable: null, title: "A" }),
      summary({ id: 3, ordering_variable: "dose", title: "M" }),
    ];
    // asc default: [dose(3), temperature(1), null(2)] → desc reverses to [2, 1, 3]
    const out = filterSort(rows, { search: "", filter: "all", sort: "variable", dir: "desc" });
    expect(out.map((r) => r.id)).toEqual([2, 1, 3]);
  });

  it("'recent' desc (default) preserves backend order; asc reverses it", () => {
    const rows = [summary({ id: 3 }), summary({ id: 1 }), summary({ id: 2 })];
    const desc = filterSort(rows, { search: "", filter: "all", sort: "recent", dir: "desc" });
    expect(desc.map((r) => r.id)).toEqual([3, 1, 2]);
    const asc = filterSort(rows, { search: "", filter: "all", sort: "recent", dir: "asc" });
    expect(asc.map((r) => r.id)).toEqual([2, 1, 3]);
  });

  it("omitting dir uses the sort's default direction (back-compat)", () => {
    const rows = [
      summary({ id: 1, member_count: 2 }),
      summary({ id: 2, member_count: 9 }),
    ];
    const out = filterSort(rows, { search: "", filter: "all", sort: "size" });
    expect(out.map((r) => r.id)).toEqual([2, 1]);
  });
});

describe("folio URL codec (FOL P2-3)", () => {
  it("no params parse to the defaults", () => {
    expect(folioControlsFromParams(new URLSearchParams(""))).toEqual({
      search: "", filter: "all", sort: "recent", dir: "desc",
    });
  });

  it("?q=&filter=&sort= round-trip into controls", () => {
    expect(
      folioControlsFromParams(new URLSearchParams("q=LL37&filter=transition&sort=size")),
    ).toEqual({ search: "LL37", filter: "transition", sort: "size", dir: "desc" });
    expect(
      folioControlsFromParams(new URLSearchParams("filter=cross&sort=variable")),
    ).toEqual({ search: "", filter: "cross", sort: "variable", dir: "asc" });
  });

  it("unknown filter/sort tokens fall back to the defaults (stale shared link still renders)", () => {
    expect(
      folioControlsFromParams(new URLSearchParams("filter=bogus&sort=nope")),
    ).toEqual({ search: "", filter: "all", sort: "recent", dir: "desc" });
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
    const controls = { search: "lipid 1-2", filter: "transition", sort: "size", dir: "desc" } as const;
    const params = writeFolioControlsToParams(new URLSearchParams(""), controls);
    expect(folioControlsFromParams(params)).toEqual(controls);
  });

  describe("dir param (F2)", () => {
    it("absent dir parses to the sort's default direction", () => {
      expect(folioControlsFromParams(new URLSearchParams("")).dir).toBe("desc"); // recent → desc
      expect(folioControlsFromParams(new URLSearchParams("sort=variable")).dir).toBe("asc");
      expect(folioControlsFromParams(new URLSearchParams("sort=size")).dir).toBe("desc");
    });

    it("an explicit non-default dir parses through", () => {
      expect(
        folioControlsFromParams(new URLSearchParams("sort=size&dir=asc")).dir,
      ).toBe("asc");
      expect(
        folioControlsFromParams(new URLSearchParams("sort=variable&dir=desc")).dir,
      ).toBe("desc");
    });

    it("an unknown dir token falls back to the sort's default", () => {
      expect(
        folioControlsFromParams(new URLSearchParams("sort=size&dir=sideways")).dir,
      ).toBe("desc");
    });

    it("dir equal to the sort default writes NO dir param", () => {
      const params = writeFolioControlsToParams(new URLSearchParams(""), {
        search: "", filter: "all", sort: "size", dir: "desc",
      });
      expect(params.get("dir")).toBe(null);
    });

    it("dir differing from the sort default writes the dir param", () => {
      const params = writeFolioControlsToParams(new URLSearchParams(""), {
        search: "", filter: "all", sort: "size", dir: "asc",
      });
      expect(params.get("dir")).toBe("asc");
    });

    it("a non-default dir round-trips", () => {
      const controls = { search: "", filter: "all", sort: "variable", dir: "desc" } as const;
      const params = writeFolioControlsToParams(new URLSearchParams(""), controls);
      expect(folioControlsFromParams(params)).toEqual(controls);
    });
  });
});
