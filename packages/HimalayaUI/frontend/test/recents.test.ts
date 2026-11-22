import { describe, it, expect, beforeEach, vi } from "vitest";
import { byRecent, MAX_RECENTS, recentRank, recordRecent } from "../src/lib/recents";

describe("recent-use cache", () => {
  beforeEach(() => localStorage.clear());

  it("keeps keys newest first, without repeats, per kind, capped", () => {
    recordRecent("label", "buffer");
    recordRecent("label", "37C");
    recordRecent("label", "buffer");
    recordRecent("comparison:1", "12");
    expect(recentRank("label")).toEqual(["buffer", "37C"]);
    expect(recentRank("comparison:1")).toEqual(["12"]);
    expect(recentRank("comparison:2")).toEqual([]);
    for (let i = 0; i < MAX_RECENTS + 5; i++) recordRecent("label", `l${i}`);
    expect(recentRank("label")).toHaveLength(MAX_RECENTS);
    expect(recentRank("label")[0]).toBe(`l${MAX_RECENTS + 4}`);
  });

  it("reads garbage or a throwing store as empty and never throws", () => {
    localStorage.setItem("himalaya-ui:recent:label", "{nope");
    expect(recentRank("label")).toEqual([]);
    localStorage.setItem("himalaya-ui:recent:label", JSON.stringify(["a", 3, null]));
    expect(recentRank("label")).toEqual(["a"]);
    const get = vi.spyOn(Storage.prototype, "getItem").mockImplementation(() => { throw new Error("blocked"); });
    const set = vi.spyOn(Storage.prototype, "setItem").mockImplementation(() => { throw new Error("blocked"); });
    expect(recentRank("label")).toEqual([]);
    expect(() => recordRecent("label", "x")).not.toThrow();
    get.mockRestore();
    set.mockRestore();
  });

  it("byRecent puts recent items first in recency order and keeps the rest in their order", () => {
    const items = ["a", "b", "c", "d"];
    expect(byRecent(items, (x) => x, ["c", "a"])).toEqual(["c", "a", "b", "d"]);
    expect(byRecent(items, (x) => x, [])).toEqual(items);
  });
});
