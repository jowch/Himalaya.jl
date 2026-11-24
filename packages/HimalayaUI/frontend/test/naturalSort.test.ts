import { describe, expect, it } from "vitest";
import { naturalCompare, naturalSortBy } from "../src/lib/naturalSort";

describe("natural sort", () => {
  it("orders digit runs by their value, so X_02A comes before X_10A", () => {
    expect(naturalSortBy(["X_10A", "X_02A", "X_1B", "X_100"], (s) => s)).toEqual(["X_1B", "X_02A", "X_10A", "X_100"]);
  });

  it("ignores case, and breaks a tie between names that read the same so the order is stable", () => {
    expect(naturalSortBy(["b", "A", "a"], (s) => s)).toEqual(["A", "a", "b"]);
    expect(naturalCompare("X_2", "X_02")).not.toBe(0);
    expect(naturalCompare("X_2", "X_02")).toBe(-naturalCompare("X_02", "X_2"));
  });

  it("sorts by the key it is given without changing the input", () => {
    const items = [{ name: "s10" }, { name: "s9" }];
    expect(naturalSortBy(items, (s) => s.name).map((s) => s.name)).toEqual(["s9", "s10"]);
    expect(items.map((s) => s.name)).toEqual(["s10", "s9"]);
  });
});
