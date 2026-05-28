import { describe, it, expect } from "vitest";
import { seriesRatio } from "../src/lib/seriesRatio";

describe("seriesRatio", () => {
  it("formats cubic Pn3m radicands as √N terms", () => {
    expect(seriesRatio("Pn3m", [1, 2, 3])).toBe("√2 : √3 : 2");
  });
  it("formats lamellar as bare integers", () => {
    expect(seriesRatio("Lamellar", [1, 2, 3])).toBe("1 : 2 : 3");
  });
  it("formats hexagonal mixed integer/radical", () => {
    expect(seriesRatio("Hexagonal", [1, 2, 3])).toBe("1 : √3 : 2");
  });
  it("collapses perfect squares to integers (√4 → 2)", () => {
    expect(seriesRatio("Im3m", [1, 2])).toBe("√2 : 2");
  });
  it("dedupes and sorts positions", () => {
    expect(seriesRatio("Lamellar", [3, 1, 1, 2])).toBe("1 : 2 : 3");
  });
  it("caps at four terms with an ellipsis", () => {
    expect(seriesRatio("Lamellar", [1, 2, 3, 4, 5])).toBe("1 : 2 : 3 : 4 …");
  });
  it("returns empty string for no positions or unknown phase", () => {
    expect(seriesRatio("Lamellar", [])).toBe("");
    expect(seriesRatio("Nonsense", [1])).toBe("");
  });
  it("ignores positions beyond the known series length", () => {
    expect(seriesRatio("Pn3m", [1, 999])).toBe("√2");
  });
});
