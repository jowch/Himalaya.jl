import { describe, it, expect } from "vitest";
import { tileColumns } from "../../src/ui/tileColumns";

describe("tileColumns", () => {
  it("keeps six across when the plot leaves little room below it", () => {
    expect(tileColumns(760, 120)).toBe(6);
    expect(tileColumns(426, 200)).toBe(6);
  });

  it("goes two across when two rows of half-width tiles fit below the plot", () => {
    expect(tileColumns(760, 900)).toBe(2);
    expect(tileColumns(426, 520)).toBe(2);
  });

  it("steps through three and four across as the free height shrinks", () => {
    expect(tileColumns(760, 600)).toBe(3);
    expect(tileColumns(760, 460)).toBe(4);
  });

  it("only uses fewer columns as the column narrows, never more", () => {
    const free = 600;
    const widths = [1000, 900, 800, 700, 600, 500, 400, 300];
    const cols = widths.map((w) => tileColumns(w, free));
    cols.slice(1).forEach((c, i) => expect(c).toBeLessThanOrEqual(cols[i]!));
  });

  it("falls back to six across before the column has a width", () => {
    expect(tileColumns(0, 2000)).toBe(6);
  });
});
