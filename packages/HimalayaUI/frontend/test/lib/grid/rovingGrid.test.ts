import { describe, it, expect } from "vitest";
import {
  nextGridCoord,
  type GridCoord,
  type GridDims,
} from "../../../src/lib/grid/rovingGrid";

// Realistic Samples grid: header row 0 + 139 data rows, 6 columns.
const samples: GridDims = { rows: 140, cols: 6 };
// Small grid for edge-clamping cases.
const small: GridDims = { rows: 4, cols: 3 };

const mid: GridCoord = { row: 2, col: 1 };

describe("nextGridCoord — arrow movement (mid-grid)", () => {
  const cases: Array<[string, GridCoord]> = [
    ["ArrowUp", { row: 1, col: 1 }],
    ["ArrowDown", { row: 3, col: 1 }],
    ["ArrowLeft", { row: 2, col: 0 }],
    ["ArrowRight", { row: 2, col: 2 }],
  ];
  for (const [key, expected] of cases) {
    it(`${key} moves by 1 on the right axis`, () => {
      expect(nextGridCoord(mid, key, samples)).toEqual(expected);
    });
  }
});

describe("nextGridCoord — edge clamping (returns same coord, NOT null)", () => {
  it("ArrowUp at row 0 stays put", () => {
    const c = { row: 0, col: 2 };
    expect(nextGridCoord(c, "ArrowUp", small)).toEqual({ row: 0, col: 2 });
  });
  it("ArrowDown at last row stays put", () => {
    const c = { row: small.rows - 1, col: 2 };
    expect(nextGridCoord(c, "ArrowDown", small)).toEqual({
      row: small.rows - 1,
      col: 2,
    });
  });
  it("ArrowLeft at col 0 stays put", () => {
    const c = { row: 2, col: 0 };
    expect(nextGridCoord(c, "ArrowLeft", small)).toEqual({ row: 2, col: 0 });
  });
  it("ArrowRight at last col stays put", () => {
    const c = { row: 2, col: small.cols - 1 };
    expect(nextGridCoord(c, "ArrowRight", small)).toEqual({
      row: 2,
      col: small.cols - 1,
    });
  });
});

describe("nextGridCoord — Home / End", () => {
  it("Home → column 0 in the same row", () => {
    expect(nextGridCoord({ row: 70, col: 4 }, "Home", samples)).toEqual({
      row: 70,
      col: 0,
    });
  });
  it("End → last column in the same row", () => {
    expect(nextGridCoord({ row: 70, col: 1 }, "End", samples)).toEqual({
      row: 70,
      col: samples.cols - 1,
    });
  });
  it("Ctrl+Home → grid origin {0,0}", () => {
    expect(
      nextGridCoord({ row: 70, col: 4 }, "Home", samples, { ctrl: true }),
    ).toEqual({ row: 0, col: 0 });
  });
  it("Ctrl+End → last cell {rows-1, cols-1}", () => {
    expect(
      nextGridCoord({ row: 5, col: 1 }, "End", samples, { ctrl: true }),
    ).toEqual({ row: samples.rows - 1, col: samples.cols - 1 });
  });
});

describe("nextGridCoord — PageUp / PageDown", () => {
  it("PageDown moves +pageSize (default 10)", () => {
    expect(nextGridCoord({ row: 20, col: 3 }, "PageDown", samples)).toEqual({
      row: 30,
      col: 3,
    });
  });
  it("PageUp moves -pageSize (default 10)", () => {
    expect(nextGridCoord({ row: 20, col: 3 }, "PageUp", samples)).toEqual({
      row: 10,
      col: 3,
    });
  });
  it("PageDown clamps to last row", () => {
    expect(nextGridCoord({ row: 135, col: 0 }, "PageDown", samples)).toEqual({
      row: samples.rows - 1,
      col: 0,
    });
  });
  it("PageUp clamps to row 0", () => {
    expect(nextGridCoord({ row: 3, col: 0 }, "PageUp", samples)).toEqual({
      row: 0,
      col: 0,
    });
  });
  it("honors a custom pageSize", () => {
    expect(
      nextGridCoord({ row: 20, col: 0 }, "PageDown", samples, { pageSize: 25 }),
    ).toEqual({ row: 45, col: 0 });
  });
});

describe("nextGridCoord — non-nav keys return null", () => {
  for (const key of ["Enter", "Tab", "a", " ", "Escape"]) {
    it(`${JSON.stringify(key)} returns null`, () => {
      expect(nextGridCoord(mid, key, samples)).toBeNull();
    });
  }
});

describe("nextGridCoord — purity", () => {
  it("does not mutate the input coord", () => {
    const input = { row: 2, col: 1 };
    nextGridCoord(input, "ArrowDown", samples);
    expect(input).toEqual({ row: 2, col: 1 });
  });
  it("returns a new object (not the same reference)", () => {
    const input = { row: 0, col: 0 };
    const out = nextGridCoord(input, "ArrowUp", samples);
    expect(out).not.toBe(input);
    expect(out).toEqual({ row: 0, col: 0 });
  });
});
