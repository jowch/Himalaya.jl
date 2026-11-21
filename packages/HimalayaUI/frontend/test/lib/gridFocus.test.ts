import { describe, it, expect } from "vitest";
import { gridStep } from "../../src/lib/gridFocus";

type Box = { left: number; top: number; width: number; height: number };
const at = (b: Box): Box => b;
const grid: Box[] = [
  at({ left: 0, top: 0, width: 100, height: 50 }),
  at({ left: 120, top: 0, width: 100, height: 50 }),
  at({ left: 0, top: 70, width: 100, height: 50 }),
  at({ left: 120, top: 70, width: 100, height: 50 }),
  at({ left: 0, top: 200, width: 100, height: 50 }),
];

describe("gridStep", () => {
  it("steps left and right in reading order, clamped at the ends", () => {
    expect(gridStep(grid, 0, "ArrowRight")).toBe(1);
    expect(gridStep(grid, 1, "ArrowRight")).toBe(2);
    expect(gridStep(grid, 4, "ArrowRight")).toBe(4);
    expect(gridStep(grid, 0, "ArrowLeft")).toBe(0);
  });

  it("steps up and down to the nearest card in the next row", () => {
    expect(gridStep(grid, 1, "ArrowDown")).toBe(3);
    expect(gridStep(grid, 3, "ArrowDown")).toBe(4);
    expect(gridStep(grid, 3, "ArrowUp")).toBe(1);
    expect(gridStep(grid, 0, "ArrowUp")).toBe(0);
    expect(gridStep(grid, 4, "ArrowDown")).toBe(4);
  });

  it("falls back to reading order when the cards have no layout", () => {
    const flat = [0, 1, 2].map(() => at({ left: 0, top: 0, width: 0, height: 0 }));
    expect(gridStep(flat, 0, "ArrowDown")).toBe(1);
    expect(gridStep(flat, 1, "ArrowUp")).toBe(0);
  });

  it("starts at the first card when nothing is focused", () => {
    expect(gridStep(grid, -1, "ArrowDown")).toBe(0);
    expect(gridStep(grid, -1, "ArrowLeft")).toBe(0);
  });
});
