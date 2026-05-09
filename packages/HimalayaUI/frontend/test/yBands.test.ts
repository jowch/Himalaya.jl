/**
 * Tests for the shared y-band layout math (extracted from MultiTracePlot.tsx
 * as part of #93 / #99 follow-up — see src/lib/comparison/yBands.ts).
 */
import { describe, it, expect } from "vitest";
import { computeYBands } from "../src/lib/comparison/yBands";

describe("computeYBands", () => {
  it("returns one band per member, summing to panelHeight", () => {
    const bands = computeYBands([1, 1, 1], 300);
    expect(bands).toHaveLength(3);
    expect(bands[0]).toEqual([0, 100]);
    expect(bands[1]).toEqual([100, 200]);
    expect(bands[2]).toEqual([200, 300]);
  });

  it("respects unequal band_height ratios", () => {
    const bands = computeYBands([1, 2, 1], 400);
    // ratio sum = 4; each ratio*100 = pixel slice.
    expect(bands[0]).toEqual([0, 100]);
    expect(bands[1]).toEqual([100, 300]);
    expect(bands[2]).toEqual([300, 400]);
  });

  it("preserves caller-supplied order (display_order)", () => {
    // Same total members but different input order yields different bands.
    const a = computeYBands([1, 2], 300);
    const b = computeYBands([2, 1], 300);
    expect(a[0]![1]).toBe(100);  // first band 100 px tall when ratio = 1
    expect(b[0]![1]).toBe(200);  // first band 200 px tall when ratio = 2
  });

  it("returns an empty array for no members", () => {
    expect(computeYBands([], 300)).toEqual([]);
  });

  it("handles a zero panel height gracefully (degenerate render frame)", () => {
    const bands = computeYBands([1, 1], 0);
    expect(bands[0]).toEqual([0, 0]);
    expect(bands[1]).toEqual([0, 0]);
  });

  it("falls back to equal slicing when every ratio is zero (total ≤ 0)", () => {
    // Without this guard, total=0 produces NaN bands. Equal slicing keeps
    // the plot well-defined even if a comparison has all-zero band_heights.
    const bands = computeYBands([0, 0, 0], 300);
    expect(bands).toEqual([[0, 100], [100, 200], [200, 300]]);
  });
});
