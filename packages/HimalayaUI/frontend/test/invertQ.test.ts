/**
 * Tests for the shared `invertQ` Observable Plot helper (Plan §Phase 6,
 * Task 6.2). This is the runtime cast pattern documented in CLAUDE.md
 * frontend gotchas — the plot element exposes a `.scale("x").invert(px)`
 * runtime method that the DOM types don't know about.
 */
import { describe, it, expect } from "vitest";
import { invertQ, applyQ } from "../src/lib/plot/invertQ";

function fakePlot(): unknown {
  return {
    scale: (name: string) => {
      if (name !== "x") return undefined;
      return {
        invert: (px: number) => px / 100,
        apply:  (q: number) => q * 100,
      };
    },
  };
}

describe("invertQ", () => {
  it("returns the q value for the given pixel", () => {
    expect(invertQ(fakePlot(), 50)).toBe(0.5);
    expect(invertQ(fakePlot(), 200)).toBe(2);
  });

  it("returns null for a plot without a scale function", () => {
    expect(invertQ(undefined, 50)).toBeNull();
    expect(invertQ(null, 50)).toBeNull();
    expect(invertQ({}, 50)).toBeNull();
  });

  it("returns null when the x-scale lacks an invert method", () => {
    const broken = { scale: () => ({ apply: (q: number) => q }) };
    expect(invertQ(broken, 50)).toBeNull();
  });

  it("returns null for non-finite pixel coordinates", () => {
    expect(invertQ(fakePlot(), NaN)).toBeNull();
    expect(invertQ(fakePlot(), Infinity)).toBeNull();
  });
});

describe("applyQ", () => {
  it("returns the pixel for the given q", () => {
    expect(applyQ(fakePlot(), 0.5)).toBe(50);
  });

  it("returns null for a plot without a scale function", () => {
    expect(applyQ(undefined, 0.5)).toBeNull();
  });
});
