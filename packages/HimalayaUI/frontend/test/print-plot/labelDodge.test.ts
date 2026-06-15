import { describe, it, expect } from "vitest";
import { dodgeX } from "../../src/print/plot/marks/PlotLabels";
import { measureTextWidth } from "../../src/lib/plot/measureText";

describe("dodgeX (width-aware, minimal movement)", () => {
  it("leaves well-spaced labels exactly at their natural x (no spurious shift)", () => {
    // centers 0/50/100, half-widths 5, gap 5 → min centre distance 15; all clear.
    expect(dodgeX([0, 50, 100], [5, 5, 5], 5)).toEqual([0, 50, 100]);
  });

  it("pushes only the overlapping label, anchoring the leftmost", () => {
    // 0 and 8 are 8 apart but need >= 5+5+5 = 15: 0 stays, 8 → 15.
    const out = dodgeX([0, 8], [5, 5], 5);
    expect(out[0]).toBe(0);
    expect(out[1]).toBe(15);
  });

  it("leaves a clear neighbour untouched while resolving a tight pair", () => {
    expect(dodgeX([0, 4, 100], [5, 5, 5], 5)).toEqual([0, 15, 100]);
  });

  it("uses each label's own width (a narrow pair needs less room than a wide one)", () => {
    expect(dodgeX([0, 4], [2, 2], 1)[1]).toBe(5); // 0+2+2+1
    expect(dodgeX([0, 4], [10, 10], 1)[1]).toBe(21); // 0+10+10+1
  });

  it("returns [] for empty input", () => {
    expect(dodgeX([], [], 5)).toEqual([]);
  });
});

describe("measureTextWidth", () => {
  it("returns 0 for empty text", () => {
    expect(measureTextWidth("", { px: 11, weight: 700, family: "monospace" })).toBe(0);
  });

  it("falls back to a monospace estimate when canvas is unavailable (jsdom)", () => {
    // jsdom has no 2D canvas context, so the estimate path runs: len × px × 0.6.
    expect(measureTextWidth("ab", { px: 10, weight: 700, family: "monospace" })).toBeCloseTo(12, 5);
    // A wider label measures wider — this is what drives width-aware dodging.
    const w2 = measureTextWidth("√2", { px: 11, weight: 700, family: "monospace" });
    const w3 = measureTextWidth("√11", { px: 11, weight: 700, family: "monospace" });
    expect(w3).toBeGreaterThan(w2);
  });
});
