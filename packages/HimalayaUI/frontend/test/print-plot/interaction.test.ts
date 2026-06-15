import { describe, it, expect } from "vitest";
import {
  hitTestPeaks,
  zoomXDomain,
  PEAK_HIT_PX,
} from "../../src/print/plot/interaction";

const toPx = (q: number) => q * 100; // q=1 → 100px

describe("hitTestPeaks", () => {
  it("returns the nearest peak within tolerance", () => {
    const peaks = [
      { id: 1, q: 1 },
      { id: 2, q: 2 },
    ];
    expect(hitTestPeaks(peaks, 105, toPx, PEAK_HIT_PX)?.id).toBe(1);
  });

  it("returns null beyond tolerance", () => {
    const peaks = [{ id: 1, q: 1 }];
    expect(hitTestPeaks(peaks, 130, toPx, PEAK_HIT_PX)).toBeNull();
  });

  it("skips optimistic placeholders (id < 0)", () => {
    const peaks = [{ id: -1, q: 1 }];
    expect(hitTestPeaks(peaks, 100, toPx, PEAK_HIT_PX)).toBeNull();
  });

  it("the later peak wins an exact tie (<=)", () => {
    const peaks = [
      { id: 1, q: 1 },
      { id: 2, q: 1 },
    ];
    expect(hitTestPeaks(peaks, 100, toPx, PEAK_HIT_PX)?.id).toBe(2);
  });
});

describe("zoomXDomain — log", () => {
  it("zooms in about the cursor, narrowing the domain", () => {
    const next = zoomXDomain({
      cursorQ: 0.1,
      deltaY: -100,
      current: [0.01, 1],
      extent: [0.01, 1],
      type: "log",
    });
    expect(next).not.toBeNull();
    expect(next![0]).toBeGreaterThanOrEqual(0.01);
    expect(next![1]).toBeLessThanOrEqual(1);
    expect(next![1] - next![0]).toBeLessThan(1 - 0.01);
  });

  it("clamps to the full extent when zooming out", () => {
    const next = zoomXDomain({
      cursorQ: 0.1,
      deltaY: 1000,
      current: [0.05, 0.2],
      extent: [0.01, 1],
      type: "log",
    });
    expect(next![0]).toBeGreaterThanOrEqual(0.01);
    expect(next![1]).toBeLessThanOrEqual(1);
  });

  it("returns null when the span would get too small", () => {
    const tiny = zoomXDomain({
      cursorQ: 0.1,
      deltaY: -100000,
      current: [0.0999, 0.1001],
      extent: [0.01, 1],
      type: "log",
    });
    expect(tiny).toBeNull();
  });
});

describe("zoomXDomain — linear", () => {
  it("narrows about the cursor in linear space", () => {
    const next = zoomXDomain({
      cursorQ: 5,
      deltaY: -100,
      current: [0, 10],
      extent: [0, 10],
      type: "linear",
    });
    expect(next![1] - next![0]).toBeLessThan(10);
  });
});
