import { describe, it, expect } from "vitest";
import {
  makeAxis,
  makeProjection,
  positiveExtent,
} from "../../src/print/plot/projection";

describe("makeAxis", () => {
  it("linear to/invert round-trips", () => {
    const a = makeAxis([0, 10], [0, 100], "linear");
    expect(a.to(5)).toBeCloseTo(50);
    expect(a.invert(50)).toBeCloseTo(5);
  });

  it("log maps decades evenly and inverts", () => {
    const a = makeAxis([0.01, 1], [0, 200], "log");
    expect(a.to(0.01)).toBeCloseTo(0);
    expect(a.to(1)).toBeCloseTo(200);
    expect(a.to(0.1)).toBeCloseTo(100);
    expect(a.invert(100)).toBeCloseTo(0.1);
  });

  it("ticks returns values inside the domain", () => {
    const a = makeAxis([0.01, 1], [0, 200], "log");
    const t = a.ticks(5);
    expect(t.length).toBeGreaterThan(0);
    for (const v of t) {
      expect(v).toBeGreaterThanOrEqual(0.01);
      expect(v).toBeLessThanOrEqual(1);
    }
  });
});

describe("makeProjection", () => {
  it("inverts the y range so the domain max sits at the top (px 0)", () => {
    const p = makeProjection({
      xDomain: [0.01, 1],
      yDomain: [1, 1000],
      plotWidth: 300,
      plotHeight: 200,
      xType: "log",
      yType: "log",
    });
    expect(p.y.to(1000)).toBeCloseTo(0);
    expect(p.y.to(1)).toBeCloseTo(200);
    expect(p.x.to(0.01)).toBeCloseTo(0);
    expect(p.x.to(1)).toBeCloseTo(300);
  });
});

describe("positiveExtent", () => {
  it("ignores non-finite and non-positive values", () => {
    expect(positiveExtent([0, -1, 2, NaN, 8])).toEqual([2, 8]);
  });

  it("falls back when there is no positive data", () => {
    expect(positiveExtent([0, -5, NaN])).toEqual([1, 10]);
  });
});
