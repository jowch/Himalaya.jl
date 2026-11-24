import { describe, it, expect } from "vitest";
import { differenceAxis, noiseScale, subtractedSeries } from "../../src/workspace/difference";

const q = Array.from({ length: 201 }, (_, i) => 0.02 + i * 0.001);

describe("difference axis", () => {
  it("takes the linear band from the noise in the upper half of the q range", () => {
    const d = q.map((x, i) => (x > 0.12 ? [-2, 0, 2][i % 3]! : 500));
    expect(noiseScale(q, d, [0.02, 0.22])).toBeCloseTo(2 * 1.4826, 3);
  });

  it("never returns a zero band", () => {
    expect(noiseScale(q, q.map(() => 0), [0.02, 0.22])).toBeGreaterThan(0);
  });

  it("is linear in Linear and q²I, and pseudolog in the log views", () => {
    const d = q.map((x) => (x > 0.1 && x < 0.102 ? 1000 : Math.sin(x * 900)));
    expect(differenceAxis("linear", q, d, [0.02, 0.22]).kind).toBe("linear");
    expect(differenceAxis("kratky", q, d, [0.02, 0.22]).kind).toBe("linear");
    const ps = differenceAxis("log", q, d, [0.02, 0.22]);
    expect(ps.kind).toBe("pseudolog");
    expect(differenceAxis("loglog", q, d, [0.02, 0.22]).kind).toBe("pseudolog");
    expect(ps.to(0)).toBe(0);
    expect(ps.to(-5)).toBeCloseTo(-ps.to(5));
    expect(ps.to(1000) / ps.to(100)).toBeLessThan(2);
  });

  it("holds the whole difference in view, the low-q swing below the first peak included", () => {
    const d = q.map((x) => (x < 0.04 ? -1e5 : x > 0.1 && x < 0.102 ? 1000 : Math.sin(x * 900)));
    for (const view of ["log", "loglog", "linear", "kratky"] as const) {
      const ax = differenceAxis(view, q, d, [0.02, 0.22]);
      expect(ax.domain[1]).toBeGreaterThanOrEqual(ax.to(1000));
      expect(ax.domain[0]).toBeLessThanOrEqual(ax.to(-1e5));
    }
  });

  it("holds a lone dip that is under 1% of the points", () => {
    const d = q.map((x) => (x > 0.15 && x < 0.1505 ? -400 : Math.sin(x * 900)));
    const ax = differenceAxis("linear", q, d, [0.02, 0.22]);
    expect(ax.domain[0]).toBeLessThanOrEqual(-400);
  });

  it("fits only the q range in view, so a zoom past the low-q swing gives the peaks the height", () => {
    const d = q.map((x) => (x < 0.04 ? -1e5 : x > 0.1 && x < 0.102 ? 1000 : Math.sin(x * 900)));
    const ax = differenceAxis("linear", q, d, [0.05, 0.22]);
    expect(ax.domain[0]).toBeGreaterThan(-200);
    expect(ax.domain[1]).toBeGreaterThanOrEqual(1000);
  });

  it("labels pseudolog ticks with real intensities, zero and powers of ten on both sides", () => {
    const d = q.map((x) => (x > 0.1 && x < 0.102 ? 1000 : x > 0.15 && x < 0.165 ? -300 : Math.sin(x * 900)));
    const ax = differenceAxis("log", q, d, [0.02, 0.22]);
    const labels = ax.ticks.map((t) => t.label);
    expect(labels).toContain("0");
    expect(labels).toContain("1000");
    expect(labels.some((l) => l.startsWith("-"))).toBe(true);
    for (const t of ax.ticks) expect(t.value).toBeGreaterThanOrEqual(ax.domain[0]);
  });
});

describe("subtracted curve series", () => {
  const sq = [0.02, 0.05, 0.1, 0.15, 0.2];
  const sd = [5000, 300, 20, -3, 0.5];

  it("is on a logarithmic y in Log and Log-log and a linear y in Linear and Kratky", () => {
    expect(subtractedSeries("log", sq, sd, [0.02, 0.2]).yType).toBe("log");
    expect(subtractedSeries("loglog", sq, sd, [0.02, 0.2]).yType).toBe("log");
    expect(subtractedSeries("linear", sq, sd, [0.02, 0.2]).yType).toBe("linear");
    expect(subtractedSeries("kratky", sq, sd, [0.02, 0.2]).yType).toBe("linear");
  });

  it("holds the whole curve in range, choosing the log floor from the positive values", () => {
    const [lo, hi] = subtractedSeries("log", sq, sd, [0.02, 0.2]).domain;
    expect(lo).toBeGreaterThan(0);
    expect(lo).toBeLessThanOrEqual(0.5);
    expect(lo).toBeGreaterThan(0.01);
    expect(hi).toBeGreaterThanOrEqual(5000);
  });

  it("holds the whole curve in a linear view, from the lowest value to the largest", () => {
    const [lo, hi] = subtractedSeries("linear", sq, sd, [0.02, 0.2]).domain;
    expect(lo).toBeLessThanOrEqual(-3);
    expect(hi).toBeGreaterThanOrEqual(5000);
  });

  it("multiplies by q² in Kratky and leaves the values alone otherwise", () => {
    expect(subtractedSeries("kratky", sq, sd, [0.02, 0.2]).I[2]).toBeCloseTo(0.1 * 0.1 * 20);
    expect(subtractedSeries("log", sq, sd, [0.02, 0.2]).I).toEqual(sd);
  });

  it("views other curves the same way and keeps them in the y range", () => {
    const other = { q: sq, d: [9e5, 4e4, 20, 3, 1] };
    const out = subtractedSeries("kratky", sq, sd, [0.02, 0.2], [other]);
    expect(out.others[0]![2]).toBeCloseTo(0.1 * 0.1 * 20);
    expect(out.domain[1]).toBeGreaterThan(subtractedSeries("kratky", sq, sd, [0.02, 0.2]).domain[1]);
  });
});
