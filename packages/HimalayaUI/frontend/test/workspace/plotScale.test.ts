import { describe, it, expect } from "vitest";
import { decadeLabel, fitYDomain, linearTicks, logTicks } from "../../src/workspace/plotScale";

describe("plotScale", () => {
  it("writes decades as powers of ten", () => {
    expect(decadeLabel(0)).toBe("10⁰");
    expect(decadeLabel(-12)).toBe("10⁻¹²");
    expect(decadeLabel(5)).toBe("10⁵");
  });

  it("labels a log intensity axis by decade when two or more decades are in view", () => {
    expect(logTicks([0.5, 2e3], 5, true).map((t) => t.label)).toEqual(["10⁰", "10¹", "10²", "10³"]);
  });

  it("thins the decades to about the tick count", () => {
    const ticks = logTicks([1e-2, 1e9], 5, true);
    expect(ticks.length).toBeLessThanOrEqual(6);
    expect(ticks[0]!.label).toBe("10⁻²");
  });

  it("falls back to plain labels inside one decade, and always for a log q axis", () => {
    expect(logTicks([20, 80], 5, true).every((t) => !t.label.startsWith("10") || t.value === 10)).toBe(true);
    expect(logTicks([20, 80], 5, true).length).toBeGreaterThanOrEqual(2);
    expect(logTicks([0.005, 0.4], 8, false).map((t) => t.label)).toContain("0.10");
  });

  it("gives linear ticks exact labels", () => {
    expect(linearTicks([0, 0.3], 8).map((t) => t.label)).toContain("0.10");
  });

  it("pads the fitted range more at the top than the bottom", () => {
    const series = { q: [0.1, 0.2, 0.3], I: [10, 100, 1000] };
    const [lo, hi] = fitYDomain(series, [], [0.1, 0.3], "log");
    expect(Math.log10(10) - Math.log10(lo)).toBeCloseTo(0.12);
    expect(Math.log10(hi) - Math.log10(1000)).toBeCloseTo(0.192);
    const [llo, lhi] = fitYDomain(series, [], [0.1, 0.3], "linear");
    expect(10 - llo).toBeCloseTo(59.4);
    expect(lhi - 1000).toBeCloseTo(95.04);
  });
});
