import { describe, it, expect } from "vitest";
import type { Peak } from "../../src/api";
import { traceGeometry } from "../../src/workspace/plot";

const trace = { q: [0.05, 0.1, 0.15, 0.2], I: [100, 10, 1000, 1], sigma: [0, 0, 0, 0] };
const peak = (id: number, q: number, excluded = false): Peak =>
  ({ id, exposure_id: 1, q, intensity: null, prominence: null, sharpness: null, source: "auto", excluded }) as Peak;

describe("traceGeometry", () => {
  it("returns null without a usable trace", () => {
    expect(traceGeometry(undefined, [], [], 400, 300)).toBeNull();
    expect(traceGeometry({ q: [0.1], I: [1], sigma: [0] }, [], [], 400, 300)).toBeNull();
  });

  it("plots log intensity so the brightest point is highest", () => {
    const g = traceGeometry(trace, [], [], 400, 300)!;
    const ys = g.path.map(([, y]) => y);
    expect(Math.min(...ys)).toBe(ys[2]);
    expect(Math.max(...ys)).toBe(ys[3]);
    expect(g.path[0]![0]).toBe(g.plot.left);
    expect(g.path[3]![0]).toBe(g.plot.right);
  });

  it("drops non-positive intensities", () => {
    const g = traceGeometry({ q: [0.05, 0.1, 0.15], I: [10, 0, 100], sigma: [0, 0, 0] }, [], [], 400, 300)!;
    expect(g.path).toHaveLength(2);
  });

  it("places peaks on the curve and keeps their excluded flag", () => {
    const g = traceGeometry(trace, [peak(1, 0.15), peak(2, 0.1, true), peak(3, 0.5)], [], 400, 300)!;
    expect(g.peaks.map((p) => p.id)).toEqual([1, 2]);
    expect(g.peaks[0]!.y).toBeCloseTo(g.path[2]![1]);
    expect(g.peaks[1]!.excluded).toBe(true);
  });

  it("keeps only comb lines inside the q range", () => {
    const g = traceGeometry(trace, [], [{ id: 7, qs: [0.01, 0.1, 0.3], color: "red", emphasis: true }], 400, 300)!;
    expect(g.combs[0]!.xs).toHaveLength(1);
  });

  it("labels x ticks inside the q range", () => {
    const g = traceGeometry(trace, [], [], 400, 300)!;
    expect(g.xTicks.length).toBeGreaterThan(1);
    for (const t of g.xTicks) {
      expect(t.x).toBeGreaterThanOrEqual(g.plot.left - 0.01);
      expect(t.x).toBeLessThanOrEqual(g.plot.right + 0.01);
    }
  });
});
