import { describe, it, expect } from "vitest";
import { differenceAxis, noiseScale } from "../../src/workspace/difference";

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
    expect(differenceAxis("linear", q, d, [0.02, 0.22], []).kind).toBe("linear");
    expect(differenceAxis("kratky", q, d, [0.02, 0.22], []).kind).toBe("linear");
    const ps = differenceAxis("log", q, d, [0.02, 0.22], [0.101]);
    expect(ps.kind).toBe("pseudolog");
    expect(differenceAxis("loglog", q, d, [0.02, 0.22], [0.101]).kind).toBe("pseudolog");
    expect(ps.to(0)).toBe(0);
    expect(ps.to(-5)).toBeCloseTo(-ps.to(5));
    expect(ps.to(1000) / ps.to(100)).toBeLessThan(2);
  });

  it("keeps a peak far above the noise inside the range and cuts the low-q swing below the first peak", () => {
    const d = q.map((x) => (x < 0.04 ? -1e5 : x > 0.1 && x < 0.102 ? 1000 : Math.sin(x * 900)));
    const ax = differenceAxis("log", q, d, [0.02, 0.22], [0.101]);
    expect(ax.domain[1]).toBeGreaterThanOrEqual(ax.to(1000));
    expect(ax.domain[0]).toBeGreaterThan(ax.to(-1e5));
  });

  it("labels pseudolog ticks with real intensities, zero and powers of ten on both sides", () => {
    const d = q.map((x) => (x > 0.1 && x < 0.102 ? 1000 : x > 0.15 && x < 0.165 ? -300 : Math.sin(x * 900)));
    const ax = differenceAxis("log", q, d, [0.02, 0.22], [0.101]);
    const labels = ax.ticks.map((t) => t.label);
    expect(labels).toContain("0");
    expect(labels).toContain("1000");
    expect(labels.some((l) => l.startsWith("-"))).toBe(true);
    for (const t of ax.ticks) expect(t.value).toBeGreaterThanOrEqual(ax.domain[0]);
  });
});
