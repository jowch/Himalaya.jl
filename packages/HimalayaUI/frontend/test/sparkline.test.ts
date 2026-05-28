import { describe, it, expect } from "vitest";
import { sparklinePath, SPARK_W, SPARK_H } from "../src/lib/plot/sparkline";

describe("sparklinePath", () => {
  it("returns null for an empty trace", () => {
    expect(sparklinePath({ q: [], I: [], sigma: [] })).toBeNull();
  });
  it("produces a path that starts with M and stays within the viewBox", () => {
    const q = Array.from({ length: 50 }, (_, i) => 0.02 + i * 0.008);
    const I = q.map((x) => Math.exp(-((x - 0.1) ** 2) / 0.0002));
    const p = sparklinePath({ q, I, sigma: q.map(() => 0) });
    expect(p).not.toBeNull();
    expect(p!.startsWith("M")).toBe(true);
    const ys = [...p!.matchAll(/[ML]\s*[\d.]+\s+([\d.]+)/g)].map((m) => Number(m[1]));
    expect(Math.min(...ys)).toBeGreaterThanOrEqual(0);
    expect(Math.max(...ys)).toBeLessThanOrEqual(SPARK_H);
  });
  it("ignores non-finite / non-positive intensities without throwing", () => {
    const p = sparklinePath({ q: [0.02, 0.04, 0.06], I: [NaN, 0, 5], sigma: [0, 0, 0] });
    expect(p).not.toBeNull();
  });
});

// Keep the geometry constants exported so consumers (ScopingSparkline) share them.
void SPARK_W;
