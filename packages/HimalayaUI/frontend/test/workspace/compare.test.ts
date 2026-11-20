import { describe, it, expect } from "vitest";
import { addCompare, compareScale, fluxRatio, traceFlux, transmittedFlux, TRACE_FRAME_AGGREGATION, clampDomain, stripRange, interpolateAt, matchScale, robustRange, scaledTrace, subtractTrace, MAX_COMPARES } from "../../src/workspace/compare";

const q = Array.from({ length: 101 }, (_, i) => 0.01 + i * 0.002);
const background = (x: number): number => 5 + 1 / x;
const bragg = (x: number): number => 400 * Math.exp(-((x - 0.1) ** 2) / 2e-5);

describe("compare", () => {
  it("interpolates linearly and is NaN outside the trace", () => {
    const t = { q: [0.1, 0.2], I: [10, 20] };
    expect(interpolateAt(t, 0.15)).toBeCloseTo(15);
    expect(interpolateAt(t, 0.1)).toBe(10);
    expect(interpolateAt(t, 0.25)).toBeNaN();
    expect(interpolateAt({ q: [], I: [] }, 0.1)).toBeNaN();
  });

  it("scales a buffer onto the sample's high-q background even when the sample has a peak", () => {
    const sample = { q, I: q.map((x) => 3 * background(x) + bragg(x)) };
    const buffer = { q, I: q.map(background) };
    const s = matchScale(sample, buffer);
    expect(s).toBeCloseTo(3, 3);
    const diff = subtractTrace(sample, buffer, s);
    const at = (x: number): number => diff.d[diff.q.findIndex((v) => Math.abs(v - x) < 1e-9)]!;
    expect(at(0.1)).toBeCloseTo(400, 0);
    expect(Math.abs(at(0.2))).toBeLessThan(1e-6);
  });

  it("falls back to a scale of 1 when the traces do not overlap", () => {
    expect(matchScale({ q: [0.1, 0.2], I: [1, 1] }, { q: [0.3, 0.4], I: [1, 1] })).toBe(1);
  });

  it("only subtracts where the reference covers the sample's q", () => {
    const diff = subtractTrace({ q: [0.1, 0.2, 0.3], I: [5, 5, 5] }, { q: [0.15, 0.35], I: [1, 1] }, 2);
    expect(diff).toEqual({ q: [0.2, 0.3], d: [3, 3] });
  });

  it("scales intensities without touching q", () => {
    expect(scaledTrace({ q: [0.1, 0.2], I: [1, 2] }, 3)).toEqual({ q: [0.1, 0.2], I: [3, 6] });
  });

  it("adds a comparison once and keeps the most recent few", () => {
    expect(addCompare([1, 2], 2)).toEqual([1, 2]);
    expect(addCompare([1, 2], 3)).toEqual([1, 2, 3]);
    const full = Array.from({ length: MAX_COMPARES }, (_, i) => i + 1);
    expect(addCompare(full, 9)).toEqual([...full.slice(1), 9]);
  });

  it("masks the sample's own peaks in the high-q window, which would otherwise over-scale the reference", () => {
    const bg = (x: number): number => 5 + 1 / x;
    const peaksAt = [0.19, 0.2];
    const sample = { q, I: q.map((x) => 1.3 * bg(x) + peaksAt.reduce((s, p) => s + 40 * Math.exp(-((x - p) ** 2) / 4e-6), 0)) };
    const buffer = { q, I: q.map(bg) };
    expect(matchScale(sample, buffer)).toBeGreaterThan(1.3 * 1.05);
    expect(matchScale(sample, buffer, peaksAt)).toBeCloseTo(1.3, 2);
  });

  it("takes the median, so one outlier in the window does not move the scale", () => {
    const buffer = { q, I: q.map(() => 1) };
    const sample = { q, I: q.map((x, i) => (i === q.length - 1 ? 1000 : 2)) };
    expect(matchScale(sample, buffer)).toBe(2);
  });

  it("sets the strip's range from the bulk of the difference, not one huge residual, and keeps zero in view", () => {
    const d = [...Array.from({ length: 200 }, (_, i) => (i % 20 === 0 ? 50 : 1)), -5000];
    const [lo, hi] = robustRange(d);
    expect(lo).toBe(0);
    expect(hi).toBe(50);
    expect(robustRange([3, 3, 3])).toEqual([0, 3]);
    expect(robustRange([])).toEqual([-1, 1]);
  });

  describe("stripRange with a realistic low-q residual", () => {
    const qs = Array.from({ length: 919 }, (_, i) => 0.0058 + i * ((0.4045 - 0.0058) / 918));
    const peakQ = [0.127, 0.18, 0.22];
    const d = qs.map(
      (x) => -4000 * Math.exp(-(x - 0.0058) / 0.008) + peakQ.reduce((a, p) => a + 20 * Math.exp(-((x - p) ** 2) / 4e-6), 0),
    );
    const residualShare = d.filter((v) => v < -20).length / d.length;

    it("the residual is more than the outlier cut can drop", () => {
      expect(residualShare).toBeGreaterThan(0.03);
      expect(robustRange(d)[0]).toBeLessThan(-100);
    });

    it("starts the range just below the first peak in view", () => {
      const [lo, hi] = stripRange(qs, d, [qs[0]!, qs[qs.length - 1]!], peakQ);
      expect(hi).toBeGreaterThan(15);
      expect(lo).toBeGreaterThan(-1);
    });

    it("skips the low-q end when no peak is in view", () => {
      const [lo] = stripRange(qs, d, [qs[0]!, 0.12], []);
      expect(lo).toBeGreaterThan(-20);
    });
  });

  it("clamps a kept zoom to the trace's q range and drops it when they do not overlap", () => {
    expect(clampDomain([0.1, 0.3], [0.05, 0.2])).toEqual([0.1, 0.2]);
    expect(clampDomain([0.3, 0.4], [0.05, 0.2])).toBeNull();
    expect(clampDomain(null, [0.05, 0.2])).toBeNull();
  });
});

describe("beam-monitor scaling", () => {
  const mon = (transmitted_counts: number | null, exposure_time: number | null = 15, counting_time: number | null = 15.8) => ({
    transmitted_counts, exposure_time, counting_time,
  });
  const withFrames = <T extends object>(t: T, ...frames: ReturnType<typeof mon>[]) => ({ ...t, frames });

  it("measures the flux through the sample during the detector exposure, not the longer counting window", () => {
    expect(transmittedFlux(mon(7900))).toBeCloseTo((7900 * 15) / 15.8, 6);
  });

  it("has no flux when a monitor is missing, zero or not a number", () => {
    expect(transmittedFlux(mon(null))).toBeNull();
    expect(transmittedFlux(mon(0))).toBeNull();
    expect(transmittedFlux(mon(7000, null))).toBeNull();
    expect(transmittedFlux(mon(7000, 15, null))).toBeNull();
    expect(transmittedFlux(mon(Number.NaN))).toBeNull();
    expect(transmittedFlux(undefined)).toBeNull();
  });

  it("adds every frame of a scan, a beam-off frame counting as no flux, and needs monitors on all of them", () => {
    expect(traceFlux([mon(7900), mon(7900)])).toEqual({ total: (2 * 7900 * 15) / 15.8, frames: 2 });
    expect(traceFlux([mon(7900), mon(-3)])).toEqual({ total: (7900 * 15) / 15.8, frames: 2 });
    expect(traceFlux([mon(7900), mon(null)])).toBeNull();
    expect(traceFlux([mon(-3)])).toBeNull();
    expect(traceFlux([])).toBeNull();
    expect(traceFlux(undefined)).toBeNull();
  });

  it("scales a raw-count reference by the ratio of transmitted flux, even when the sample scatters at high q", () => {
    const sample = withFrames({ q, I: q.map((x) => 1.2 * background(x) + 40 * x) }, mon(6000));
    const buffer = withFrames({ q, I: q.map((x) => 2 * background(x)) }, mon(10000));
    const r = compareScale(sample, buffer, []);
    expect(r.source).toBe("monitor");
    expect(r.scale).toBeCloseTo(0.6, 9);
    expect(Math.min(...subtractTrace(sample, buffer, r.scale).d)).toBeGreaterThan(0);
  });

  it("puts exposures of different length on one scale", () => {
    expect(compareScale(withFrames({ q, I: q }, mon(7000, 30, 30.8)), withFrames({ q, I: q }, mon(7000, 15, 15.8)), []).scale)
      .toBeCloseTo((30 / 30.8) / (15 / 15.8), 9);
  });

  it("uses the scan's total flux when both traces hold the same number of frames", () => {
    const r = compareScale(withFrames({ q, I: q }, mon(6000), mon(2000)), withFrames({ q, I: q }, mon(4000), mon(4000)), []);
    expect(r.source).toBe("monitor");
    expect(r.scale).toBeCloseTo(1, 9);
  });

  it("falls back to matching while it is unconfirmed whether a scan's frames are summed or averaged and the counts differ", () => {
    expect(TRACE_FRAME_AGGREGATION).toBe("unconfirmed");
    const sample = withFrames({ q, I: q.map((x) => 3 * background(x)) }, mon(7000), mon(7000));
    const buffer = withFrames({ q, I: q.map(background) }, mon(7000));
    expect(fluxRatio({ total: 2, frames: 2 }, { total: 1, frames: 1 })).toBeNull();
    const r = compareScale(sample, buffer, []);
    expect(r.source).toBe("matched");
    expect(r.scale).toBeCloseTo(3, 3);
  });

  it("falls back to the peak-masked high-q match when either trace lacks monitors", () => {
    const sample = withFrames({ q, I: q.map((x) => 3 * background(x)) }, mon(7000));
    const buffer = withFrames({ q, I: q.map(background) }, mon(null));
    const r = compareScale(sample, buffer, []);
    expect(r.source).toBe("matched");
    expect(r.scale).toBeCloseTo(3, 3);
    expect(compareScale({ q, I: sample.I }, withFrames({ q, I: buffer.I }, mon(7000)), []).source).toBe("matched");
  });
});
