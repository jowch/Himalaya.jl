import { describe, it, expect } from "vitest";
import { buildTraceExportSpec } from "../../src/lib/figure-export/adapters/traceAdapter";
import { TRACE_DIMS } from "../../src/lib/figure-export/presets";
import type { Peak, IndexEntry, Trace } from "../../src/api";

const trace: Trace = {
  q: [0.1, 0.2, 0.3, 0.4, 0.5],
  I: [10, 8, 6, 4, 2],
  sigma: [1, 1, 1, 1, 1],
};
const peakAuto: Peak = {
  id: 1, exposure_id: 100, q: 0.18, intensity: 9, prominence: 5,
  sharpness: 1.0, source: "auto", excluded: false,
};
const peakManual: Peak = {
  ...peakAuto, id: 2, q: 0.32, source: "manual",
};
const peakExcluded: Peak = {
  ...peakAuto, id: 3, q: 0.42, excluded: true,
};
const indexPn3m: IndexEntry = {
  id: 10, exposure_id: 100, phase: "Pn3m", basis: 0.18, score: 0.9,
  r_squared: 0.99, lattice_d: 100, ngc: 1.5, status: "candidate", kind: "auto",
  inputs_hash: "h", peaks: [], predicted_q: [0.18, 0.36],
};

describe("buildTraceExportSpec", () => {
  it("produces an ExportSpec at TRACE_DIMS", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [peakAuto], activeGroupIndices: [indexPn3m],
      experimentName: "JC23", sampleName: "Sample 4", exposureLabel: "Exposure 7",
      xDomain: null, yDomain: null, xType: "log",
    });
    expect(spec.width).toBe(TRACE_DIMS.width);
    expect(spec.height).toBe(TRACE_DIMS.height);
  });

  it("primary title is 'experiment · sample · exposure'", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [],
      experimentName: "JC23", sampleName: "Sample 4", exposureLabel: "Exposure 7",
      xDomain: null, yDomain: null, xType: "log",
    });
    expect(spec.title.primary).toContain("JC23");
    expect(spec.title.primary).toContain("Sample 4");
    expect(spec.title.primary).toContain("Exposure 7");
  });

  it("legend includes 'auto' row when only auto peaks present", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [peakAuto], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label.toLowerCase());
    expect(labels.some((l) => l.includes("auto"))).toBe(true);
  });

  it("legend adds 'manual' row when manual peaks present", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [peakAuto, peakManual], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label.toLowerCase());
    expect(labels.some((l) => l.includes("manual"))).toBe(true);
  });

  it("legend adds 'excluded' row when excluded peaks present", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [peakAuto, peakExcluded], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label.toLowerCase());
    expect(labels.some((l) => l.includes("excluded"))).toBe(true);
  });

  it("legend adds one row per active-group index phase", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [indexPn3m],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label);
    expect(labels.some((l) => l.includes("Pn3m"))).toBe(true);
  });

  it("plot.x.type honours xType arg", () => {
    const specLog = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    const specLin = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "linear",
    });
    expect((specLog.plot.x as { type?: string }).type).toBe("log");
    expect((specLin.plot.x as { type?: string }).type).toBe("linear");
  });

  it("does NOT set title/caption/figure on plot (would break renderer)", () => {
    const spec = buildTraceExportSpec({
      trace, peaks: [], activeGroupIndices: [],
      experimentName: "E", sampleName: "S", exposureLabel: "X",
      xDomain: null, yDomain: null, xType: "log",
    });
    expect((spec.plot as { title?: unknown }).title).toBeUndefined();
    expect((spec.plot as { caption?: unknown }).caption).toBeUndefined();
    expect((spec.plot as { figure?: unknown }).figure).toBeUndefined();
  });
});
