import { describe, it, expect } from "vitest";
import type { Peak, Trace } from "../../src/api";
import { toTraceModel } from "../../src/print/plot/traceModel";

function peak(over: Partial<Peak>): Peak {
  return {
    id: 1, exposure_id: 1, q: 0.1, intensity: 100, prominence: null,
    sharpness: null, source: "auto", excluded: false, ...over,
  };
}

describe("toTraceModel", () => {
  const trace: Trace = { q: [0.1, 0.2], I: [1, 2], sigma: [0, 0] };

  it("maps each Peak to a PlotPeak (id/q/source/excluded) and passes phase through", () => {
    const peaks = [
      peak({ id: 5, q: 0.11, intensity: 50, source: "auto", excluded: false }),
      peak({ id: 6, q: 0.22, intensity: null, source: "manual", excluded: true }),
    ];
    const tm = toTraceModel(trace, peaks, "Im3m");
    expect(tm.trace).toBe(trace);
    expect(tm.phase).toBe("Im3m");
    expect(tm.peaks).toHaveLength(2);
    expect(tm.peaks[0]).toMatchObject({ id: 5, q: 0.11, intensity: 50, source: "auto", excluded: false });
    expect(tm.peaks[1]).toMatchObject({ id: 6, q: 0.22, source: "manual", excluded: true });
  });

  it("passes a null phase through unchanged", () => {
    expect(toTraceModel(trace, [], null).phase).toBeNull();
  });

  it("anchors a null-intensity (curation) peak to the trace curve at its q (nearest sample)", () => {
    // trace = q:[0.1,0.2] I:[1,2]; a manually-added peak at q=0.1 has no backend
    // intensity → nearest sample is q=0.1 → I=1. Without this it would drop to
    // baselineI in PlotPeaks and render flat on the axis.
    const tm = toTraceModel(trace, [peak({ id: 1, q: 0.1, intensity: null })], null);
    expect(tm.peaks[0]).toMatchObject({ intensity: 1 });
  });

  it("omits intensity (does not set it to undefined) when there is no curve to anchor to", () => {
    const empty: Trace = { q: [], I: [], sigma: [] };
    const tm = toTraceModel(empty, [peak({ id: 1, intensity: null })], null);
    expect("intensity" in tm.peaks[0]!).toBe(false);
  });
});

