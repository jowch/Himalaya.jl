import { describe, it, expect } from "vitest";
import { buildFocusFigureRow } from "../../src/print/export/focusFigureRow";
import type { Trace, Peak, IndexEntry } from "../../src/api";

const trace: Trace = {
  q: [0.1, 0.2, 0.3, 0.4],
  I: [10, 100, 50, 20],
  sigma: [1, 1, 1, 1],
};

function peak(id: number, q: number, intensity: number | null, extra?: Partial<Peak>): Peak {
  return {
    id,
    q,
    intensity,
    source: "auto",
    excluded: false,
    ...extra,
  } as Peak;
}

function index(id: number, phase: string, peakIds: number[]): IndexEntry {
  return {
    id,
    exposure_id: 1,
    phase,
    basis: 0.1,
    score: 0.9,
    r_squared: 0.99,
    lattice_d: 12,
    ngc: 2,
    status: "candidate",
    kind: "auto",
    inputs_hash: null,
    peaks: peakIds.map((pid) => ({ peak_id: pid, ratio_position: 1, residual: 0, q_observed: 0.2 })),
    predicted_q: [],
  } as IndexEntry;
}

describe("buildFocusFigureRow", () => {
  it("anchors only the peaks claimed by the active assignment, sorted by q", () => {
    const peaks = [peak(1, 0.2, 100), peak(2, 0.3, 50), peak(3, 0.4, 20)];
    const row = buildFocusFigureRow({
      trace,
      peaks,
      activeIndices: [index(10, "Im3m", [3, 1])], // claims peaks 1 and 3, NOT 2
      phase: "Im3m",
      label: "frame 1",
    });
    expect(row.anchors.map((a) => a.id)).toEqual([1, 3]); // ascending q
    expect(row.phase).toBe("Im3m");
    expect(row.anchors.every((a) => a.phase === "Im3m")).toBe(true);
  });

  it("places a claimed peak with no measured intensity onto the curve", () => {
    const peaks = [peak(1, 0.2, null)]; // manual/curation peak, no intensity
    const row = buildFocusFigureRow({
      trace,
      peaks,
      activeIndices: [index(10, "Im3m", [1])],
      phase: "Im3m",
      label: "frame 1",
    });
    // traceIntensityAt(0.2) on the trace = exact sample 100
    expect(row.anchors[0]!.intensity).toBe(100);
  });

  it("yields no anchors and a neutral row when there is no active assignment", () => {
    const peaks = [peak(1, 0.2, 100)];
    const row = buildFocusFigureRow({
      trace,
      peaks,
      activeIndices: [],
      phase: null,
      label: "frame 1",
    });
    expect(row.anchors).toEqual([]);
    expect(row.phase).toBeNull();
    expect(row.trace).toBe(trace);
    expect(row.bandHeight).toBe(1);
  });

  it("ignores claimed peak ids that are not in the observed peak set", () => {
    const peaks = [peak(1, 0.2, 100)];
    const row = buildFocusFigureRow({
      trace,
      peaks,
      activeIndices: [index(10, "Im3m", [1, 99])], // 99 doesn't exist
      phase: "Im3m",
      label: "frame 1",
    });
    expect(row.anchors.map((a) => a.id)).toEqual([1]);
  });
});
