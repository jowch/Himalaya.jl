/**
 * MemberHeatmapLayer tests (#208 — render-core heatmap representation).
 *
 * Mark factories invoke mocked `Plot.rect` and we inspect the captured
 * args. Two layers of contract:
 *
 *   - `buildMemberHeatmapMarks` shape: one `Plot.rect` mark per member, with
 *     `HEATMAP_BIN_COUNT` cells, each carrying a `color-mix` fill that lerps
 *     the resolved group/phase hue against the paper plate by intensity.
 *
 *   - `binPeaksOnly`: pure helper, peak-on-decay should bin to a tall single
 *     bin (the peak) and zero (or near-zero) elsewhere after the rolling-min
 *     baseline subtraction.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import * as Plot from "@observablehq/plot";
import type { SeriesMember } from "../src/api";
import {
  buildMemberHeatmapMarks,
  binPeaksOnly,
  HEATMAP_BIN_COUNT,
} from "../src/components/MemberHeatmapLayer";
import { phaseColor } from "../src/phases";
import { COMPARE_PALETTE } from "../src/lib/comparison/coloring";

vi.mock("@observablehq/plot", () => ({
  rect: vi.fn((data: unknown, opts: unknown) => ({ _kind: "rect", data, opts })),
}));

beforeEach(() => {
  (Plot.rect as unknown as { mockClear: () => void }).mockClear();
});

function makeMember(over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id: 1,
    series_id: 100,
    exposure_id: 42,
    display_order: 0,
    band_height: 1,
    y_offset: 0,
    normalization: "qwindow",
    color_override: null,
    label_override: null,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot: {
      effective_peaks: [],
      confirmed_index: {
        id: 7, phase: "Pn3m", lattice_d: 12, r_squared: 0.99, ngc: -1.51,
        peak_ids: [],
      },
      analysis_inputs_hash: "abc123",
    },
    is_stale: false,
    created_by: null,
    created_at: null,
    ...over,
  };
}

const trace = {
  q: [0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9],
  I: [10, 20, 30, 50, 80, 60, 40, 30, 20, 10],
  sigma: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
};

describe("buildMemberHeatmapMarks", () => {
  it("emits no marks when the trace is missing or empty", () => {
    expect(buildMemberHeatmapMarks({
      member: makeMember(),
      trace: undefined,
      yBand: [0, 100],
      qDomain: [0, 1],
    })).toEqual([]);
    expect(buildMemberHeatmapMarks({
      member: makeMember(),
      trace: { q: [], I: [], sigma: [] },
      yBand: [0, 100],
      qDomain: [0, 1],
    })).toEqual([]);
  });

  it("emits cells + an outer keyline rect (HEATMAP_BIN_COUNT cells across the q-domain)", () => {
    buildMemberHeatmapMarks({
      member: makeMember(),
      trace,
      yBand: [0, 50],
      qDomain: [0, 1],
    });
    // Two rect marks now: the binned cells + the framing keyline (R3-Y06).
    expect(Plot.rect).toHaveBeenCalledTimes(2);
    const calls = (Plot.rect as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    // The cells call is the one WITHOUT `fill:"none"` (the keyline has it).
    const cellsCall = calls.find(
      (c) => (c[1] as { fill?: unknown }).fill !== "none",
    )!;
    const cells = cellsCall[0] as Array<{ x1: number; x2: number; fill: string; memberId: number }>;
    expect(cells).toHaveLength(HEATMAP_BIN_COUNT);
    // First cell starts at qDomain[0], last cell ends at qDomain[1].
    expect(cells[0]!.x1).toBeCloseTo(0, 6);
    expect(cells[HEATMAP_BIN_COUNT - 1]!.x2).toBeCloseTo(1, 6);
    // Each cell carries a `color-mix(...)` fill string.
    for (const c of cells) {
      expect(c.fill).toMatch(/^color-mix\(in oklab,/);
    }
    // memberId tag survives so a downstream hover layer can identify the row.
    expect(new Set(cells.map((c) => c.memberId))).toEqual(new Set([1]));
  });

  it("emits an outer hair keyline rect framing each row (R3-Y06)", () => {
    const marks = buildMemberHeatmapMarks({
      member: makeMember(),
      trace,
      yBand: [10, 30],
      qDomain: [0.05, 0.9],
      groupingMode: "byPhase",
      allMembers: [makeMember()],
      sampleIdFor: () => 1,
    });
    // Two rect marks: the binned cells + the framing keyline.
    expect(marks.length).toBe(2);
    const calls = (Plot.rect as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    // The keyline is the call with fill:none + a hair stroke.
    const keyline = calls.find(
      (c) => (c[1] as { fill?: unknown }).fill === "none",
    );
    expect(keyline).toBeDefined();
    const opts = keyline![1] as { stroke: string; strokeWidth: number };
    expect(opts.stroke).toBe("var(--color-hair)");
    expect(opts.strokeWidth).toBe(1);
  });

  it("colours each cell against the resolved group/phase hue", () => {
    // byPhase + an indexed member → phase hue (`Pn3m`).
    buildMemberHeatmapMarks({
      member: makeMember(),
      trace,
      yBand: [0, 50],
      qDomain: [0, 1],
      groupingMode: "byPhase",
      allMembers: [makeMember()],
      sampleIdFor: () => 1,
    });
    const cells = ((Plot.rect as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![0] as Array<{ fill: string }>);
    // The fill string embeds the phase colour as its "ink" side.
    const pn3m = phaseColor("Pn3m");
    expect(cells.every((c) => c.fill.includes(pn3m))).toBe(true);
  });

  it("colours rows via COMPARE_PALETTE in bySample mode (paper-tuned hue)", () => {
    buildMemberHeatmapMarks({
      member: makeMember(),
      trace,
      yBand: [0, 50],
      qDomain: [0, 1],
      groupingMode: "bySample",
      allMembers: [makeMember()],
      sampleIdFor: () => 7,
    });
    const cells = ((Plot.rect as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![0] as Array<{ fill: string }>);
    // Some palette entry is used as the ink — every cell's fill must
    // reference that exact entry.
    const matched = COMPARE_PALETTE.some((p) =>
      cells.every((c) => c.fill.includes(p))
    );
    expect(matched).toBe(true);
  });

  it("y-band envelope is inset (rows don't fuse)", () => {
    buildMemberHeatmapMarks({
      member: makeMember(),
      trace,
      yBand: [0, 50],
      qDomain: [0, 1],
    });
    const opts = (Plot.rect as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![1] as { y1: number; y2: number };
    expect(opts.y1).toBeGreaterThan(0);
    expect(opts.y2).toBeLessThan(50);
    expect(opts.y2).toBeGreaterThan(opts.y1);
  });
});

describe("binPeaksOnly", () => {
  it("returns a zero-array when the trace is empty", () => {
    const bins = binPeaksOnly([], [], 0, 1, 10);
    expect(bins).toEqual(new Array(10).fill(0));
  });

  it("a single sharp peak survives baseline subtraction", () => {
    // Build a synthetic peak-on-decay: gentle decay 0.05→0.95 q, single
    // sharp spike at q=0.5 with amplitude ~10× the local decay value.
    const qs: number[] = [];
    const Is: number[] = [];
    for (let i = 0; i < 200; i++) {
      const q = 0.05 + (i / 199) * 0.9;
      qs.push(q);
      const decay = 0.5 / Math.max(q, 0.01);
      const peak = Math.exp(-((q - 0.5) ** 2) / (2 * 0.005 ** 2)) * 30;
      Is.push(decay + peak);
    }
    const bins = binPeaksOnly(qs, Is, 0, 1, 50);
    // The peak bin should be substantially above the global mean.
    const peakBin = Math.floor((0.5 - 0) / (1 - 0) * 50);
    const mean = bins.reduce((a, b) => a + b, 0) / bins.length;
    expect(bins[peakBin]!).toBeGreaterThan(mean * 5);
  });
});
