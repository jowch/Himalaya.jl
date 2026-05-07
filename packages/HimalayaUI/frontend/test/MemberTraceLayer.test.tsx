/**
 * MemberTraceLayer tests (Plan §Phase 6, Task 6.1).
 *
 * MemberTraceLayer is a *mark factory* — it returns an array of Observable
 * Plot marks for one comparison member, not its own SVG. The component
 * doesn't render visible DOM under JSDOM; we test the data flow by reading
 * the data passed to mocked `Plot.line` / `Plot.dot` / `Plot.text` factories.
 *
 * `mockReturnValueOnce` etc. is avoided — we let the mocks accumulate calls
 * across one render and inspect the args by call-index, since the mark
 * factory invokes them in a deterministic order.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import * as Plot from "@observablehq/plot";
import type { ComparisonMember } from "../src/api";
import { buildMemberMarks } from "../src/components/MemberTraceLayer";
import { phaseColor } from "../src/phases";

vi.mock("@observablehq/plot", () => ({
  line: vi.fn((data: unknown, opts: unknown) => ({ _kind: "line", data, opts })),
  dot:  vi.fn((data: unknown, opts: unknown) => ({ _kind: "dot",  data, opts })),
  text: vi.fn((data: unknown, opts: unknown) => ({ _kind: "text", data, opts })),
}));

beforeEach(() => {
  (Plot.line as unknown as { mockClear: () => void }).mockClear();
  (Plot.dot  as unknown as { mockClear: () => void }).mockClear();
  (Plot.text as unknown as { mockClear: () => void }).mockClear();
});

// Convenience builders.
function makeMember(over: Partial<ComparisonMember> = {}): ComparisonMember {
  return {
    id: 1,
    comparison_id: 100,
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
      effective_peaks: [
        { id: 11, q: 0.30, intensity: 50, sharpness: 1, source: "auto" },
        { id: 12, q: 0.50, intensity: 80, sharpness: 1, source: "auto" },
      ],
      confirmed_index: {
        id: 7, phase: "Pn3m", lattice_d: 12, r_squared: 0.99, ngc: -1.51,
        peak_ids: [11, 12],
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
  q: [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7],
  I: [10, 20, 30, 50, 80, 60, 40],
  sigma: [0, 0, 0, 0, 0, 0, 0],
};

describe("buildMemberMarks (Phase 6.1)", () => {
  it("emits a line mark with one point per trace sample", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
    });
    expect(Plot.line).toHaveBeenCalledTimes(1);
    const lineCall = (Plot.line as unknown as { mock: { calls: unknown[][] } }).mock.calls[0]!;
    const lineData = lineCall[0] as Array<{ q: number; y: number }>;
    expect(lineData).toHaveLength(trace.q.length);
    // q values preserved on the line data.
    expect(lineData[0]!.q).toBe(0.1);
    expect(lineData[lineData.length - 1]!.q).toBe(0.7);
  });

  it("emits a dot mark with one entry per visible (not hidden) snapshot peak", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
    });
    const dotCalls = (Plot.dot as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    expect(dotCalls.length).toBeGreaterThanOrEqual(1);
    // Aggregate all peak data emitted (we may split black vs phase-colored).
    const peakRows: Array<{ q: number; peakId: number }> = [];
    for (const call of dotCalls) {
      const data = call[0] as Array<{ q: number; peakId: number }>;
      for (const row of data) peakRows.push(row);
    }
    expect(peakRows.map((r) => r.q).sort()).toEqual([0.30, 0.50]);
  });

  it("omits hidden peaks from the dot mark", () => {
    const m = makeMember();
    buildMemberMarks({
      member: m,
      trace,
      yBand: [0, 100],
      peakDisplay: { hidden: [11], labeled: [] },
    });
    const dotCalls = (Plot.dot as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    const peakIds: number[] = [];
    for (const call of dotCalls) {
      for (const row of call[0] as Array<{ peakId: number }>) peakIds.push(row.peakId);
    }
    expect(peakIds).not.toContain(11);
    expect(peakIds).toContain(12);
  });

  it("emits text labels only for peaks listed in peak_display.labeled", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
      peakDisplay: { hidden: [], labeled: [12] },
    });
    expect(Plot.text).toHaveBeenCalled();
    const textData = (Plot.text as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![0] as Array<{ q: number; peakId: number }>;
    expect(textData).toHaveLength(1);
    expect(textData[0]!.peakId).toBe(12);
  });

  it("does not emit a text mark when no labels are configured", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
      peakDisplay: { hidden: [], labeled: [] },
    });
    expect(Plot.text).not.toHaveBeenCalled();
  });

  it("renders peaks BLACK by default (no highlight)", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
    });
    // Each Plot.dot call carries an `opts.fill` (or per-row color); we accept
    // either a fill at the mark level or per-row data tagged with color.
    const dotCalls = (Plot.dot as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    let sawBlack = false;
    let sawNonBlack = false;
    for (const call of dotCalls) {
      const opts = call[1] as { fill?: string | ((d: unknown) => string) } | undefined;
      const data = call[0] as Array<{ color?: string }>;
      const fillIsString = typeof opts?.fill === "string";
      const blackOpt = fillIsString && (opts!.fill === "black" || opts!.fill === "#000" || opts!.fill === "#000000");
      if (blackOpt) sawBlack = true;
      for (const row of data) {
        if (row.color) {
          if (row.color === "black" || row.color === "#000" || row.color === "#000000") sawBlack = true;
          else sawNonBlack = true;
        }
      }
    }
    expect(sawBlack).toBe(true);
    expect(sawNonBlack).toBe(false);
  });

  it("recolors snapshot.confirmed_index peaks to phase color when highlightedIndexId matches", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
      highlightedIndexId: 7,
    });
    const dotCalls = (Plot.dot as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    const expected = phaseColor("Pn3m");
    let sawExpected = false;
    for (const call of dotCalls) {
      const opts = call[1] as { fill?: string | ((d: unknown) => string) } | undefined;
      const data = call[0] as Array<{ peakId: number; color?: string }>;
      // either the mark-level fill is the phase color, or per-row color is.
      if (typeof opts?.fill === "string" && opts.fill === expected) sawExpected = true;
      for (const row of data) {
        if (row.color === expected && (row.peakId === 11 || row.peakId === 12)) {
          sawExpected = true;
        }
      }
    }
    expect(sawExpected).toBe(true);
  });

  it("does NOT recolor when highlightedIndexId does not match snapshot.confirmed_index.id", () => {
    buildMemberMarks({
      member: makeMember(),
      trace,
      yBand: [0, 100],
      highlightedIndexId: 999,
    });
    const dotCalls = (Plot.dot as unknown as { mock: { calls: unknown[][] } }).mock.calls;
    const expected = phaseColor("Pn3m");
    let sawExpected = false;
    for (const call of dotCalls) {
      const opts = call[1] as { fill?: string } | undefined;
      const data = call[0] as Array<{ color?: string }>;
      if (opts?.fill === expected) sawExpected = true;
      for (const row of data) if (row.color === expected) sawExpected = true;
    }
    expect(sawExpected).toBe(false);
  });

  it("produces no peak/label marks when snapshot is null (orphan member)", () => {
    buildMemberMarks({
      member: makeMember({ snapshot: null }),
      trace,
      yBand: [0, 100],
    });
    expect(Plot.dot).not.toHaveBeenCalled();
    expect(Plot.text).not.toHaveBeenCalled();
    // Line mark may still be emitted if a trace was passed.
  });

  it("produces no marks at all when trace is undefined", () => {
    buildMemberMarks({
      member: makeMember(),
      trace: undefined,
      yBand: [0, 100],
    });
    expect(Plot.line).not.toHaveBeenCalled();
    // Peaks/labels still depend on the trace for y-positioning, so they are
    // also suppressed.
    expect(Plot.dot).not.toHaveBeenCalled();
    expect(Plot.text).not.toHaveBeenCalled();
  });
});
