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
import type { SeriesMember } from "../src/api";
import { buildMemberMarks } from "../src/components/MemberTraceLayer";
import { phaseColor } from "../src/phases";
import {
  COMPARE_PALETTE,
  ORPHAN_FALLBACK,
  colorFor,
} from "../src/lib/comparison/coloring";

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

// ── Phase 9 gap-fix: line stroke wired to grouping mode ────────────────────
//
// Until this batch, `MemberTraceLayer` rendered the line stroke as
// `member.color_override ?? "var(--color-fg)"`. The grouping toggle changed
// the picker palette + (already-implemented) hover phase coloring, but the
// line color stayed the same regardless of mode. These tests pin the wired
// behaviour: the stroke is `colorFor(member, mode, palette, ctx)`.
function getLineStroke(): string | undefined {
  const lineCalls = (Plot.line as unknown as { mock: { calls: unknown[][] } }).mock.calls;
  const last = lineCalls[lineCalls.length - 1];
  if (!last) return undefined;
  const opts = last[1] as { stroke?: string } | undefined;
  return opts?.stroke;
}

describe("buildMemberMarks line stroke (Phase 9 grouping-mode wiring)", () => {
  const trivialResolver = () => null;

  it("bySample: two members with the same sample_id render the same stroke", () => {
    const a = makeMember({ id: 1, exposure_id: 10, display_order: 0, snapshot: null });
    const b = makeMember({ id: 2, exposure_id: 11, display_order: 1, snapshot: null });
    const sampleIdFor = (m: SeriesMember) =>
      m.exposure_id === 10 ? 7 : m.exposure_id === 11 ? 7 : null;
    const ctx = { allMembers: [a, b], sampleIdFor };

    buildMemberMarks({
      member: a,
      trace,
      yBand: [0, 100],
      groupingMode: "bySample",
      allMembers: [a, b],
      sampleIdFor,
    });
    const strokeA = getLineStroke();
    (Plot.line as unknown as { mockClear: () => void }).mockClear();

    buildMemberMarks({
      member: b,
      trace,
      yBand: [0, 100],
      groupingMode: "bySample",
      allMembers: [a, b],
      sampleIdFor,
    });
    const strokeB = getLineStroke();

    expect(strokeA).toBeDefined();
    expect(strokeA).toBe(strokeB);
    // And both match the colorFor library's resolution exactly.
    expect(strokeA).toBe(colorFor(a, "bySample", COMPARE_PALETTE, ctx));
  });

  it("byPhase: two members with the same confirmed_index.phase render the same stroke", () => {
    const a = makeMember({ id: 1, exposure_id: 10, display_order: 0 });
    const b = makeMember({ id: 2, exposure_id: 11, display_order: 1 });

    buildMemberMarks({
      member: a,
      trace,
      yBand: [0, 100],
      groupingMode: "byPhase",
      allMembers: [a, b],
      sampleIdFor: trivialResolver,
    });
    const strokeA = getLineStroke();
    (Plot.line as unknown as { mockClear: () => void }).mockClear();

    buildMemberMarks({
      member: b,
      trace,
      yBand: [0, 100],
      groupingMode: "byPhase",
      allMembers: [a, b],
      sampleIdFor: trivialResolver,
    });
    const strokeB = getLineStroke();

    expect(strokeA).toBe(strokeB);
    expect(strokeA).toBe(phaseColor("Pn3m"));
  });

  it("distinct: each member gets a unique stroke", () => {
    const a = makeMember({ id: 1, exposure_id: 10, display_order: 0 });
    const b = makeMember({ id: 2, exposure_id: 11, display_order: 1 });

    buildMemberMarks({
      member: a,
      trace,
      yBand: [0, 100],
      groupingMode: "distinct",
      allMembers: [a, b],
      sampleIdFor: trivialResolver,
    });
    const strokeA = getLineStroke();
    (Plot.line as unknown as { mockClear: () => void }).mockClear();

    buildMemberMarks({
      member: b,
      trace,
      yBand: [0, 100],
      groupingMode: "distinct",
      allMembers: [a, b],
      sampleIdFor: trivialResolver,
    });
    const strokeB = getLineStroke();

    expect(strokeA).not.toBe(strokeB);
  });

  it("color_override always wins regardless of grouping mode", () => {
    const m = makeMember({ color_override: "#abcdef" });
    for (const mode of ["bySample", "byPhase", "distinct"] as const) {
      (Plot.line as unknown as { mockClear: () => void }).mockClear();
      buildMemberMarks({
        member: m,
        trace,
        yBand: [0, 100],
        groupingMode: mode,
        allMembers: [m],
        sampleIdFor: () => 42,
      });
      expect(getLineStroke()).toBe("#abcdef");
    }
  });

  it("suppresses peak anchors for a form-factor member (real trace, no Bragg marks) — E-7", () => {
    const m = makeMember({
      snapshot: {
        effective_peaks: [
          { id: 11, q: 0.30, intensity: 50, sharpness: 1, source: "auto" },
        ],
        confirmed_index: null,
        confirmed_phases: [],
        assignment_state: "form_factor",
        analysis_inputs_hash: "abc",
      },
    });
    buildMemberMarks({ member: m, trace, yBand: [0, 100] });
    // The trace line still renders…
    expect(Plot.line).toHaveBeenCalled();
    // …but NO peak dot/glyph marks (no anchors on a form-factor trace).
    expect(Plot.dot).not.toHaveBeenCalled();
    // …and a form-factor member keeps a FULL-opacity trace (structured
    // broad-shoulder scattering is real signal worth reading at full strength).
    const ffOpts = (Plot.line as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![1] as { strokeOpacity?: number };
    expect(ffOpts.strokeOpacity).toBe(1);
  });

  it("suppresses peak anchors for a null member (featureless) — E-7", () => {
    const m = makeMember({
      snapshot: {
        effective_peaks: [
          { id: 11, q: 0.30, intensity: 50, sharpness: 1, source: "auto" },
        ],
        confirmed_index: null,
        confirmed_phases: [],
        assignment_state: "null",
        analysis_inputs_hash: "abc",
      },
    });
    buildMemberMarks({ member: m, trace, yBand: [0, 100] });
    expect(Plot.line).toHaveBeenCalled();
    expect(Plot.dot).not.toHaveBeenCalled();
    // A null member ("nothing interesting") reads DE-EMPHASIZED — a dimmed real
    // trace, visually distinct from a full-opacity form_factor member and from
    // an indexed member. This is the three-state distinction ON the waterfall,
    // not just in the phase-strip / member rows.
    const nullOpts = (Plot.line as unknown as { mock: { calls: unknown[][] } })
      .mock.calls[0]![1] as { strokeOpacity?: number };
    expect(nullOpts.strokeOpacity).toBe(0.4);
  });

  it("orphan member (no confirmed_index) under byPhase falls back to the orphan gray", () => {
    const m = makeMember({ snapshot: null });
    buildMemberMarks({
      member: m,
      trace,
      yBand: [0, 100],
      groupingMode: "byPhase",
      allMembers: [m],
      sampleIdFor: () => null,
    });
    expect(getLineStroke()).toBe(ORPHAN_FALLBACK);
  });
});
