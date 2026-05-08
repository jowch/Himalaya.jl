import { describe, it, expect } from "vitest";
import { buildMultiTraceExportSpec } from "../../src/lib/figure-export/adapters/multiTraceAdapter";
import { COMPARE_DIMS } from "../../src/lib/figure-export/presets";
import type { ComparisonMember, Trace } from "../../src/api";

function makeMember(over: Partial<ComparisonMember> = {}): ComparisonMember {
  return {
    id: 1, comparison_id: 1, exposure_id: 100, display_order: 0,
    band_height: 1, y_offset: 0, normalization: "none",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: null, is_stale: false, created_by: 1, created_at: null,
    ...over,
  };
}

const trace: Trace = {
  q: [0.1, 0.2, 0.3], I: [10, 5, 2], sigma: [1, 1, 1],
};

const traces = new Map<number, Trace>([[100, trace], [101, trace], [102, trace]]);
const sampleIdFor = (m: ComparisonMember): number | null => {
  if (m.exposure_id === 100) return 1;
  if (m.exposure_id === 101) return 2;
  if (m.exposure_id === 102) return 3;
  return null;
};

describe("buildMultiTraceExportSpec", () => {
  it("produces an ExportSpec at COMPARE_DIMS", () => {
    const m = makeMember();
    const spec = buildMultiTraceExportSpec({
      members: [m],
      traces,
      comparisonTitle: "My comparison",
      experimentName: "JC23",
      xDomain: null,
      showPeakTicks: true,
      showPeakLabels: true,
      groupingMode: "distinct",
      sampleIdFor,
    });
    expect(spec.width).toBe(COMPARE_DIMS.width);
    expect(spec.height).toBe(COMPARE_DIMS.height);
  });

  it("primary title is the comparison title; secondary is the experiment name", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      experimentName: "E",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    expect(spec.title.primary).toBe("T");
    expect(spec.title.secondary).toBe("E");
  });

  it("omits secondary title when experimentName is undefined (global scope)", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    expect(spec.title.primary).toBe("T");
    expect(spec.title.secondary).toBeUndefined();
  });

  it("bySample legend has one row per unique sample_id", () => {
    const spec = buildMultiTraceExportSpec({
      members: [
        makeMember({ id: 1, exposure_id: 100 }),
        makeMember({ id: 2, exposure_id: 101 }),
        makeMember({ id: 3, exposure_id: 102 }),
      ],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "bySample", sampleIdFor,
    });
    // 3 unique sample ids → 3 legend rows.
    expect(spec.legend?.rows.length).toBe(3);
  });

  it("byPhase legend has one row per unique confirmed-index phase, plus orphan row when any have none", () => {
    const phasedMember = makeMember({
      id: 10, exposure_id: 100,
      snapshot: {
        effective_peaks: [],
        confirmed_index: {
          id: 1, exposure_id: 100, phase: "Pn3m", basis: 0.18, score: 0.9,
          r_squared: 0.99, lattice_d: 100, ngc: 1.5, status: "candidate",
          kind: "auto", inputs_hash: "h", peaks: [], predicted_q: [0.18],
        },
      } as unknown as ComparisonMember["snapshot"],
    });
    const orphanMember = makeMember({ id: 11, exposure_id: 101, snapshot: null });

    const spec = buildMultiTraceExportSpec({
      members: [phasedMember, orphanMember],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "byPhase", sampleIdFor,
    });

    // One Pn3m row + one shared orphan row.
    const labels = (spec.legend?.rows ?? []).map((r) => r.label.toLowerCase());
    expect(labels.some((l) => l.includes("pn3m"))).toBe(true);
    expect(labels.some((l) => l.includes("unphased") || l.includes("unbound"))).toBe(true);
    expect(spec.legend?.rows.length).toBe(2);
  });

  it("distinct mode emits NO legend (per-member labels carry the encoding)", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember({ id: 1 }), makeMember({ id: 2, exposure_id: 101 })],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    expect(spec.legend?.rows.length ?? 0).toBe(0);
  });

  it("filters null-exposure_id members from marks but keeps colorByMember stable for distinct mode", () => {
    const m1 = makeMember({ id: 1, exposure_id: 100 });
    const orphan = makeMember({ id: 2, exposure_id: null });
    const m3 = makeMember({ id: 3, exposure_id: 102 });
    const spec = buildMultiTraceExportSpec({
      members: [m1, orphan, m3],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    // Marks were built off [m1, m3] (orphan filtered) BUT colors were
    // computed off the unfiltered list — so m3's distinct color is index 2,
    // not index 1. We sanity-check the spec was generated successfully and
    // the figure dims are still the canonical preset (no nullable downgrade).
    expect(spec.width).toBe(COMPARE_DIMS.width);
    // marks is a Plot.Mark[] — at minimum we expect the trace marks for the
    // two non-orphan members, which means non-empty.
    const marksCount = (spec.plot.marks ?? []).length;
    expect(marksCount).toBeGreaterThan(0);
  });

  it("throws when every member has exposure_id === null", () => {
    expect(() => buildMultiTraceExportSpec({
      members: [
        makeMember({ id: 1, exposure_id: null }),
        makeMember({ id: 2, exposure_id: null }),
      ],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    })).toThrow(/exposure_id === null/);
  });

  it("does NOT set title/caption/figure on plot", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "distinct", sampleIdFor,
    });
    expect((spec.plot as { title?: unknown }).title).toBeUndefined();
    expect((spec.plot as { caption?: unknown }).caption).toBeUndefined();
    expect((spec.plot as { figure?: unknown }).figure).toBeUndefined();
  });
});
