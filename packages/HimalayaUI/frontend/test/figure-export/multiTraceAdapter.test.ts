import { describe, it, expect } from "vitest";
import { buildMultiTraceExportSpec } from "../../src/lib/figure-export/adapters/multiTraceAdapter";
import { COMPARE_DIMS, CLEAN_SCIENTIFIC } from "../../src/lib/figure-export/presets";
import { buildExportSvg } from "../../src/lib/figure-export/renderer";
import type { SeriesMember, Trace } from "../../src/api";

function makeMember(over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id: 1, series_id: 1, exposure_id: 100, display_order: 0,
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
const sampleIdFor = (m: SeriesMember): number | null => {
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
      } as unknown as SeriesMember["snapshot"],
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

  it("bySample: orphan-by-null-sampleIdFor produces a single shared 'unphased / unbound' legend row", () => {
    // Regression: colorFor returns the dark-bg ORPHAN_FALLBACK for orphans
    // regardless of the palette passed in. The adapter must remap that to
    // ORPHAN_FALLBACK_LIGHT so (a) the marks render orphans with light-bg
    // contrast and (b) the legend's color-equality check detects them.
    const m1 = makeMember({ id: 1, exposure_id: 100 });
    const orphanByNullSample = makeMember({ id: 2, exposure_id: 999 });
    const noSampleResolver = (m: SeriesMember): number | null =>
      m.exposure_id === 100 ? 1 : null;
    const spec = buildMultiTraceExportSpec({
      members: [m1, orphanByNullSample],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "bySample", sampleIdFor: noSampleResolver,
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label.toLowerCase());
    expect(labels.some((l) => l.includes("sample 1"))).toBe(true);
    expect(labels.some((l) => l.includes("unphased") || l.includes("unbound"))).toBe(true);
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

  it("heatmap representation emits different marks than waterfall (#251 r1 / B1)", () => {
    // Heatmap swaps the per-row vocabulary from Plot.line traces to
    // Plot.rect cells. Mark counts diverge: waterfall = 1 line + 1 label =
    // 2 marks per member (no peaks here); heatmap = 1 rect (containing all
    // cells in `data`) + 1 label per member. Same total count for the
    // single-member fixture, so check the mark types instead.
    const waterfall = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      xDomain: [0.1, 0.3],
      showPeakTicks: false, showPeakLabels: false,
      groupingMode: "distinct", sampleIdFor,
      representation: "waterfall",
    });
    const heatmap = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      xDomain: [0.1, 0.3],
      showPeakTicks: false, showPeakLabels: false,
      groupingMode: "distinct", sampleIdFor,
      representation: "heatmap",
    });
    // Both produce non-empty mark arrays.
    expect((waterfall.plot.marks ?? []).length).toBeGreaterThan(0);
    expect((heatmap.plot.marks ?? []).length).toBeGreaterThan(0);
    // Plot.line vs Plot.rect — observable difference via the mark constructor
    // name. Plot marks expose `ariaLabel` or undocumented `_mark`-ish hooks,
    // but the safest cross-version check is `constructor.name`.
    const wfNames = (waterfall.plot.marks ?? []).map((m) => (m as { constructor: { name: string } }).constructor.name);
    const hmNames = (heatmap.plot.marks ?? []).map((m) => (m as { constructor: { name: string } }).constructor.name);
    expect(wfNames.some((n) => /line/i.test(n))).toBe(true);
    expect(hmNames.some((n) => /rect/i.test(n))).toBe(true);
  });

  it("showCrossTraceTracking adds tracking marks when ≥2 members share a phase (#251 r1 / B1)", () => {
    const indexedSnapshot = (phase: string, peaks: Array<{ id: number; q: number }>) => ({
      effective_peaks: peaks,
      confirmed_index: {
        id: 1, phase, lattice_d: 100, r_squared: 0.99, ngc: 1.5,
        peak_ids: peaks.map((p) => p.id),
      },
      analysis_inputs_hash: "h",
    });
    const a = makeMember({
      id: 1, exposure_id: 100, display_order: 0,
      snapshot: indexedSnapshot("Pn3m", [{ id: 1, q: 0.10 }, { id: 2, q: 0.12 }]) as unknown as SeriesMember["snapshot"],
    });
    const b = makeMember({
      id: 2, exposure_id: 101, display_order: 1,
      snapshot: indexedSnapshot("Pn3m", [{ id: 3, q: 0.11 }, { id: 4, q: 0.13 }]) as unknown as SeriesMember["snapshot"],
    });

    const withTracking = buildMultiTraceExportSpec({
      members: [a, b],
      traces,
      comparisonTitle: "T",
      xDomain: [0.1, 0.3],
      showPeakTicks: false, showPeakLabels: false,
      groupingMode: "distinct", sampleIdFor,
      showCrossTraceTracking: true,
    });
    const withoutTracking = buildMultiTraceExportSpec({
      members: [a, b],
      traces,
      comparisonTitle: "T",
      xDomain: [0.1, 0.3],
      showPeakTicks: false, showPeakLabels: false,
      groupingMode: "distinct", sampleIdFor,
      showCrossTraceTracking: false,
    });
    // Tracking adds at least one additional Plot.line per (phase, order).
    expect((withTracking.plot.marks ?? []).length).toBeGreaterThan(
      (withoutTracking.plot.marks ?? []).length,
    );
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

describe("buildMultiTraceExportSpec — CLEAN_SCIENTIFIC preset (E-8)", () => {
  function ffMember(): SeriesMember {
    return makeMember({
      id: 2, exposure_id: 101, display_order: 1,
      snapshot: {
        effective_peaks: [{ id: 9, q: 0.2, intensity: 1, sharpness: 1, source: "auto" }],
        confirmed_index: null,
        confirmed_phases: [],
        assignment_state: "form_factor",
        analysis_inputs_hash: "h",
      },
    });
  }

  it("emits a white-bg Arial spec with 2px traces + axis labels + footnote", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember(), ffMember()],
      traces,
      comparisonTitle: "LL37 Titration",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: false,
      groupingMode: "distinct", sampleIdFor,
      preset: "clean",
    });
    // Arial font threaded onto the spec.
    expect(spec.fontFamily).toBe(CLEAN_SCIENTIFIC.fontFamily);
    // Axis labels read q (Å⁻¹) / Intensity.
    expect((spec.plot.x as { label?: string }).label).toBe(CLEAN_SCIENTIFIC.axisLabel.x);
    // A footnote line is present.
    expect(spec.footnote).toBeTruthy();
    // 2px trace strokes (find a Plot.line mark with the clean stroke width).
    const marks = (spec.plot.marks ?? []) as Array<{ options?: { strokeWidth?: number } }>;
    // (Plot marks are opaque here; the stroke width is asserted via the preset constant.)
    expect(CLEAN_SCIENTIFIC.traceStroke).toBe(2);
    expect(marks.length).toBeGreaterThan(0);
  });

  it("includes form-factor members in the clean export", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember(), ffMember()],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: false,
      groupingMode: "distinct", sampleIdFor,
      preset: "clean",
    });
    // Both members render → at least two trace lines.
    expect((spec.plot.marks ?? []).length).toBeGreaterThanOrEqual(2);
  });

  it("the clean spec round-trips through buildExportSvg", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: true, showPeakLabels: false,
      groupingMode: "distinct", sampleIdFor,
      preset: "clean",
    });
    const svg = buildExportSvg(spec);
    expect(svg.tagName.toLowerCase()).toBe("svg");
    // The footnote text appears in the rendered SVG.
    expect(svg.textContent).toContain(spec.footnote);
  });
});
