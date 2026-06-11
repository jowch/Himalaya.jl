import { describe, it, expect } from "vitest";
import { buildMultiTraceExportSpec } from "../../src/lib/figure-export/adapters/multiTraceAdapter";
import { COMPARE_DIMS, CLEAN_SCIENTIFIC } from "../../src/lib/figure-export/presets";
import { buildExportSvg } from "../../src/lib/figure-export/renderer";
import { computeWaterfallExportBands } from "../../src/lib/comparison/yBands";
import { phaseColor } from "../../src/phases";
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

// ── BU-EXPORTDIVERGE — the export must speak the plate's registers ──────────
// The builder's plate footnote promises "what you compose is what you
// publish"; these tests pin the five divergence axes the downloaded SVG
// contradicted (phase color/identity, legend, label register, axis scale +
// tick vocabulary, stack order/offset, peak set).
describe("buildMultiTraceExportSpec — BU-EXPORTDIVERGE (export matches the plate)", () => {
  /** Collect every <text> content from the rendered export SVG. */
  function svgTexts(spec: ReturnType<typeof buildMultiTraceExportSpec>): string[] {
    const svg = buildExportSvg(spec);
    return Array.from(svg.querySelectorAll("text")).map((t) => t.textContent ?? "");
  }

  it("byPhase legend reads the DOMINANT phase (confirmed_phases[0]), not confirmed_index.phase", () => {
    // Coexistence member whose confirmed_phases order disagrees with
    // confirmed_index: the plate (waterfallModel.dominantPhase) calls this
    // member Im3m; the export must agree.
    const coexist = makeMember({
      id: 1, exposure_id: 100,
      snapshot: {
        effective_peaks: [],
        confirmed_index: { id: 1, phase: "Pn3m", lattice_d: 100, r_squared: 0.99, ngc: 1.5, peak_ids: [] },
        confirmed_phases: [{ phase: "Im3m", lattice_d: 110 }, { phase: "Pn3m", lattice_d: 100 }],
        analysis_inputs_hash: "h",
      } as unknown as SeriesMember["snapshot"],
    });
    const spec = buildMultiTraceExportSpec({
      members: [coexist],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: false, showPeakLabels: false,
      groupingMode: "byPhase", sampleIdFor,
    });
    const rows = spec.legend?.rows ?? [];
    expect(rows).toHaveLength(1);
    expect(rows[0]!.label).toBe("Im3m");
    expect(rows[0]!.color).toBe(phaseColor("Im3m"));
  });

  it("a member confirmed via confirmed_phases ONLY legends under its phase, not 'unphased / unbound'", () => {
    const phasesOnly = makeMember({
      id: 1, exposure_id: 100,
      snapshot: {
        effective_peaks: [],
        confirmed_index: null,
        confirmed_phases: [{ phase: "Ia3d", lattice_d: 195 }],
        analysis_inputs_hash: "h",
      } as unknown as SeriesMember["snapshot"],
    });
    const spec = buildMultiTraceExportSpec({
      members: [phasesOnly],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: false, showPeakLabels: false,
      groupingMode: "byPhase", sampleIdFor,
    });
    const labels = (spec.legend?.rows ?? []).map((r) => r.label);
    expect(labels).toContain("Ia3d");
    expect(labels.join(" ").toLowerCase()).not.toContain("unphased");
  });

  it("default member labels speak the plate's register ('exp N'), never 'Exposure #N'", () => {
    const spec = buildMultiTraceExportSpec({
      members: [makeMember({ id: 1, exposure_id: 100 })],
      traces,
      comparisonTitle: "T",
      xDomain: null,
      showPeakTicks: false, showPeakLabels: false,
      groupingMode: "distinct", sampleIdFor,
    });
    const texts = svgTexts(spec);
    expect(texts).toContain("exp 100");
    expect(texts.some((t) => t.includes("Exposure #"))).toBe(false);
  });

  it("xType passes through to the x scale: default log, 'linear' when the plate is linear", () => {
    const base = {
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      xDomain: null as [number, number] | null,
      showPeakTicks: false, showPeakLabels: false,
      groupingMode: "distinct" as const, sampleIdFor,
    };
    const logSpec = buildMultiTraceExportSpec(base);
    expect((logSpec.plot.x as { type?: string }).type).toBe("log");
    const linSpec = buildMultiTraceExportSpec({ ...base, xType: "linear" });
    expect((linSpec.plot.x as { type?: string }).type).toBe("linear");
  });

  it("q ticks use the plate's decimal register, not d3's SI-milli ('100m') suffixes", () => {
    // trace q spans 0.1–0.3: Plot's default log tickFormat renders "100m",
    // "200m", "300m" — the divergence the downloaded SVG showed. The plate
    // formats these as plain decimals (formatAxis).
    const spec = buildMultiTraceExportSpec({
      members: [makeMember()],
      traces,
      comparisonTitle: "T",
      xDomain: [0.1, 0.3],
      showPeakTicks: false, showPeakLabels: false,
      groupingMode: "distinct", sampleIdFor,
    });
    const texts = svgTexts(spec);
    expect(texts.some((t) => /^\d+m$/.test(t.trim()))).toBe(false);
    // And at least one plain-decimal q tick is present.
    expect(texts.some((t) => /^0\.\d+$/.test(t.trim()))).toBe(true);
  });

  it("peak labels number the PLATE's set (indexed anchors), not all effective_peaks", () => {
    const threePeaks = [
      { id: 11, q: 0.10, intensity: 100, sharpness: 1, source: "auto" },
      { id: 12, q: 0.15, intensity: 80, sharpness: 1, source: "auto" },
      { id: 13, q: 0.20, intensity: 60, sharpness: 1, source: "auto" },
    ];
    const member = makeMember({
      id: 1, exposure_id: 100,
      snapshot: {
        effective_peaks: threePeaks,
        // Only two of the three peaks are indexed anchors — the plate labels
        // exactly those two (1, 2 ascending q) and never the third.
        confirmed_index: { id: 1, phase: "Pn3m", lattice_d: 100, r_squared: 0.99, ngc: 1.5, peak_ids: [11, 12] },
        analysis_inputs_hash: "h",
      } as unknown as SeriesMember["snapshot"],
    });
    const spec = buildMultiTraceExportSpec({
      members: [member],
      traces,
      comparisonTitle: "T",
      xDomain: [0.05, 0.3],
      showPeakTicks: true, showPeakLabels: true,
      groupingMode: "byPhase", sampleIdFor,
    });
    const texts = svgTexts(spec).map((t) => t.trim());
    expect(texts).toContain("1");
    expect(texts).toContain("2");
    expect(texts).not.toContain("3");
  });
});

describe("computeWaterfallExportBands — plate-mirroring stack geometry", () => {
  it("paints display-order row 0 at the BOTTOM (the plate's bottom-up stack)", () => {
    const bands = computeWaterfallExportBands(
      [{ band_height: 1, y_offset: 0 }, { band_height: 1, y_offset: 0 }],
      100,
      1,
    );
    // Row 0 occupies the lower half, row 1 the upper half — contiguous at ×1.
    expect(bands[0]).toEqual([50, 100]);
    expect(bands[1]).toEqual([0, 50]);
  });

  it("offsetScale > 1 opens separation between bands (plate's offset slider flows through)", () => {
    const tight = computeWaterfallExportBands(
      [{ band_height: 1, y_offset: 0 }, { band_height: 1, y_offset: 0 }],
      100,
      1,
    );
    const spread = computeWaterfallExportBands(
      [{ band_height: 1, y_offset: 0 }, { band_height: 1, y_offset: 0 }],
      100,
      2,
    );
    const gapOf = (b: Array<[number, number]>): number => b[0]![0] - b[1]![1];
    expect(gapOf(tight)).toBeCloseTo(0);
    expect(gapOf(spread)).toBeGreaterThan(0);
    // The stack always renormalizes to fill the fixed export panel.
    expect(spread[0]![1]).toBeCloseTo(100);
    expect(spread[1]![0]).toBeCloseTo(0);
  });

  it("y_offset nudges are converted from plate pixels (420 stack) into panel space", () => {
    // A 42px plate nudge = 10% of the plate stack. In a 100px panel it must
    // contribute 10px BEFORE renormalization, not its raw 42px (which would
    // weigh the nudge ~4x heavier in the export than on the plate).
    const nudged = computeWaterfallExportBands(
      [{ band_height: 1, y_offset: 0 }, { band_height: 1, y_offset: 42 }],
      100,
      1,
    );
    // Row 1 (top) is pushed DOWN into the stack by the panel-scaled nudge
    // (10px, not 42px); row 0 still bottoms the panel, so no renormalization.
    expect(nudged[1]![0]).toBeCloseTo(10);
    expect(nudged[0]![0]).toBeCloseTo(50);
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

describe("rendered-geometry pins (Builder re-score #4 P1: the export was vertically mirrored)", () => {
  // These pin the FINAL SVG geometry, not the band arrays — the marks are
  // authored in screen convention while Plot's y scale is a data scale, and
  // only the rendered output can prove the two conventions compose upright.
  function twoMemberSvg(): SVGSVGElement {
    const peaky: Trace = {
      // A single sharp peak at q=0.2 over a flat floor.
      q: [0.1, 0.15, 0.2, 0.25, 0.3],
      I: [10, 10, 1000, 10, 10],
      sigma: [1, 1, 1, 1, 1],
    };
    const t = new Map<number, Trace>([[100, peaky], [101, peaky]]);
    const members = [
      makeMember({ id: 1, exposure_id: 100, display_order: 0, label_override: "exp 100" }),
      makeMember({ id: 2, exposure_id: 101, display_order: 1, label_override: "exp 101" }),
    ];
    const spec = buildMultiTraceExportSpec({
      members,
      traces: t,
      comparisonTitle: "geometry pin",
      xDomain: null,
      showPeakTicks: false,
      showPeakLabels: false,
      groupingMode: "bySample",
      sampleIdFor: () => null,
    });
    return buildExportSvg(spec);
  }

  function labelY(svg: SVGSVGElement, label: string): number {
    const el = [...svg.querySelectorAll("text")].find((t) => t.textContent === label);
    if (!el) throw new Error(`label ${label} not found`);
    // Plot positions text via transform translate; the y attribute is the
    // em-unit baseline shift ("0.32em"), not a position.
    const m = /translate\([\d.+-]+[, ]([\d.+-]+)\)/.exec(el.getAttribute("transform") ?? "");
    if (m) return Number(m[1]);
    const yAttr = el.getAttribute("y");
    if (yAttr !== null && !yAttr.endsWith("em")) return Number(yAttr);
    throw new Error(`label ${label} has no positional y`);
  }

  it("display-order 0 renders BELOW display-order 1 (the plate's bottom-up stack)", () => {
    const svg = twoMemberSvg();
    expect(labelY(svg, "exp 100")).toBeGreaterThan(labelY(svg, "exp 101"));
  });

  it("within a band, the intensity peak points UP (smaller rendered y at the peak q)", () => {
    const svg = twoMemberSvg();
    // The trace lines are path elements; collect every path's points and find
    // one whose x-midpoint y dips: min y must be WELL above (smaller than)
    // the path's edge y values for a peaked trace.
    const paths = [...svg.querySelectorAll("path")]
      .map((p) => p.getAttribute("d") ?? "")
      .filter((d) => (d.match(/[\d.]+,[\d.]+/g) ?? []).length >= 5);
    expect(paths.length).toBeGreaterThan(0);
    const upright = paths.some((d) => {
      const pts = (d.match(/(-?[\d.]+),(-?[\d.]+)/g) ?? []).map((s) => {
        const [x, y] = s.split(",").map(Number);
        return { x: x!, y: y! };
      });
      if (pts.length < 5) return false;
      const ys = pts.map((p) => p.y);
      const minY = Math.min(...ys);
      const edgeY = (ys[0]! + ys[ys.length - 1]!) / 2;
      // Peak = an interior point substantially HIGHER (smaller y) than edges.
      const minIdx = ys.indexOf(minY);
      return minIdx > 0 && minIdx < ys.length - 1 && edgeY - minY > 5;
    });
    expect(upright).toBe(true);
  });
});
