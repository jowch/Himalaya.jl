/**
 * labelDodge tests (Plan §Phase 8, Task 8.2).
 *
 * The dodge layout function takes a sorted list of labeled-peak positions
 * (in q + pixel-y space) plus a width estimate and emits dodged label
 * positions + per-label leader-line endpoints. Output guarantees:
 *
 *   - Sparse peaks: label q stays exactly equal to peak q (no leader).
 *   - Crowded peaks: labels spread horizontally to avoid collisions; the
 *     dodged label q ≠ peak q ⇒ a leader line is needed.
 *   - Order preserved: the i-th input peak maps to the i-th output entry.
 *
 * The function operates in pixel space for the dodge math (so widths are
 * comparable) and accepts a q→pixel + pixel→q pair for converting back to
 * the data coordinate system. Tests use the identity scale (px = q*100)
 * and the inverse (q = px/100) so the assertions are easy to read.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import * as Plot from "@observablehq/plot";
import { layoutPeakLabels } from "../src/lib/plot/labelDodge";
import { buildMemberMarks } from "../src/components/MemberTraceLayer";

vi.mock("@observablehq/plot", () => ({
  line: vi.fn((data: unknown, opts: unknown) => ({ _kind: "line", data, opts })),
  dot:  vi.fn((data: unknown, opts: unknown) => ({ _kind: "dot",  data, opts })),
  text: vi.fn((data: unknown, opts: unknown) => ({ _kind: "text", data, opts })),
  link: vi.fn((data: unknown, opts: unknown) => ({ _kind: "link", data, opts })),
}));

beforeEach(() => {
  (Plot.line as unknown as { mockClear: () => void }).mockClear();
  (Plot.dot  as unknown as { mockClear: () => void }).mockClear();
  (Plot.text as unknown as { mockClear: () => void }).mockClear();
  (Plot.link as unknown as { mockClear: () => void }).mockClear();
});

const toPx = (q: number) => q * 100;
const fromPx = (px: number) => px / 100;

interface Peak { q: number; y: number; peakId: number }

function makePeaks(qs: number[]): Peak[] {
  return qs.map((q, i) => ({ q, y: 50, peakId: i + 1 }));
}

describe("layoutPeakLabels (Phase 8.2)", () => {
  it("sparse peaks: labels sit directly above each triangle (no dodge)", () => {
    const peaks = makePeaks([0.10, 0.50, 0.90]);
    const out = layoutPeakLabels(peaks, { toPx, fromPx, labelWidthPx: 30 });
    expect(out).toHaveLength(3);
    for (const o of out) {
      expect(o.qLabel).toBe(o.qPeak);
    }
  });

  it("crowded peaks: dodge spreads labels so no two overlap in pixel space", () => {
    // Three peaks all within 20 pixels of each other (q=0.30, 0.32, 0.34
    // → px=30, 32, 34) should not all be able to share a 30 px wide label.
    const peaks = makePeaks([0.30, 0.32, 0.34]);
    const out = layoutPeakLabels(peaks, { toPx, fromPx, labelWidthPx: 30 });
    expect(out).toHaveLength(3);

    // Verify pixel-space non-overlap pairwise.
    const xs = out.map((d) => toPx(d.qLabel)).sort((a, b) => a - b);
    for (let i = 1; i < xs.length; i++) {
      const gap = xs[i]! - xs[i - 1]!;
      expect(gap).toBeGreaterThanOrEqual(30);
    }
  });

  it("emits leader lines (qLabel ≠ qPeak) for displaced labels", () => {
    const peaks = makePeaks([0.30, 0.32, 0.34]);
    const out = layoutPeakLabels(peaks, { toPx, fromPx, labelWidthPx: 30 });
    // At least one (and likely all) labels must have moved off-center.
    const moved = out.filter((d) => d.qLabel !== d.qPeak);
    expect(moved.length).toBeGreaterThan(0);
  });

  it("preserves input order (i-th input → i-th output)", () => {
    const peaks: Peak[] = [
      { q: 0.10, y: 50, peakId: 100 },
      { q: 0.30, y: 50, peakId: 200 },
      { q: 0.50, y: 50, peakId: 300 },
    ];
    const out = layoutPeakLabels(peaks, { toPx, fromPx, labelWidthPx: 30 });
    expect(out.map((d) => d.peakId)).toEqual([100, 200, 300]);
  });

  it("leader-line endpoints connect to the original triangle position", () => {
    const peaks = makePeaks([0.30, 0.32]);
    const out = layoutPeakLabels(peaks, { toPx, fromPx, labelWidthPx: 30 });
    for (let i = 0; i < out.length; i++) {
      expect(out[i]!.qPeak).toBe(peaks[i]!.q);
      expect(out[i]!.yPeak).toBe(peaks[i]!.y);
    }
  });

  it("includes a label string derived from the peak q-value", () => {
    const peaks = makePeaks([0.123456]);
    const out = layoutPeakLabels(peaks, { toPx, fromPx, labelWidthPx: 30 });
    expect(out[0]!.label.length).toBeGreaterThan(0);
    // q=0.123456.toPrecision(3) = "0.123" — current format.
    expect(out[0]!.label).toBe("0.123");
  });

  it("places label y above the peak y by a stable offset", () => {
    const peaks: Peak[] = [{ q: 0.30, y: 100, peakId: 1 }];
    const out = layoutPeakLabels(peaks, { toPx, fromPx, labelWidthPx: 30 });
    expect(out[0]!.yLabel).toBeLessThan(100);
  });

  it("returns empty array for empty input", () => {
    const out = layoutPeakLabels([], { toPx, fromPx, labelWidthPx: 30 });
    expect(out).toEqual([]);
  });

  it("handles a single peak (returns it unchanged)", () => {
    const peaks: Peak[] = [{ q: 0.42, y: 60, peakId: 9 }];
    const out = layoutPeakLabels(peaks, { toPx, fromPx, labelWidthPx: 30 });
    expect(out).toHaveLength(1);
    expect(out[0]!.qLabel).toBe(0.42);
    expect(out[0]!.qPeak).toBe(0.42);
  });

  it("integration: buildMemberMarks emits Plot.link marks for crowded labels with xScale", () => {
    const member = {
      id: 1, series_id: 1, exposure_id: 42, display_order: 0,
      band_height: 1, y_offset: 0, normalization: "qwindow" as const,
      color_override: null, label_override: null,
      q_window_min: null, q_window_max: null, peak_display: null,
      snapshot: {
        effective_peaks: [
          { id: 11, q: 0.30, intensity: 50, sharpness: 1, source: "auto" as const },
          { id: 12, q: 0.31, intensity: 50, sharpness: 1, source: "auto" as const },
          { id: 13, q: 0.32, intensity: 50, sharpness: 1, source: "auto" as const },
        ],
        confirmed_index: null,
        analysis_inputs_hash: "h",
      },
      is_stale: false, created_by: null, created_at: null,
    };
    const trace = {
      q: [0.10, 0.20, 0.30, 0.31, 0.32, 0.40],
      I: [1, 2, 3, 4, 5, 6],
      sigma: [0, 0, 0, 0, 0, 0],
    };

    buildMemberMarks({
      member,
      trace,
      yBand: [0, 100],
      peakDisplay: { hidden: [], labeled: [11, 12, 13] },
      xScale: { toPx, fromPx },
    });

    // Crowded labels at q=0.30,0.31,0.32 (px=30,31,32) all want a 32-px-wide
    // slot, so the dodge MUST push them apart and emit at least one link.
    expect((Plot.text as unknown as { mock: { calls: unknown[] } }).mock.calls.length).toBeGreaterThan(0);
    expect((Plot.link as unknown as { mock: { calls: unknown[] } }).mock.calls.length).toBeGreaterThan(0);
  });

  it("integration: buildMemberMarks does NOT emit a link mark for sparse labels", () => {
    const member = {
      id: 1, series_id: 1, exposure_id: 42, display_order: 0,
      band_height: 1, y_offset: 0, normalization: "qwindow" as const,
      color_override: null, label_override: null,
      q_window_min: null, q_window_max: null, peak_display: null,
      snapshot: {
        effective_peaks: [
          { id: 11, q: 0.10, intensity: 50, sharpness: 1, source: "auto" as const },
          { id: 12, q: 0.50, intensity: 50, sharpness: 1, source: "auto" as const },
          { id: 13, q: 0.90, intensity: 50, sharpness: 1, source: "auto" as const },
        ],
        confirmed_index: null,
        analysis_inputs_hash: "h",
      },
      is_stale: false, created_by: null, created_at: null,
    };
    const trace = { q: [0.10, 0.50, 0.90], I: [1, 2, 3], sigma: [0, 0, 0] };

    buildMemberMarks({
      member,
      trace,
      yBand: [0, 100],
      peakDisplay: { hidden: [], labeled: [11, 12, 13] },
      xScale: { toPx, fromPx },
    });

    expect((Plot.text as unknown as { mock: { calls: unknown[] } }).mock.calls.length).toBeGreaterThan(0);
    expect((Plot.link as unknown as { mock: { calls: unknown[] } }).mock.calls.length).toBe(0);
  });

  it("dodge spreads labels symmetrically around the cluster center", () => {
    // Two peaks crowded at q=0.40, 0.42 should one shift left, one right.
    const peaks = makePeaks([0.40, 0.42]);
    const out = layoutPeakLabels(peaks, { toPx, fromPx, labelWidthPx: 30 });
    expect(out).toHaveLength(2);
    const center = (peaks[0]!.q + peaks[1]!.q) / 2;
    // Output should be sorted by qLabel for the dodge to make sense; verify
    // one is less than (or equal) to the cluster center and the other is
    // greater than (or equal). Since the dodge is symmetric the average of
    // the dodged positions should remain near the original center.
    const avg = (out[0]!.qLabel + out[1]!.qLabel) / 2;
    expect(Math.abs(avg - center)).toBeLessThan(0.05);
  });
});
