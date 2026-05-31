/**
 * SeriesTrackingOverlay (Plan E, Task E-2/E-3) — the SVG overlay that paints
 * the (phase,order) migration tracks + hollow ghost-ring carets and exposes the
 * hover-anchor handles. Renders into PlotSurface's overlay <svg>.
 */
import { describe, it, expect, vi } from "vitest";
import { render, fireEvent } from "@testing-library/react";
import type { SeriesMember } from "../src/api";
import { buildAnchorMap } from "../src/lib/series/anchors";
import { SeriesTrackingOverlay } from "../src/components/SeriesTrackingOverlay";
import type { PlotOverlayContext } from "../src/components/PlotSurface";

function member(
  id: number, phase: string | null, lattice_d: number, observedQs: number[],
): SeriesMember {
  const peakIds = observedQs.map((_, i) => id * 100 + i);
  return {
    id, series_id: 1, exposure_id: id * 10, display_order: id - 1,
    band_height: 1, y_offset: 0, normalization: "qwindow",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: phase
      ? {
          effective_peaks: observedQs.map((q, i) => ({
            id: peakIds[i]!, q, intensity: 50, sharpness: 1, source: "auto" as const,
          })),
          confirmed_index: {
            id: id * 1000, phase, lattice_d, r_squared: 0.99, ngc: -1.5, peak_ids: peakIds,
          },
          analysis_inputs_hash: "h",
        }
      : { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
    is_stale: false, created_by: null, created_at: null,
  };
}

// Identity-ish ctx: q→px is q*1000, y passthrough.
const ctx: PlotOverlayContext = {
  applyQ: (q) => q * 1000,
  invertQ: (px) => px / 1000,
  applyY: (v) => v,
  width: 400, height: 200,
  margins: { left: 50, right: 14, top: 8, bottom: 32 },
  hitTest: () => null,
};

const a = 200;
const qOf = (rad: number) => (2 * Math.PI * Math.sqrt(rad)) / a;

function renderOverlay(extra: Partial<React.ComponentProps<typeof SeriesTrackingOverlay>> = {}) {
  const m1 = member(1, "Pn3m", a, [qOf(2), qOf(4)]);          // absents order 1
  const m2 = member(2, "Pn3m", a, [qOf(2), qOf(3), qOf(4)]);  // full
  const members = [m1, m2];
  const anchorMap = buildAnchorMap(members);
  const yBands: Array<[number, number]> = [[0, 100], [100, 200]];
  return render(
    <svg>
      <SeriesTrackingOverlay
        ctx={ctx}
        members={members}
        anchorMap={anchorMap}
        yBands={yBands}
        {...extra}
      />
    </svg>,
  );
}

describe("SeriesTrackingOverlay", () => {
  it("draws one polyline per (phase,order) track with ≥2 carriers", () => {
    const { container } = renderOverlay();
    const lines = container.querySelectorAll('[data-role="series-track-line"]');
    // Pn3m orders 0,1,2 each have two carriers → 3 polylines.
    expect(lines.length).toBe(3);
  });

  it("plants a hollow caret ghost ring at the predicted-q of an absent order when tracked", () => {
    const { container } = renderOverlay({ trackedKey: "Pn3m:1" });
    const carets = container.querySelectorAll('[data-shape="caret"]');
    // m1 absents order 1 → exactly one ghost caret on the tracked order-1 line.
    expect(carets.length).toBe(1);
  });

  it("renders interactive anchor hit targets for observed peaks", () => {
    const { container } = renderOverlay();
    const hits = container.querySelectorAll('[data-role="series-anchor-hit"]');
    // Observed beads: m1 has 2 (orders 0,2), m2 has 3 (orders 0,1,2) → 5.
    expect(hits.length).toBe(5);
  });

  it("invokes onTrack with the (phase,order) key when an anchor is hovered", () => {
    const onTrack = vi.fn();
    const { container } = renderOverlay({ onTrack });
    const hit = container.querySelector('[data-role="series-anchor-hit"]')!;
    fireEvent.mouseEnter(hit);
    expect(onTrack).toHaveBeenCalledWith(expect.stringMatching(/^Pn3m:\d$/));
  });

  it("highlights only the tracked key's polyline when trackedKey is set", () => {
    const { container } = renderOverlay({ trackedKey: "Pn3m:0" });
    const tracked = container.querySelectorAll('[data-tracked="true"]');
    expect(tracked.length).toBe(1);
  });
});
