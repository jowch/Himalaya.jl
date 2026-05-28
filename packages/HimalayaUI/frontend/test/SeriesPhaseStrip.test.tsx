import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import type { SeriesMember } from "../src/api";
import { SeriesPhaseStrip } from "../src/components/SeriesPhaseStrip";

function member(phase: string | null, order: number): SeriesMember {
  return {
    id: order + 1, series_id: 1, exposure_id: order + 100, display_order: order,
    band_height: 1, y_offset: 0, normalization: "max", color_override: null,
    label_override: null, q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: phase === null ? null : {
      effective_peaks: [{ id: 1, q: 0.04, intensity: 1, sharpness: 1, source: "auto" as const }],
      confirmed_index: { id: order + 1, phase, lattice_d: 100, r_squared: 0.99, ngc: 3, peak_ids: [1] },
      analysis_inputs_hash: "h",
    },
    is_stale: false, created_by: 1, created_at: null,
  };
}

describe("SeriesPhaseStrip", () => {
  it("renders one segment cell per member", () => {
    render(<SeriesPhaseStrip members={[member("Pn3m", 0), member("Lamellar", 1), member("Im3m", 2)]} />);
    expect(screen.getAllByTestId("ps-seg")).toHaveLength(3);
  });

  it("captions a transition as 'first → last'", () => {
    render(<SeriesPhaseStrip members={[member("Pn3m", 0), member("Lamellar", 1)]} />);
    const cap = screen.getByTestId("ps-cap");
    expect(cap).toHaveTextContent("Pn3m");
    expect(cap).toHaveTextContent("Lamellar");
    expect(cap).toHaveTextContent("→");
  });

  it("captions a single phase as 'phase throughout'", () => {
    render(<SeriesPhaseStrip members={[member("Pn3m", 0), member("Pn3m", 1)]} />);
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/throughout/i);
  });

  it("captions an all-unindexed strip as no clear phase", () => {
    render(<SeriesPhaseStrip members={[member(null, 0), member(null, 1)]} />);
    expect(screen.getByTestId("ps-cap")).toHaveTextContent(/no clear phase/i);
  });
});
