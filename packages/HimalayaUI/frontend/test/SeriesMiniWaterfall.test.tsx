import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import type { SeriesMember } from "../src/api";
import { SeriesMiniWaterfall } from "../src/components/SeriesMiniWaterfall";

function member(phase: string | null, peaks: number[], order: number): SeriesMember {
  return {
    id: order + 1, series_id: 1, exposure_id: order + 100, display_order: order,
    band_height: 1, y_offset: 0, normalization: "max", color_override: null,
    label_override: null, q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: {
      effective_peaks: peaks.map((q, i) => ({ id: i + 1, q, intensity: 1, sharpness: 1, source: "auto" as const })),
      confirmed_index: phase === null ? null : {
        id: order + 1, phase, lattice_d: 100, r_squared: 0.99, ngc: 3, peak_ids: peaks.map((_, i) => i + 1),
      },
      analysis_inputs_hash: "h",
    },
    is_stale: false, created_by: 1, created_at: null,
  };
}

describe("SeriesMiniWaterfall", () => {
  it("renders an svg tagged with the member (trace) count", () => {
    render(<SeriesMiniWaterfall members={[member("Pn3m", [0.04], 0), member("Lamellar", [0.1], 1)]} />);
    const fig = screen.getByTestId("series-mini-waterfall");
    expect(fig.tagName.toLowerCase()).toBe("svg");
    expect(fig).toHaveAttribute("data-trace-count", "2");
  });

  it("draws one waterfall line path per member", () => {
    const { container } = render(
      <SeriesMiniWaterfall members={[member("Pn3m", [0.04], 0), member("Lamellar", [0.1], 1), member("Im3m", [0.05], 2)]} />,
    );
    expect(container.querySelectorAll("[data-testid='wf-line']")).toHaveLength(3);
  });

  it("renders without throwing for zero members", () => {
    render(<SeriesMiniWaterfall members={[]} />);
    expect(screen.getByTestId("series-mini-waterfall")).toHaveAttribute("data-trace-count", "0");
  });
});
