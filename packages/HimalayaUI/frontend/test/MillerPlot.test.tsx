import { describe, it, expect, vi } from "vitest";
import { render } from "@testing-library/react";
import type { IndexEntry } from "../src/api";
import { MillerPlot, toScatterData } from "../src/components/MillerPlot";

vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    el.setAttribute("data-testid", "miller-svg");
    return el;
  }),
  dot:   vi.fn(() => ({ _kind: "dot" })),
  linearRegressionY: vi.fn(() => ({ _kind: "linreg" })),
  text:  vi.fn(() => ({ _kind: "text" })),
}));

describe("toScatterData", () => {
  it("maps each peak to (ratio, q_observed, phase)", () => {
    const ix: IndexEntry = {
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.7071, 0.866, 1.0],
      peaks: [
        { peak_id: 10, ratio_position: 1, residual: 0.001, q_observed: 0.71 },
        { peak_id: 11, ratio_position: 3, residual: 0.002, q_observed: 1.01 },
      ],
    };
    const rows = toScatterData([ix]);
    expect(rows).toHaveLength(2);
    expect(rows[0]).toMatchObject({ ratio: 0.7071 / 0.5, q: 0.71, phase: "Pn3m" });
    expect(rows[1]).toMatchObject({ ratio: 1.0 / 0.5, q: 1.01, phase: "Pn3m" });
  });

  it("skips peaks whose ratio_position is out of bounds", () => {
    const ix: IndexEntry = {
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.7],
      peaks: [{ peak_id: 10, ratio_position: 99, residual: 0, q_observed: 1.0 }],
    };
    expect(toScatterData([ix])).toEqual([]);
  });
});

describe("<MillerPlot>", () => {
  it("renders a container with testid 'miller-plot'", () => {
    const { container } = render(<MillerPlot indices={[]} />);
    expect(container.querySelector('[data-testid="miller-plot"]')).not.toBeNull();
  });

  it("invokes Plot.dot and Plot.linearRegressionY when indices have peaks", async () => {
    const Plot = await import("@observablehq/plot");
    (Plot.dot as unknown as { mockClear: () => void }).mockClear();
    (Plot.linearRegressionY as unknown as { mockClear: () => void }).mockClear();
    const ix: IndexEntry = {
      id: 1, exposure_id: 1, phase: "Pn3m", basis: 0.5, score: 1,
      r_squared: 0.99, lattice_d: 12, status: "candidate",
      predicted_q: [0.7, 1.0],
      peaks: [
        { peak_id: 10, ratio_position: 1, residual: 0, q_observed: 0.71 },
        { peak_id: 11, ratio_position: 2, residual: 0, q_observed: 1.0 },
      ],
    };
    render(<MillerPlot indices={[ix]} />);
    expect(Plot.dot).toHaveBeenCalled();
    expect(Plot.linearRegressionY).toHaveBeenCalled();
  });
});
