import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { IndexingPlot, type IndexSeries } from "../../src/ui/IndexingPlot";
import { CURVATURE, R_SQUARED } from "../../src/ui/explanations";

const series = (phase: string): IndexSeries => ({
  key: phase, phase, color: "#888", lattice: 100, r2: 0.9987,
  fit: { points: [{ ratio: 1.41, predicted: 0.089, observed: 0.089 }], slope: 0.063, r2: 0.9987, landed: 1 },
});

describe("IndexingPlot legend", () => {
  it("explains κ and R² as terms a keyboard reaches", () => {
    render(<IndexingPlot series={[series("Pn3m")]} />);
    expect(screen.getByTestId("workspace-indexing-curvature").closest("[tabindex]")).toHaveAccessibleDescription(CURVATURE);
    expect(screen.getByText("R² 0.9987").closest("[tabindex]")).toHaveAccessibleDescription(R_SQUARED);
  });

  it("shows no κ for a phase without a minimal surface", () => {
    render(<IndexingPlot series={[series("Lamellar")]} />);
    expect(screen.queryByTestId("workspace-indexing-curvature")).toBeNull();
  });
});
