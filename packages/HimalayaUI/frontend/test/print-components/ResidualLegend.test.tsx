import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ResidualLegend } from "../../src/print/components/ResidualLegend";

describe("<ResidualLegend>", () => {
  it("renders the three residual-vocabulary labels", () => {
    render(<ResidualLegend />);
    expect(screen.getByText("in tolerance")).toBeInTheDocument();
    expect(screen.getByText("out of tolerance")).toBeInTheDocument();
    expect(screen.getByText("off scale")).toBeInTheDocument();
  });

  it("renders one glyph per entry, matching the chart's shape vocabulary", () => {
    const { container } = render(<ResidualLegend />);
    const glyphs = container.querySelectorAll('[data-role="resid-legend-glyph"]');
    expect(glyphs.length).toBe(3);
    expect(glyphs[0].querySelector("[data-shape]")!.getAttribute("data-shape")).toBe("dot-filled");
    expect(glyphs[1].querySelector("[data-shape]")!.getAttribute("data-shape")).toBe("dot-hollow");
    expect(glyphs[2].querySelector("[data-shape]")!.getAttribute("data-shape")).toBe("chevron");
  });

  it("names the band and track, derived from the chart constants", () => {
    render(<ResidualLegend />);
    // ±2.2% band / ±3% track come from RESID_BAND / RESID_DOMAIN — if the chart's
    // constants change, this copy follows automatically and this pin updates with it.
    expect(screen.getByText("band ±2.2% · track ±3%")).toBeInTheDocument();
  });
});
