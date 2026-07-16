import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { CombLegend } from "../../src/print/components/CombLegend";

describe("<CombLegend>", () => {
  it("renders all four labels by default", () => {
    render(<CombLegend />);
    expect(screen.getByText(/observed/i)).toBeInTheDocument();
    expect(screen.getByText(/manual/i)).toBeInTheDocument();
    expect(screen.getByText(/predicted, absent/i)).toBeInTheDocument();
    expect(screen.getByText(/excluded/i)).toBeInTheDocument();
  });

  it("renders one [data-role='peak-glyph'] per shown item by default (4)", () => {
    const { container } = render(<CombLegend />);
    expect(
      container.querySelectorAll('[data-role="peak-glyph"]').length,
    ).toBe(4);
  });

  it("respects items — subset shows only those items", () => {
    const { container } = render(
      <CombLegend items={["auto", "excluded"]} />,
    );
    expect(
      container.querySelectorAll('[data-role="peak-glyph"]').length,
    ).toBe(2);
    expect(screen.getByText(/observed/i)).toBeInTheDocument();
    expect(screen.getByText(/excluded/i)).toBeInTheDocument();
    expect(screen.queryByText(/manual/i)).not.toBeInTheDocument();
    expect(screen.queryByText(/predicted, absent/i)).not.toBeInTheDocument();
  });

  it("renders items in the order given", () => {
    const { container } = render(
      <CombLegend items={["manual", "auto"]} />,
    );
    const glyphs = container.querySelectorAll('[data-role="peak-glyph"]');
    expect(glyphs.length).toBe(2);
    // First item is manual → diamond
    const firstShape = glyphs[0].querySelector("[data-shape]")!;
    expect(firstShape.getAttribute("data-shape")).toBe("diamond");
    // Second item is auto → triangle-down
    const secondShape = glyphs[1].querySelector("[data-shape]")!;
    expect(secondShape.getAttribute("data-shape")).toBe("triangle-down");
  });
});
