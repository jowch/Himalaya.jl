import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ReadingRow } from "../../src/print/components/ReadingRow";

describe("<ReadingRow> rendering", () => {
  it("renders the root [data-testid=reading-row]", () => {
    render(<ReadingRow phase="Pn3m" span="1:0 → 1:0.75" lattice="a 205 → 195 Å" />);
    expect(screen.getByTestId("reading-row")).toBeInTheDocument();
  });

  it("renders a circle Swatch carrying the phase", () => {
    const { container } = render(
      <ReadingRow phase="Pn3m" span="1:0 → 1:0.75" lattice="a 205 → 195 Å" />
    );
    const swatch = container.querySelector("[data-swatch]");
    expect(swatch).toBeInTheDocument();
    expect(swatch).toHaveAttribute("data-phase", "Pn3m");
    expect(swatch).toHaveAttribute("data-shape", "circle");
  });

  it("renders the phase name inside [data-phase-label]", () => {
    const { container } = render(
      <ReadingRow phase="Pn3m" span="1:0 → 1:0.75" lattice="a 205 → 195 Å" />
    );
    const label = container.querySelector("[data-phase-label]");
    expect(label).toBeInTheDocument();
    expect(label).toHaveTextContent("Pn3m");
    expect(label).toHaveAttribute("data-phase", "Pn3m");
  });

  it("renders the span text", () => {
    render(<ReadingRow phase="Pn3m" span="1:0 → 1:0.75" lattice="a 205 → 195 Å" />);
    expect(screen.getByText("1:0 → 1:0.75")).toBeInTheDocument();
  });

  it("renders the lattice text", () => {
    render(<ReadingRow phase="Pn3m" span="1:0 → 1:0.75" lattice="a 205 → 195 Å" />);
    expect(screen.getByText("a 205 → 195 Å")).toBeInTheDocument();
  });
});
