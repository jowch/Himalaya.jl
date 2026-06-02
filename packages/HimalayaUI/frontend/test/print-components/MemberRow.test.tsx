import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { MemberRow } from "../../src/print/components/MemberRow";

describe("<MemberRow> rendering", () => {
  it("renders the root [data-testid=member-row]", () => {
    render(<MemberRow name="trace-01" phase="Pn3m" />);
    expect(screen.getByTestId("member-row")).toBeInTheDocument();
  });

  it("renders the grip handle (⋮⋮ glyph, aria-hidden)", () => {
    const { container } = render(<MemberRow name="trace-01" phase="Pn3m" />);
    const grip = container.querySelector("[data-testid='grip-handle']");
    expect(grip).toBeInTheDocument();
    expect(grip).toHaveAttribute("aria-hidden", "true");
    expect(grip).toHaveTextContent("⋮⋮");
  });

  it("renders a square Swatch carrying the phase", () => {
    const { container } = render(<MemberRow name="trace-01" phase="Pn3m" />);
    const swatch = container.querySelector("[data-swatch]");
    expect(swatch).toBeInTheDocument();
    expect(swatch).toHaveAttribute("data-phase", "Pn3m");
    expect(swatch).toHaveAttribute("data-shape", "square");
  });

  it("renders the name text", () => {
    render(<MemberRow name="trace-01" phase="Pn3m" />);
    expect(screen.getByText("trace-01")).toBeInTheDocument();
  });

  it("renders the sub text when passed", () => {
    render(<MemberRow name="trace-01" sub="S-014 · 1:0.5" phase="Pn3m" />);
    expect(screen.getByText("S-014 · 1:0.5")).toBeInTheDocument();
  });

  it("omits the sub line when not passed", () => {
    render(<MemberRow name="trace-01" phase="Pn3m" />);
    expect(screen.queryByText("S-014 · 1:0.5")).not.toBeInTheDocument();
  });

  it("renders the PhaseChip reflecting the phase", () => {
    const { container } = render(<MemberRow name="trace-01" phase="Pn3m" />);
    const chip = container.querySelector("[data-testid='phase-chip']");
    expect(chip).toBeInTheDocument();
    expect(chip).toHaveTextContent("Pn3m");
  });

  it("forwards coexistWith to the PhaseChip", () => {
    const { container } = render(
      <MemberRow name="trace-01" phase="Pn3m" coexistWith={["Im3m"]} />
    );
    const chip = container.querySelector("[data-testid='phase-chip']");
    expect(chip).toHaveAttribute("data-coexist", "true");
    expect(chip).toHaveTextContent("Pn3m + Im3m");
  });
});
