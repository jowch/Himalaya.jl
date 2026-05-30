import { render, screen, cleanup } from "@testing-library/react";
import { afterEach, describe, it, expect } from "vitest";
import { PhaseChip } from "../../src/components/ui/PhaseChip";

afterEach(cleanup);

describe("<PhaseChip>", () => {
  it("always renders the phase name as text (the colour-blind second channel)", () => {
    render(<PhaseChip phase="Pn3m" />);
    expect(screen.getByText("Pn3m")).toBeInTheDocument();
  });

  it("renders the name for every variant × size combination", () => {
    const combos = [
      ["tint", "sm"],
      ["tint", "md"],
      ["solid", "sm"],
      ["solid", "md"],
    ] as const;
    for (const [variant, size] of combos) {
      const { unmount } = render(
        <PhaseChip phase="Im3m" variant={variant} size={size} />,
      );
      expect(screen.getByText("Im3m")).toBeInTheDocument();
      unmount();
    }
  });

  it("defaults to data-testid=phase-chip and exposes data-variant/data-size", () => {
    render(<PhaseChip phase="Ia3d" />);
    const chip = screen.getByTestId("phase-chip");
    expect(chip).toHaveTextContent("Ia3d");
    expect(chip).toHaveAttribute("data-variant", "tint");
    expect(chip).toHaveAttribute("data-size", "sm");
  });

  it("reflects an explicit variant and size on the data-attributes", () => {
    render(<PhaseChip phase="Square" variant="solid" size="md" />);
    const chip = screen.getByTestId("phase-chip");
    expect(chip).toHaveAttribute("data-variant", "solid");
    expect(chip).toHaveAttribute("data-size", "md");
  });

  it("lets a consumer override the testid via spread", () => {
    render(<PhaseChip phase="Hexagonal" data-testid="member-meta-phase-chip" />);
    expect(screen.getByTestId("member-meta-phase-chip")).toHaveTextContent(
      "Hexagonal",
    );
  });

  it("renders an unknown phase (phaseColor FALLBACK) without throwing", () => {
    render(<PhaseChip phase="Bogus" />);
    expect(screen.getByText("Bogus")).toBeInTheDocument();
  });

  it("appends a placement className to the element", () => {
    render(<PhaseChip phase="Lamellar" className="ml-2" />);
    expect(screen.getByTestId("phase-chip").className).toContain("ml-2");
  });
});
