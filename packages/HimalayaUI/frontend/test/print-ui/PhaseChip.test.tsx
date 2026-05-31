import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { PhaseChip } from "../../src/print/ui/PhaseChip";

describe("<PhaseChip> coexistence", () => {
  it("renders single phase by default", () => {
    render(<PhaseChip phase="Pn3m" />);
    expect(screen.getByTestId("phase-chip").textContent).toBe("Pn3m");
  });
  it("renders two-phase label with dominant first", () => {
    render(<PhaseChip phase="Pn3m" coexistWith="Lamellar" />);
    const chip = screen.getByTestId("phase-chip");
    expect(chip.textContent).toBe("Pn3m + Lam");
    expect(chip.getAttribute("data-coexist")).toBe("true");
  });
});
