import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { BonnetBadge } from "../../src/print/ui/BonnetBadge";

describe("BonnetBadge", () => {
  it('renders the text "Bonnet"', () => {
    render(<BonnetBadge />);
    expect(screen.getByText("Bonnet")).toBeTruthy();
    expect(screen.getByTestId("bonnet-badge").textContent).toContain("Bonnet");
  });

  it('has data-testid="bonnet-badge"', () => {
    render(<BonnetBadge />);
    expect(screen.getByTestId("bonnet-badge")).toBeTruthy();
  });

  it("renders an aria-hidden glyph with textContent ⭙", () => {
    render(<BonnetBadge />);
    const glyph = screen
      .getByTestId("bonnet-badge")
      .querySelector('[aria-hidden="true"]');
    expect(glyph).toBeTruthy();
    expect(glyph?.textContent).toBe("⭙");
  });
});
