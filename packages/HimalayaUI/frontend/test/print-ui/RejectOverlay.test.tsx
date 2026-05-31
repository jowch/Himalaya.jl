import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { RejectOverlay } from "../../src/print/ui/RejectOverlay";

describe("RejectOverlay", () => {
  it('renders an aria-hidden svg', () => {
    render(<RejectOverlay />);
    const el = screen.getByTestId("reject-overlay");
    expect(el).not.toBeNull();
    expect(el.getAttribute("aria-hidden")).toBe("true");
  });

  it('has data-testid="reject-overlay"', () => {
    render(<RejectOverlay />);
    expect(screen.getByTestId("reject-overlay")).not.toBeNull();
  });

  it("contains exactly two line elements", () => {
    render(<RejectOverlay />);
    expect(screen.getByTestId("reject-overlay").querySelectorAll("line").length).toBe(2);
  });
});
