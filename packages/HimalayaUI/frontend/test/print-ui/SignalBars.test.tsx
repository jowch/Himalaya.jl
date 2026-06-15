import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { SignalBars } from "../../src/print/ui/SignalBars";

describe("SignalBars", () => {
  it("renders exactly 5 bars", () => {
    render(<SignalBars value={4} max={5} />);
    expect(screen.getByTestId("signal-bars").querySelectorAll("i").length).toBe(5);
  });

  it("reflects the active count from value/max (4 of 5)", () => {
    render(<SignalBars value={4} max={5} />);
    expect(
      screen.getByTestId("signal-bars").querySelectorAll('[data-on="true"]').length,
    ).toBe(4);
  });

  it("reflects the active count from value/max (1 of 5)", () => {
    render(<SignalBars value={1} max={5} />);
    expect(
      screen.getByTestId("signal-bars").querySelectorAll('[data-on="true"]').length,
    ).toBe(1);
  });

  it('has role="img" with aria-label "Signal 4 of 5"', () => {
    render(<SignalBars value={4} max={5} />);
    const el = screen.getByRole("img", { name: "Signal 4 of 5" });
    expect(el).toBeTruthy();
    expect(el.getAttribute("data-testid")).toBe("signal-bars");
  });

  it("clamps above max (value 10, max 5 → 5 active)", () => {
    render(<SignalBars value={10} max={5} />);
    expect(
      screen.getByTestId("signal-bars").querySelectorAll('[data-on="true"]').length,
    ).toBe(5);
  });

  it("clamps below zero (value -3 → 0 active)", () => {
    render(<SignalBars value={-3} max={5} />);
    expect(
      screen.getByTestId("signal-bars").querySelectorAll('[data-on="true"]').length,
    ).toBe(0);
  });

  it("defaults max to 5 when omitted (value 3 → 3 active)", () => {
    render(<SignalBars value={3} />);
    expect(
      screen.getByTestId("signal-bars").querySelectorAll('[data-on="true"]').length,
    ).toBe(3);
  });
});
