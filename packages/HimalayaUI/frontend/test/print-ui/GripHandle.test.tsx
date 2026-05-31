import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { GripHandle } from "../../src/print/ui/GripHandle";

describe("GripHandle", () => {
  it("renders the ⋮⋮ text", () => {
    render(<GripHandle />);
    expect(screen.getByTestId("grip-handle").textContent).toBe("⋮⋮");
  });

  it('is aria-hidden', () => {
    render(<GripHandle />);
    expect(screen.getByTestId("grip-handle").getAttribute("aria-hidden")).toBe("true");
  });

  it('has data-testid="grip-handle"', () => {
    render(<GripHandle />);
    expect(screen.getByTestId("grip-handle")).not.toBeNull();
  });
});
