import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { KbKey } from "../../src/print/ui/KbKey";

describe("KbKey", () => {
  it("renders the children text", () => {
    render(<KbKey>⌘K</KbKey>);
    expect(screen.getByTestId("kbkey").textContent).toBe("⌘K");
  });

  it("is a <kbd> element", () => {
    render(<KbKey>X</KbKey>);
    expect(screen.getByTestId("kbkey").tagName.toLowerCase()).toBe("kbd");
  });

  it('has data-testid="kbkey"', () => {
    render(<KbKey>X</KbKey>);
    expect(screen.getByTestId("kbkey")).toBeTruthy();
  });
});
