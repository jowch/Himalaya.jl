import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { HintText } from "../../src/print/ui/HintText";

describe("HintText", () => {
  it("renders its children text", () => {
    render(<HintText>No matches found.</HintText>);
    expect(screen.getByText("No matches found.")).toBeTruthy();
  });

  it("renders as a paragraph element", () => {
    render(<HintText>hint</HintText>);
    const el = screen.getByText("hint");
    expect(el.tagName).toBe("P");
  });

  it("forwards a placement className to the root element", () => {
    render(<HintText className="mt-4">placement test</HintText>);
    const el = screen.getByText("placement test");
    expect(el.className).toContain("mt-4");
  });
});
