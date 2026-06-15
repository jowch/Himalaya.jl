import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { NoticePill } from "../../src/print/ui/NoticePill";

describe("NoticePill", () => {
  it("renders its children text", () => {
    render(<NoticePill tone="new">+2 new match</NoticePill>);
    expect(screen.getByTestId("notice-pill").textContent).toBe("+2 new match");
  });

  it('has data-testid="notice-pill"', () => {
    render(<NoticePill tone="new">+1 new match</NoticePill>);
    expect(screen.getByTestId("notice-pill")).toBeTruthy();
  });

  it('tone="new" sets data-tone="new"', () => {
    render(<NoticePill tone="new">+2 new match</NoticePill>);
    expect(screen.getByTestId("notice-pill").getAttribute("data-tone")).toBe("new");
  });

  it('tone="draft" sets data-tone="draft"', () => {
    render(<NoticePill tone="draft">Draft</NoticePill>);
    expect(screen.getByTestId("notice-pill").getAttribute("data-tone")).toBe("draft");
  });

  it("forwards a placement className", () => {
    render(<NoticePill tone="new" className="ml-4">+1</NoticePill>);
    expect(screen.getByTestId("notice-pill").className).toContain("ml-4");
  });
});
