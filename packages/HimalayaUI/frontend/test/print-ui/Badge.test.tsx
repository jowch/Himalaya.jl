import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Badge } from "../../src/print/ui/Badge";

describe("Badge", () => {
  it("renders its children text", () => {
    render(<Badge>3</Badge>);
    expect(screen.getByTestId("badge").textContent).toBe("3");
  });

  it('has data-testid="badge"', () => {
    render(<Badge>1</Badge>);
    expect(screen.getByTestId("badge")).toBeTruthy();
  });

  it("forwards a placement className", () => {
    render(<Badge className="ml-4">1</Badge>);
    expect(screen.getByTestId("badge").className).toContain("ml-4");
  });

  it("forwards arbitrary rest props (e.g. data-*)", () => {
    render(<Badge data-foo="bar">1</Badge>);
    expect(screen.getByTestId("badge").getAttribute("data-foo")).toBe("bar");
  });
});
