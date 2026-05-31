import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Card } from "../../src/print/ui/Card";

describe("<Card> draft variant", () => {
  it("is not draft by default", () => {
    render(<Card data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-draft")).toBe(null);
  });
  it("marks data-draft when draft", () => {
    render(<Card draft data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-draft")).toBe("true");
  });
  it("still honors elevated alongside the existing API", () => {
    render(<Card elevated data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-elevated")).toBe("true");
  });
});

describe("<Card> padding", () => {
  it("has no data-padding by default", () => {
    render(<Card data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-padding")).toBe(null);
  });
  it("sets data-padding when padding given", () => {
    render(<Card padding="md" data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-padding")).toBe("md");
  });
});
