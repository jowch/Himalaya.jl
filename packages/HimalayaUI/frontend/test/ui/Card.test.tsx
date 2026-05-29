import { render, screen } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Card } from "../../src/components/ui/Card";

describe("Card", () => {
  it("renders children inside a div by default", () => {
    render(<Card data-testid="c">hello</Card>);
    const el = screen.getByTestId("c");
    expect(el.tagName).toBe("DIV");
    expect(el).toHaveTextContent("hello");
  });

  it("is flat by default — no data-elevated attribute", () => {
    render(<Card data-testid="c">x</Card>);
    expect(screen.getByTestId("c")).not.toHaveAttribute("data-elevated");
  });

  it("elevated sets data-elevated=\"true\"", () => {
    render(<Card data-testid="c" elevated>x</Card>);
    expect(screen.getByTestId("c")).toHaveAttribute("data-elevated", "true");
  });

  it("passes placement className through to the root", () => {
    render(<Card data-testid="c" className="mb-5 max-w-[760px]">x</Card>);
    expect(screen.getByTestId("c").className).toContain("mb-5");
    expect(screen.getByTestId("c").className).toContain("max-w-[760px]");
  });

  it("renders as a button when as=\"button\" and forwards button props", () => {
    const onClick = vi.fn();
    render(<Card as="button" data-testid="c" onClick={onClick}>go</Card>);
    const el = screen.getByTestId("c");
    expect(el.tagName).toBe("BUTTON");
    el.click();
    expect(onClick).toHaveBeenCalledOnce();
  });

  it("renders as li / section when requested", () => {
    render(<Card as="li" data-testid="li">x</Card>);
    render(<Card as="section" data-testid="sec">y</Card>);
    expect(screen.getByTestId("li").tagName).toBe("LI");
    expect(screen.getByTestId("sec").tagName).toBe("SECTION");
  });
});
