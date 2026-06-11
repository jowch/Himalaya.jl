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

describe("<Card> selected", () => {
  it("has no data-selected by default", () => {
    render(<Card data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-selected")).toBe(null);
  });
  it("sets data-selected='true' when selected=true", () => {
    render(<Card selected data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-selected")).toBe("true");
  });
  it("does not set data-selected when selected=false", () => {
    render(<Card selected={false} data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-selected")).toBe(null);
  });
  it("composes with as='button'", () => {
    render(<Card as="button" selected data-testid="c">x</Card>);
    const el = screen.getByTestId("c");
    expect(el.tagName.toLowerCase()).toBe("button");
    expect(el.getAttribute("data-selected")).toBe("true");
  });
  it("does not set data-selected when selected is absent (existing flat card)", () => {
    render(<Card elevated data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-selected")).toBe(null);
  });
});

describe("<Card> polymorphic `as`", () => {
  it("renders a <div> by default", () => {
    render(<Card data-testid="c">x</Card>);
    expect(screen.getByTestId("c").tagName.toLowerCase()).toBe("div");
  });
  it("renders the given tag for as='li'", () => {
    render(<Card as="li" data-testid="c">x</Card>);
    expect(screen.getByTestId("c").tagName.toLowerCase()).toBe("li");
  });
  it("renders the given tag for as='section'", () => {
    render(<Card as="section" data-testid="c">x</Card>);
    expect(screen.getByTestId("c").tagName.toLowerCase()).toBe("section");
  });
  it("renders the given tag for as='article' (the folio card root)", () => {
    render(<Card as="article" data-testid="c">x</Card>);
    expect(screen.getByTestId("c").tagName.toLowerCase()).toBe("article");
  });
});

describe("<Card> interactive (door affordance)", () => {
  it("is not interactive by default", () => {
    render(<Card data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-interactive")).toBe(null);
  });
  it("sets data-interactive='true' when interactive", () => {
    render(<Card interactive data-testid="c">x</Card>);
    expect(screen.getByTestId("c").getAttribute("data-interactive")).toBe("true");
  });
  it("composes with elevated + as='article' (the lifted folio door)", () => {
    render(<Card as="article" elevated interactive data-testid="c">x</Card>);
    const el = screen.getByTestId("c");
    expect(el.tagName.toLowerCase()).toBe("article");
    expect(el.getAttribute("data-elevated")).toBe("true");
    expect(el.getAttribute("data-interactive")).toBe("true");
  });
});
