import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Kicker } from "../../src/print/ui/Kicker";

describe("Kicker", () => {
  it("renders its children text", () => {
    render(<Kicker>Integration</Kicker>);
    expect(screen.getByText("Integration")).toBeInTheDocument();
  });

  it("exposes the tone via data-tone (accent)", () => {
    render(<Kicker tone="accent">Folio</Kicker>);
    expect(screen.getByText("Folio")).toHaveAttribute("data-tone", "accent");
  });

  it("defaults tone to faint", () => {
    render(<Kicker>Notes</Kicker>);
    expect(screen.getByText("Notes")).toHaveAttribute("data-tone", "faint");
  });

  it("renders a heading when as='h3'", () => {
    render(<Kicker as="h3" tone="accent">Folio</Kicker>);
    expect(screen.getByRole("heading", { name: "Folio" })).toBeInTheDocument();
  });

  it("renders no heading by default (inline label)", () => {
    render(<Kicker>Sort</Kicker>);
    expect(screen.queryByRole("heading")).not.toBeInTheDocument();
  });

  it("renders a span when as='span'", () => {
    render(<Kicker as="span">Offset</Kicker>);
    expect(screen.getByText("Offset").tagName).toBe("SPAN");
  });

  it("forwards data-testid", () => {
    render(<Kicker data-testid="focus-plot-kicker">Integration</Kicker>);
    expect(screen.getByTestId("focus-plot-kicker")).toBeInTheDocument();
  });

  it("forwards aria-hidden", () => {
    render(<Kicker aria-hidden="true" data-testid="heatmap-axis-title">x →</Kicker>);
    expect(screen.getByTestId("heatmap-axis-title")).toHaveAttribute("aria-hidden", "true");
  });

  it("accepts placement-only className (margin) and applies it", () => {
    render(<Kicker className="mb-2">samples screened</Kicker>);
    expect(screen.getByText("samples screened")).toHaveClass("mb-2");
  });
});
