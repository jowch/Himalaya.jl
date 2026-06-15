import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Wordmark } from "../../src/print/ui/Wordmark";

describe("<Wordmark>", () => {
  it("renders the brand text", () => {
    render(<Wordmark>Himalaya</Wordmark>);
    expect(screen.getByTestId("wordmark").textContent).toContain("Himalaya");
  });
  it("renders an optional tail in a data-role=tail span", () => {
    render(<Wordmark tail="SAXS">Himalaya</Wordmark>);
    const tail = screen.getByTestId("wordmark").querySelector('[data-role="tail"]');
    expect(tail).not.toBeNull();
    expect(tail?.textContent).toContain("SAXS");
  });
  it("has no tail span when tail omitted", () => {
    render(<Wordmark>Himalaya</Wordmark>);
    expect(screen.getByTestId("wordmark").querySelector('[data-role="tail"]')).toBeNull();
  });
});
