import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ScreenedMark } from "../../src/print/ui/ScreenedMark";

describe("ScreenedMark", () => {
  it('has role="img" with aria-label="Screened" when done', () => {
    render(<ScreenedMark done />);
    const el = screen.getByTestId("screened-mark");
    expect(el.getAttribute("role")).toBe("img");
    expect(el.getAttribute("aria-label")).toBe("Screened");
  });

  it('has role="img" with aria-label="Not screened" when not done', () => {
    render(<ScreenedMark done={false} />);
    const el = screen.getByTestId("screened-mark");
    expect(el.getAttribute("role")).toBe("img");
    expect(el.getAttribute("aria-label")).toBe("Not screened");
  });

  it('sets data-done="true" only when done', () => {
    render(<ScreenedMark done />);
    expect(screen.getByTestId("screened-mark").getAttribute("data-done")).toBe("true");
  });

  it("omits data-done when not done", () => {
    render(<ScreenedMark done={false} />);
    expect(screen.getByTestId("screened-mark").getAttribute("data-done")).toBeNull();
  });

  it("renders the check svg (and its path) only when done", () => {
    render(<ScreenedMark done />);
    const svg = screen.getByTestId("screened-mark").querySelector("svg");
    expect(svg).not.toBeNull();
    expect(svg?.querySelector("path")).not.toBeNull();
  });

  it("renders no svg when not done", () => {
    render(<ScreenedMark done={false} />);
    expect(screen.getByTestId("screened-mark").querySelector("svg")).toBeNull();
  });

  it("forwards a placement className", () => {
    render(<ScreenedMark done className="ml-4" />);
    expect(screen.getByTestId("screened-mark").className).toContain("ml-4");
  });
});
