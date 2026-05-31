import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { CheckCircle } from "../../src/print/ui/CheckCircle";

describe("CheckCircle", () => {
  it('has role="img" with aria-label="Selected" when checked', () => {
    render(<CheckCircle checked />);
    const el = screen.getByTestId("check-circle");
    expect(el.getAttribute("role")).toBe("img");
    expect(el.getAttribute("aria-label")).toBe("Selected");
  });

  it('has role="img" with aria-label="Not selected" when not checked', () => {
    render(<CheckCircle checked={false} />);
    const el = screen.getByTestId("check-circle");
    expect(el.getAttribute("role")).toBe("img");
    expect(el.getAttribute("aria-label")).toBe("Not selected");
  });

  it("lets a custom label prop override the aria-label", () => {
    render(<CheckCircle checked label="Screened" />);
    expect(screen.getByTestId("check-circle").getAttribute("aria-label")).toBe("Screened");
  });

  it('sets data-checked="true" only when checked', () => {
    render(<CheckCircle checked />);
    expect(screen.getByTestId("check-circle").getAttribute("data-checked")).toBe("true");
  });

  it("omits data-checked when not checked", () => {
    render(<CheckCircle checked={false} />);
    expect(screen.getByTestId("check-circle").getAttribute("data-checked")).toBeNull();
  });

  it("renders the check svg (and its path) only when checked", () => {
    render(<CheckCircle checked />);
    const svg = screen.getByTestId("check-circle").querySelector("svg");
    expect(svg).not.toBeNull();
    expect(svg?.querySelector("path")).not.toBeNull();
  });

  it("renders no svg when not checked", () => {
    render(<CheckCircle checked={false} />);
    expect(screen.getByTestId("check-circle").querySelector("svg")).toBeNull();
  });

  it("forwards a placement className", () => {
    render(<CheckCircle checked className="ml-4" />);
    expect(screen.getByTestId("check-circle").className).toContain("ml-4");
  });
});
