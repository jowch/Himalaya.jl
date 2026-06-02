import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Dot } from "../../src/print/ui/Dot";

describe("<Dot> (print)", () => {
  it("is an img with the label when labelled", () => {
    const { container } = render(<Dot label="kept" tone="success" />);
    const dot = container.querySelector("span")!;
    expect(dot.getAttribute("role")).toBe("img");
    expect(dot.getAttribute("aria-label")).toBe("kept");
    expect(dot.getAttribute("data-tone")).toBe("success");
  });

  it("drops role=img and aria-label when decorative (aria-hidden)", () => {
    const { container } = render(<Dot tone="accent" size="xs" aria-hidden />);
    const dot = container.querySelector("span")!;
    expect(dot.getAttribute("role")).toBeNull();
    expect(dot.getAttribute("aria-label")).toBeNull();
    expect(dot.getAttribute("aria-hidden")).toBe("true");
  });

  it("renders the md size class (9px circle)", () => {
    const { container } = render(<Dot tone="accent" size="md" aria-hidden />);
    const dot = container.querySelector("span")!;
    expect(dot.className).toContain("h-[9px]");
    expect(dot.className).toContain("w-[9px]");
  });

  it("adds a plate-colored border ring when bordered", () => {
    const { container } = render(<Dot tone="accent" size="md" bordered aria-hidden />);
    const dot = container.querySelector("span")!;
    expect(dot.className).toContain("border-[1.5px]");
    expect(dot.className).toContain("border-plate");
  });

  it("has no border ring by default", () => {
    const { container } = render(<Dot tone="accent" size="md" aria-hidden />);
    const dot = container.querySelector("span")!;
    expect(dot.className).not.toContain("border-plate");
  });
});
