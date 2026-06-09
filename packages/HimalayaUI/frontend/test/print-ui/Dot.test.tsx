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

  it("reflects the md size on data-size", () => {
    const { container } = render(<Dot tone="accent" size="md" aria-hidden />);
    expect(container.querySelector("span")!).toHaveAttribute("data-size", "md");
  });

  it("marks data-bordered='true' when bordered", () => {
    const { container } = render(<Dot tone="accent" size="md" bordered aria-hidden />);
    expect(container.querySelector("span")!).toHaveAttribute("data-bordered", "true");
  });

  it("has no data-bordered attribute by default", () => {
    const { container } = render(<Dot tone="accent" size="md" aria-hidden />);
    expect(container.querySelector("span")!).not.toHaveAttribute("data-bordered");
  });
  it("passes through arbitrary props (data-testid, title)", () => {
    const { container } = render(
      <Dot tone="accent" size="md" aria-hidden data-testid="d" title="hover me" />,
    );
    const dot = container.querySelector("span")!;
    expect(dot.getAttribute("data-testid")).toBe("d");
    expect(dot.getAttribute("title")).toBe("hover me");
  });
});
