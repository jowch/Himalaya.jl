import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Dot } from "../../src/components/ui/Dot";

describe("<Dot>", () => {
  it("is an img with the label when labelled", () => {
    const { container } = render(<Dot label="kept" tone="success" />);
    const dot = container.querySelector("span")!;
    expect(dot.getAttribute("role")).toBe("img");
    expect(dot.getAttribute("aria-label")).toBe("kept");
    expect(dot.getAttribute("data-tone")).toBe("success");
  });

  it("reflects the tone on data-tone", () => {
    const { container } = render(<Dot label="x" tone="accent" />);
    expect(container.querySelector("span")!.getAttribute("data-tone")).toBe("accent");
  });

  it("drops role=img and aria-label when decorative (aria-hidden)", () => {
    const { container } = render(<Dot tone="accent" size="xs" aria-hidden />);
    const dot = container.querySelector("span")!;
    expect(dot.getAttribute("role")).toBeNull();
    expect(dot.getAttribute("aria-label")).toBeNull();
    expect(dot.getAttribute("aria-hidden")).toBe("true");
  });

  it("passes through data-testid and title", () => {
    const { container } = render(<Dot tone="accent" data-testid="stale-dot" title="Has stale members" aria-hidden />);
    const dot = container.querySelector("span")!;
    expect(dot.getAttribute("data-testid")).toBe("stale-dot");
    expect(dot.getAttribute("title")).toBe("Has stale members");
  });
});
