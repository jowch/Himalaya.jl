import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { PageFrame, PAGE_WIDTHS } from "../../src/print/components/PageFrame";

describe("PageFrame", () => {
  it("renders children inside a centered frame", () => {
    render(<PageFrame width="loupe"><p>body</p></PageFrame>);
    const frame = screen.getByTestId("page-frame");
    expect(frame).toHaveTextContent("body");
    expect(frame.className).toContain("mx-auto");
  });
  it("applies the named max-width class for each surface", () => {
    for (const [name, cls] of Object.entries(PAGE_WIDTHS)) {
      const { unmount } = render(<PageFrame width={name as keyof typeof PAGE_WIDTHS}>x</PageFrame>);
      expect(screen.getByTestId("page-frame").className).toContain(cls);
      unmount();
    }
  });
  it("forwards a placement className", () => {
    render(<PageFrame width="sheet" className="py-7">x</PageFrame>);
    expect(screen.getByTestId("page-frame").className).toContain("py-7");
  });
});
