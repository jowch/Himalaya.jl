import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ProgressBar } from "../../src/print/ui/ProgressBar";

describe("<ProgressBar>", () => {
  it("exposes role=progressbar with aria value attributes", () => {
    render(<ProgressBar value={8} total={12} label="Screened 8 of 12" />);
    const bar = screen.getByRole("progressbar");
    expect(bar.getAttribute("aria-valuenow")).toBe("8");
    expect(bar.getAttribute("aria-valuemin")).toBe("0");
    expect(bar.getAttribute("aria-valuemax")).toBe("12");
    expect(bar.getAttribute("aria-label")).toBe("Screened 8 of 12");
  });

  it("sets the fill width to the value/total ratio", () => {
    render(<ProgressBar value={8} total={12} />);
    const fill = screen.getByTestId("progressbar").querySelector("span");
    expect(fill?.style.width).toMatch(/^66\.6/);
  });

  it("clamps the fill width at 100% when value exceeds total", () => {
    render(<ProgressBar value={20} total={12} />);
    const fill = screen.getByTestId("progressbar").querySelector("span");
    expect(fill?.style.width).toBe("100%");
  });

  it("renders 0% (no divide-by-zero) when total is 0", () => {
    render(<ProgressBar value={5} total={0} />);
    const fill = screen.getByTestId("progressbar").querySelector("span");
    expect(fill?.style.width).toBe("0%");
  });

  it("omits aria-label when no label is provided", () => {
    render(<ProgressBar value={1} total={2} />);
    expect(screen.getByRole("progressbar").getAttribute("aria-label")).toBe(null);
  });
});
