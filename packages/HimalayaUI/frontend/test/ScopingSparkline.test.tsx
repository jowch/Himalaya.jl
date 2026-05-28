import { describe, it, expect } from "vitest";
import { render } from "@testing-library/react";
import { ScopingSparkline } from "../src/components/ScopingSparkline";

const trace = (() => {
  const q = Array.from({ length: 40 }, (_, i) => 0.02 + i * 0.01);
  return { q, I: q.map((x) => Math.exp(-((x - 0.1) ** 2) / 0.0002)), sigma: q.map(() => 0) };
})();

describe("ScopingSparkline", () => {
  it("renders an svg path when a trace is present", () => {
    const { container } = render(<ScopingSparkline trace={trace} phase="Pn3m" />);
    expect(container.querySelector("path")).not.toBeNull();
  });
  it("renders an empty frame (no path) when the trace is undefined", () => {
    const { container, getByTestId } = render(<ScopingSparkline trace={undefined} phase={null} />);
    expect(getByTestId("scoping-sparkline-empty")).toBeInTheDocument();
    expect(container.querySelector("path")).toBeNull();
  });
});
