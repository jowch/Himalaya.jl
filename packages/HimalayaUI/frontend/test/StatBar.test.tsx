// test/StatBar.test.tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { StatBar } from "../src/print/ui/StatBar";

describe("StatBar (Phase E1 primitive)", () => {
  it("renders one cell per stat with value + caption", () => {
    render(
      <StatBar
        aria-label="Experiment stats"
        stats={[
          { key: "loads", caption: "Loads", value: "13" },
          { key: "samples", caption: "Samples", value: "170" },
          { key: "exposures", caption: "Exposures", value: "682" },
          { key: "sessions", caption: "Sessions", value: "4" },
        ]}
      />,
    );
    const cells = screen.getAllByTestId("statbar-cell");
    expect(cells).toHaveLength(4);
    expect(cells[0]).toHaveTextContent("Loads");
    expect(cells[0]).toHaveTextContent("13");
  });

  it("marks a pending cell for assistive tech + a data flag", () => {
    render(
      <StatBar
        aria-label="Experiment stats"
        stats={[{ key: "span", caption: "Span", value: "pending", pending: true }]}
      />,
    );
    expect(screen.getByTestId("statbar-cell")).toHaveAttribute("data-pending", "true");
  });

  it("exposes the aria-label on the band", () => {
    render(<StatBar aria-label="X" stats={[{ key: "a", caption: "A", value: "1" }]} />);
    expect(screen.getByTestId("statbar")).toHaveAttribute("aria-label", "X");
  });
});
