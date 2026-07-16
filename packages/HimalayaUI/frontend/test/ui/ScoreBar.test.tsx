import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ScoreBar } from "../../src/print/ui/ScoreBar";

describe("<ScoreBar>", () => {
  const fill = () => document.querySelector("[data-score-bar]") as HTMLElement;

  it("renders width proportional to value", () => {
    render(<ScoreBar value={0.75} phase="Pn3m" />);
    expect(fill().style.width).toBe("75%");
  });

  it("clamps value above 1 to 100%", () => {
    render(<ScoreBar value={1.4} phase="Im3m" />);
    expect(fill().style.width).toBe("100%");
  });

  it("clamps a negative value to 0%", () => {
    render(<ScoreBar value={-0.2} phase="Ia3d" />);
    expect(fill().style.width).toBe("0%");
  });

  it("renders an empty bar for value 0", () => {
    render(<ScoreBar value={0} phase="Lamellar" />);
    expect(fill().style.width).toBe("0%");
  });

  it("exposes the phase via data-phase so color derives internally", () => {
    render(<ScoreBar value={0.5} phase="Hexagonal" />);
    expect(fill().getAttribute("data-phase")).toBe("Hexagonal");
  });
});
