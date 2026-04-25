import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ScoreBar } from "../../src/components/ui/ScoreBar";

describe("<ScoreBar>", () => {
  it("renders a bar with width proportional to score", () => {
    render(<ScoreBar score={0.75} color="#d97706" />);
    const bar = document.querySelector("[data-score-bar]") as HTMLElement;
    expect(bar).not.toBeNull();
    expect(bar.style.width).toBe("75%");
  });

  it("clamps score to 100% at maximum", () => {
    render(<ScoreBar score={1.0} color="#22c55e" />);
    const bar = document.querySelector("[data-score-bar]") as HTMLElement;
    expect(bar.style.width).toBe("100%");
  });

  it("renders an empty bar for score of 0", () => {
    render(<ScoreBar score={0} color="#22c55e" />);
    const bar = document.querySelector("[data-score-bar]") as HTMLElement;
    expect(bar.style.width).toBe("0%");
  });

  it("applies the given color as background", () => {
    render(<ScoreBar score={0.5} color="rgb(59, 130, 246)" />);
    const bar = document.querySelector("[data-score-bar]") as HTMLElement;
    expect(bar.style.background).toBe("rgb(59, 130, 246)");
  });
});
