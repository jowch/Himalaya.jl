import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { AutogroupCard } from "../src/components/AutogroupCard";

describe("AutogroupCard", () => {
  it("reads N samples and the ordering variable into the body", () => {
    render(
      <AutogroupCard sampleCount={6} orderingVariable="LL37 : lipid ratio" onAdjust={vi.fn()} />,
    );
    expect(screen.getByTestId("autogroup-card")).toHaveTextContent("6 samples");
    expect(screen.getByTestId("autogroup-card")).toHaveTextContent("LL37 : lipid ratio");
  });

  it("falls back to a placeholder when ordering variable is null", () => {
    render(<AutogroupCard sampleCount={3} orderingVariable={null} onAdjust={vi.fn()} />);
    expect(screen.getByTestId("autogroup-card")).toHaveTextContent("their names");
  });

  // M2 "stop the controls lying": the "Confirm series" button was a no-op on a
  // read-only surface (the series is already persisted; SeriesRecipeEditor owns
  // real save/commit). Removed. "Adjust grouping" is the one real follow-up.
  it("offers Adjust grouping and no longer renders a no-op Confirm", () => {
    const onAdjust = vi.fn();
    render(<AutogroupCard sampleCount={6} orderingVariable="dose" onAdjust={onAdjust} />);
    expect(screen.queryByTestId("autogroup-confirm")).toBeNull();
    fireEvent.click(screen.getByTestId("autogroup-adjust"));
    expect(onAdjust).toHaveBeenCalled();
  });
});
