import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { AutogroupCard } from "../src/components/AutogroupCard";

describe("AutogroupCard", () => {
  it("reads N samples and the ordering variable into the body", () => {
    render(
      <AutogroupCard sampleCount={6} orderingVariable="LL37 : lipid ratio" onConfirm={vi.fn()} onAdjust={vi.fn()} />,
    );
    expect(screen.getByTestId("autogroup-card")).toHaveTextContent("6 samples");
    expect(screen.getByTestId("autogroup-card")).toHaveTextContent("LL37 : lipid ratio");
  });

  it("falls back to a placeholder when ordering variable is null", () => {
    render(<AutogroupCard sampleCount={3} orderingVariable={null} onConfirm={vi.fn()} onAdjust={vi.fn()} />);
    expect(screen.getByTestId("autogroup-card")).toHaveTextContent("their names");
  });

  it("fires onConfirm and onAdjust", () => {
    const onConfirm = vi.fn();
    const onAdjust = vi.fn();
    render(<AutogroupCard sampleCount={6} orderingVariable="dose" onConfirm={onConfirm} onAdjust={onAdjust} />);
    fireEvent.click(screen.getByTestId("autogroup-confirm"));
    fireEvent.click(screen.getByTestId("autogroup-adjust"));
    expect(onConfirm).toHaveBeenCalled();
    expect(onAdjust).toHaveBeenCalled();
  });
});
