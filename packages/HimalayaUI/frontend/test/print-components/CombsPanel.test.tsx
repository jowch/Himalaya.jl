import { render, screen, fireEvent } from "@testing-library/react";
import { CombsPanel } from "../../src/print/components/CombsPanel";
import { PN3M, IM3M, LEFTOVER } from "../../src/print/comb/comb.fixtures";

const base = {
  assigned: [PN3M, IM3M],
  leftover: LEFTOVER,
  view: "comb" as const,
  onViewChange: () => {},
};

describe("CombsPanel", () => {
  it("renders the label and the comb-view toggle", () => {
    render(<CombsPanel {...base} />);
    expect(screen.getByTestId("combs-panel")).toBeInTheDocument();
    expect(screen.getByText("Reflections · comb")).toBeInTheDocument();
    expect(screen.getByText("comb")).toBeInTheDocument();
    expect(screen.getByText("indexing space")).toBeInTheDocument();
  });
  it("renders the comb body and the legend in comb view", () => {
    render(<CombsPanel {...base} />);
    expect(screen.getByTestId("combs-body")).toBeInTheDocument();
  });
  it("calls onViewChange when the toggle is used", () => {
    const onViewChange = vi.fn();
    render(<CombsPanel {...base} onViewChange={onViewChange} />);
    fireEvent.click(screen.getByText("indexing space"));
    expect(onViewChange).toHaveBeenCalledWith("resid");
  });
  it("reflects the resid view via data-view", () => {
    render(<CombsPanel {...base} view="resid" />);
    expect(screen.getByTestId("combs-panel").dataset.view).toBe("resid");
  });
  it("forwards a placement-only className", () => {
    render(<CombsPanel {...base} className="h-full" />);
    expect(screen.getByTestId("combs-panel").className).toContain("h-full");
  });
});
