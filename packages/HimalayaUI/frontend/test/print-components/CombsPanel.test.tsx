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
  it("derives the label from the view: indexing space names itself", () => {
    render(<CombsPanel {...base} view="resid" />);
    expect(screen.getByText("Reflections · indexing space")).toBeInTheDocument();
    expect(screen.queryByText("Reflections · comb")).not.toBeInTheDocument();
  });
  it("an explicit label prop overrides the view-derived label", () => {
    render(<CombsPanel {...base} view="resid" label="Custom label" />);
    expect(screen.getByText("Custom label")).toBeInTheDocument();
    expect(screen.queryByText("Reflections · indexing space")).not.toBeInTheDocument();
  });
  it("renders the comb body and the comb legend in comb view", () => {
    render(<CombsPanel {...base} />);
    expect(screen.getByTestId("combs-body")).toBeInTheDocument();
    // Comb vocabulary present; residual vocabulary absent.
    expect(screen.getByText("predicted, absent")).toBeInTheDocument();
    expect(screen.queryByText("in tolerance")).not.toBeInTheDocument();
    expect(screen.queryByText(/band ±/)).not.toBeInTheDocument();
  });
  it("renders the residual legend (not the comb legend) in resid view", () => {
    render(<CombsPanel {...base} view="resid" />);
    // The three real residual entries plus the band/track note.
    expect(screen.getByText("in tolerance")).toBeInTheDocument();
    expect(screen.getByText("out of tolerance")).toBeInTheDocument();
    expect(screen.getByText("off scale")).toBeInTheDocument();
    expect(screen.getByText("band ±2.2% · track ±3%")).toBeInTheDocument();
    // Comb vocabulary absent.
    expect(screen.queryByText("predicted, absent")).not.toBeInTheDocument();
    expect(screen.queryByText("observed")).not.toBeInTheDocument();
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
