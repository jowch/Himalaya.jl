// test/print-components/PanelHeader.test.tsx
import { render, screen } from "@testing-library/react";
import { PanelHeader } from "../../src/print/components/PanelHeader";

describe("PanelHeader", () => {
  it("renders the label text", () => {
    render(<PanelHeader label="Detector image" />);
    expect(screen.getByText("Detector image")).toBeInTheDocument();
  });
  it("renders a right-side tools slot when children are given", () => {
    render(
      <PanelHeader label="Reflections — comb">
        <button data-testid="tool">x</button>
      </PanelHeader>,
    );
    expect(screen.getByTestId("panel-header")).toBeInTheDocument();
    expect(screen.getByTestId("tool")).toBeInTheDocument();
  });
  it("omits the tools slot when no children", () => {
    render(<PanelHeader label="Detector image" />);
    expect(screen.queryByTestId("panel-header-tools")).not.toBeInTheDocument();
  });
  it("forwards a placement-only className", () => {
    render(<PanelHeader label="X" className="mb-5" />);
    expect(screen.getByTestId("panel-header").className).toContain("mb-5");
  });
});
