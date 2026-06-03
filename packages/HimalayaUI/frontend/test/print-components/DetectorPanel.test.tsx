import { render, screen } from "@testing-library/react";
import { DetectorPanel } from "../../src/print/components/DetectorPanel";

describe("DetectorPanel", () => {
  it("renders the default label and the detector frame", () => {
    render(<DetectorPanel src={null} />);
    expect(screen.getByTestId("detector-panel")).toBeInTheDocument();
    expect(screen.getByText("Detector image")).toBeInTheDocument();
    expect(screen.getByTestId("detector-frame")).toBeInTheDocument();
  });
  it("renders a custom label", () => {
    render(<DetectorPanel src={null} label="Real source" />);
    expect(screen.getByText("Real source")).toBeInTheDocument();
  });
  it("renders the header tools slot (exposure strip)", () => {
    render(<DetectorPanel src={null} tools={<div data-testid="expo">e</div>} />);
    expect(screen.getByTestId("expo")).toBeInTheDocument();
  });
  it("renders a hint line when provided", () => {
    render(<DetectorPanel src={null} hint="The real source." />);
    expect(screen.getByText("The real source.")).toBeInTheDocument();
  });
  it("forwards a placement-only className", () => {
    render(<DetectorPanel src={null} className="h-full" />);
    expect(screen.getByTestId("detector-panel").className).toContain("h-full");
  });
});
