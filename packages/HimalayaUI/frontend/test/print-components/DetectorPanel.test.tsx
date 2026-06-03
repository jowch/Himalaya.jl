import { render, screen, fireEvent } from "@testing-library/react";
import { DetectorPanel } from "../../src/print/components/DetectorPanel";

const RINGS = [
  { q: 0.045, color: "var(--color-pn3m)" },
  { q: 0.064, color: "var(--color-pn3m)" },
];

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
  it("overlays one phase ring per provided reflection", () => {
    const { container } = render(<DetectorPanel src={null} rings={RINGS} />);
    expect(screen.getByTestId("detector-rings")).toBeInTheDocument();
    expect(container.querySelectorAll('[data-role="det-ring"]')).toHaveLength(2);
  });
  it("renders no ring overlay when no rings are provided", () => {
    render(<DetectorPanel src={null} />);
    expect(screen.queryByTestId("detector-rings")).not.toBeInTheDocument();
  });
  it("lights the ring whose q matches hoveredQ (the triple-link)", () => {
    const { container } = render(
      <DetectorPanel src={null} rings={RINGS} hoveredQ={0.045} />,
    );
    const hot = container.querySelector('[data-role="det-ring"][data-hot="true"]');
    expect(hot).not.toBeNull();
    expect(hot?.getAttribute("data-ring-q")).toBe("0.045");
  });
  it("places rings on a measured off-center beamCenter override", () => {
    const { container } = render(
      <DetectorPanel src={null} rings={RINGS} beamCenter={{ x: 0.43, y: 0.2 }} />,
    );
    const ring = container.querySelector('[data-role="ring-sharp"]');
    expect(ring?.getAttribute("cx")).toBe("0.43");
    expect(ring?.getAttribute("cy")).toBe("0.2");
  });
  it("fires onHoverQ with the ring's q when a ring is hovered", () => {
    const onHoverQ = vi.fn();
    const { container } = render(
      <DetectorPanel src={null} rings={RINGS} onHoverQ={onHoverQ} />,
    );
    const hit = container.querySelector('[data-role="ring-hit"]');
    fireEvent.mouseEnter(hit!);
    expect(onHoverQ).toHaveBeenCalledWith(0.045);
  });
});
