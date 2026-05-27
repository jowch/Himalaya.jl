import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { DetectorRingOverlay } from "../src/components/DetectorRingOverlay";
import { useAppState } from "../src/state";

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ hoveredQ: undefined });
});

const PEAK_QS = [0.045, 0.103, 0.206];

describe("DetectorRingOverlay", () => {
  it("renders one ring per peak q", () => {
    render(<DetectorRingOverlay peakQs={PEAK_QS} />);
    expect(screen.getAllByTestId(/^detector-ring-q-/)).toHaveLength(3);
  });

  it("lights the ring matching hoveredQ", () => {
    render(<DetectorRingOverlay peakQs={PEAK_QS} />);
    act(() => useAppState.getState().setHoveredQ(0.103));
    const hot = screen.getByTestId("detector-ring-q-0.103");
    expect(hot).toHaveAttribute("data-hot", "true");
    expect(screen.getByTestId("detector-ring-q-0.045"))
      .toHaveAttribute("data-hot", "false");
  });

  it("sets hoveredQ when a ring is hovered, clears on leave", () => {
    render(<DetectorRingOverlay peakQs={PEAK_QS} />);
    const ring = screen.getByTestId("detector-ring-q-0.206");
    fireEvent.mouseEnter(ring);
    expect(useAppState.getState().hoveredQ).toBe(0.206);
    fireEvent.mouseLeave(ring);
    expect(useAppState.getState().hoveredQ).toBeUndefined();
  });

  it("applies the rotation transform in landscape, none in portrait", () => {
    const { rerender } = render(
      <DetectorRingOverlay peakQs={PEAK_QS} orient="portrait" />,
    );
    const svgP = screen.getByTestId("detector-ring-overlay");
    expect(svgP.getAttribute("data-orient")).toBe("portrait");
    rerender(<DetectorRingOverlay peakQs={PEAK_QS} orient="landscape" />);
    const svgL = screen.getByTestId("detector-ring-overlay");
    expect(svgL.getAttribute("data-orient")).toBe("landscape");
    expect(svgL.getAttribute("style") ?? "").toContain("rotate(90deg)");
  });
});
