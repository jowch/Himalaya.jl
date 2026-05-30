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
    // The wide transparent hit-ring carries the hover handlers (the thin
    // visible ring is pointer-events:none so the pointer reliably lands).
    const hit = screen.getByTestId("detector-ring-hit-0.206");
    fireEvent.mouseEnter(hit);
    expect(useAppState.getState().hoveredQ).toBe(0.206);
    fireEvent.mouseLeave(hit);
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

describe("DetectorRingOverlay — phase-colour rings (Plan D-6)", () => {
  beforeEach(() => { useAppState.setState({ hoveredQ: undefined }); });

  it("colours rings by their assigned phase", () => {
    render(<DetectorRingOverlay rings={[
      { q: 0.045, color: "rgb(10, 20, 30)" },
      { q: 0.055, color: "rgb(40, 50, 60)" },
    ]} />);
    const r1 = screen.getByTestId("detector-ring-q-0.045");
    expect(r1.getAttribute("stroke")).toBe("rgb(10, 20, 30)");
    const r2 = screen.getByTestId("detector-ring-q-0.055");
    expect(r2.getAttribute("stroke")).toBe("rgb(40, 50, 60)");
  });

  it("renders concentric sets under coexistence (two phases, distinct colours)", () => {
    render(<DetectorRingOverlay rings={[
      { q: 0.045, color: "rgb(10, 20, 30)" },  // phase A
      { q: 0.064, color: "rgb(40, 50, 60)" },  // phase B
    ]} />);
    expect(screen.getAllByTestId(/^detector-ring-q-/)).toHaveLength(2);
    expect(screen.getByTestId("detector-ring-q-0.045").getAttribute("stroke"))
      .toBe("rgb(10, 20, 30)");
    expect(screen.getByTestId("detector-ring-q-0.064").getAttribute("stroke"))
      .toBe("rgb(40, 50, 60)");
  });

  it("renders a hollow ghost ring for a predicted-but-absent order", () => {
    render(<DetectorRingOverlay rings={[
      { q: 0.045, color: "rgb(10, 20, 30)" },
      { q: 0.0636, color: "rgb(10, 20, 30)", ghost: true },
    ]} />);
    const ghost = screen.getByTestId("detector-ring-q-0.0636");
    expect(ghost).toHaveAttribute("data-ghost", "true");
    // ghost is dashed (hollow) — distinguishes it from a solid observed ring.
    expect(ghost.getAttribute("stroke-dasharray")).toBeTruthy();
  });

  it("q-link lights the phase-coloured ring matching hoveredQ", () => {
    render(<DetectorRingOverlay rings={[
      { q: 0.045, color: "rgb(10, 20, 30)" },
      { q: 0.055, color: "rgb(40, 50, 60)" },
    ]} />);
    act(() => useAppState.getState().setHoveredQ(0.055));
    expect(screen.getByTestId("detector-ring-q-0.055")).toHaveAttribute("data-hot", "true");
    expect(screen.getByTestId("detector-ring-q-0.045")).toHaveAttribute("data-hot", "false");
  });
});
