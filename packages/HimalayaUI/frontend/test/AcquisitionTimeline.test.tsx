import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { AcquisitionTimeline, type AcqSession } from "../src/print/components/AcquisitionTimeline";

const SESSIONS: AcqSession[] = [
  { label: "Apr 12", loadFrameCounts: [42, 30, 55] },
  { label: "Apr 13", loadFrameCounts: [48, 52] },
];

describe("AcquisitionTimeline", () => {
  it("renders one bar per load across all sessions", () => {
    render(<AcquisitionTimeline sessions={SESSIONS} />);
    expect(screen.getAllByTestId("acq-bar")).toHaveLength(5);
  });
  it("renders a session label per session", () => {
    render(<AcquisitionTimeline sessions={SESSIONS} />);
    expect(screen.getByText(/APR 12/i)).toBeInTheDocument();
    expect(screen.getByText(/APR 13/i)).toBeInTheDocument();
  });
  it("renders nothing meaningful (empty frame) for no sessions", () => {
    render(<AcquisitionTimeline sessions={[]} />);
    expect(screen.queryAllByTestId("acq-bar")).toHaveLength(0);
  });
  it("SVG carries an explicit height attribute to prevent container-fill stretch", () => {
    render(<AcquisitionTimeline sessions={SESSIONS} />);
    // The SVG must have a fixed height attribute so it does not scale to fill
    // the grid cell height (walkthrough bug: 132px viewBox stretched to ~1300px).
    const svg = screen.getByRole("img", { name: /acquisition timeline/i });
    expect(svg).toHaveAttribute("height");
    const h = Number(svg.getAttribute("height"));
    expect(h).toBeGreaterThan(0);
    expect(h).toBeLessThanOrEqual(200); // sanity: must be compact, not viewport-filling
  });
});
