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
});
