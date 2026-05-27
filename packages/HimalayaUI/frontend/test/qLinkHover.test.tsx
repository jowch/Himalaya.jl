import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { DetectorRingOverlay } from "../src/components/DetectorRingOverlay";
import { QNumInput } from "../src/components/PlotCard";
import { useAppState } from "../src/state";

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ hoveredQ: undefined });
});

// The full PlotCard<->TraceViewer<->detector round-trip is exercised by the
// Playwright mocked spec (e2e/qlink.spec.ts), which has real layout + Plot
// scales. Here we assert the Zustand q-link contract at the seam: a ring
// hover lights the matching ring, and a focused QNumInput is undisturbed by
// hoveredQ churn (the focus-gate non-interference).
describe("q-link hover contract", () => {
  it("ring hover -> hoveredQ -> matching ring lights", () => {
    render(<DetectorRingOverlay peakQs={[0.045, 0.103]} />);
    fireEvent.mouseEnter(screen.getByTestId("detector-ring-q-0.045"));
    expect(useAppState.getState().hoveredQ).toBe(0.045);
    expect(screen.getByTestId("detector-ring-q-0.045"))
      .toHaveAttribute("data-hot", "true");
  });

  it("hoveredQ churn does NOT disturb a focused QNumInput (focus-gate)", () => {
    const onCommit = vi.fn();
    // QNumInput sets data-testid={testId} directly (PlotCard.tsx), so the
    // query key IS the testId — mirrors test/PlotCard.test.tsx.
    render(<QNumInput value={0.1} onCommit={onCommit} testId="qmin" />);
    const input = screen.getByTestId("qmin") as HTMLInputElement;
    input.focus();
    fireEvent.change(input, { target: { value: "0.1234" } });
    // q-link activity fires while the user edits:
    act(() => useAppState.getState().setHoveredQ(0.3));
    act(() => useAppState.getState().setHoveredQ(undefined));
    // draft preserved — hoveredQ never feeds QNumInput's value:
    expect(input.value).toBe("0.1234");
    expect(onCommit).not.toHaveBeenCalled();
  });
});
