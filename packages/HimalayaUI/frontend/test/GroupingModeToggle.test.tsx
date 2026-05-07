/**
 * GroupingModeToggle tests (Plan §Phase 9, Task 9.2).
 *
 * Spec selectors:
 *   container `data-testid="grouping-mode"` with `data-mode={"bySample"|"byPhase"|"distinct"}`
 *   per-option `data-value` + `data-active` + `role="radio"` + `aria-checked`.
 *
 * Behaviour:
 *   - Default reads from Zustand `groupingMode` (defaults to "bySample").
 *   - Click dispatches `setGroupingMode`.
 *   - Three options render: "By sample" / "By phase" / "Distinct".
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { GroupingModeToggle } from "../src/components/GroupingModeToggle";
import { useAppState } from "../src/state";

describe("GroupingModeToggle", () => {
  beforeEach(() => {
    // Reset to default before every test.
    useAppState.setState({ groupingMode: "bySample" });
  });

  it("renders the testid + active-mode reflection", () => {
    render(<GroupingModeToggle />);
    const root = screen.getByTestId("grouping-mode");
    expect(root).toBeInTheDocument();
    expect(root).toHaveAttribute("data-mode", "bySample");
  });

  it("renders three options labelled correctly", () => {
    render(<GroupingModeToggle />);
    expect(screen.getByText("By sample")).toBeInTheDocument();
    expect(screen.getByText("By phase")).toBeInTheDocument();
    expect(screen.getByText("Distinct")).toBeInTheDocument();
  });

  it("clicking 'By phase' updates Zustand and reflects on the container", () => {
    render(<GroupingModeToggle />);
    fireEvent.click(screen.getByText("By phase"));
    expect(useAppState.getState().groupingMode).toBe("byPhase");
    expect(screen.getByTestId("grouping-mode")).toHaveAttribute("data-mode", "byPhase");
  });

  it("clicking 'Distinct' updates Zustand", () => {
    render(<GroupingModeToggle />);
    fireEvent.click(screen.getByText("Distinct"));
    expect(useAppState.getState().groupingMode).toBe("distinct");
  });

  it("active option carries data-active='true' and aria-checked='true'", () => {
    useAppState.setState({ groupingMode: "byPhase" });
    render(<GroupingModeToggle />);
    const byPhaseBtn = screen.getByText("By phase").closest("button")!;
    expect(byPhaseBtn).toHaveAttribute("data-active", "true");
    expect(byPhaseBtn).toHaveAttribute("aria-checked", "true");

    const distinctBtn = screen.getByText("Distinct").closest("button")!;
    expect(distinctBtn).toHaveAttribute("data-active", "false");
    expect(distinctBtn).toHaveAttribute("aria-checked", "false");
  });
});
