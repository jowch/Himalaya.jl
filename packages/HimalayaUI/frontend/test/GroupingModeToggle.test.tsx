/**
 * GroupingModeToggle tests (Plan §Phase 9, Task 9.2; updated Compare UX C-4).
 *
 * Spec selectors:
 *   container `data-testid="grouping-mode"` with `data-mode={"bySample"|"byPhase"|"distinct"}`
 *   per-option `data-value` + `data-active` + `role="radio"` + `aria-checked`.
 *
 * Behaviour:
 *   - Reflects the `mode` prop on the container's `data-mode`.
 *   - Click dispatches `onChange` with the selected value.
 *   - Three options render: "By sample" / "By phase" / "Distinct".
 *
 * C-4: The component is now prop-driven (mode + onChange) rather than
 * reading/writing Zustand directly. The parent resolves the effective mode
 * via effectiveGroupingMode(draft, comparison).
 */
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { GroupingModeToggle } from "../src/components/GroupingModeToggle";

describe("GroupingModeToggle", () => {
  it("renders the testid + active-mode reflection", () => {
    render(<GroupingModeToggle mode="bySample" onChange={() => {}} />);
    const root = screen.getByTestId("grouping-mode");
    expect(root).toBeInTheDocument();
    expect(root).toHaveAttribute("data-mode", "bySample");
  });

  it("renders three options labelled correctly", () => {
    render(<GroupingModeToggle mode="bySample" onChange={() => {}} />);
    expect(screen.getByText("By sample")).toBeInTheDocument();
    expect(screen.getByText("By phase")).toBeInTheDocument();
    expect(screen.getByText("Distinct")).toBeInTheDocument();
  });

  it("clicking 'By phase' calls onChange with 'byPhase'", () => {
    const onChange = vi.fn();
    render(<GroupingModeToggle mode="bySample" onChange={onChange} />);
    fireEvent.click(screen.getByText("By phase"));
    expect(onChange).toHaveBeenCalledWith("byPhase");
  });

  it("clicking 'Distinct' calls onChange with 'distinct'", () => {
    const onChange = vi.fn();
    render(<GroupingModeToggle mode="bySample" onChange={onChange} />);
    fireEvent.click(screen.getByText("Distinct"));
    expect(onChange).toHaveBeenCalledWith("distinct");
  });

  it("active option carries data-active='true' and aria-checked='true'", () => {
    render(<GroupingModeToggle mode="byPhase" onChange={() => {}} />);
    const byPhaseBtn = screen.getByText("By phase").closest("button")!;
    expect(byPhaseBtn).toHaveAttribute("data-active", "true");
    expect(byPhaseBtn).toHaveAttribute("aria-checked", "true");

    const distinctBtn = screen.getByText("Distinct").closest("button")!;
    expect(distinctBtn).toHaveAttribute("data-active", "false");
    expect(distinctBtn).toHaveAttribute("aria-checked", "false");
  });

  it("reflects mode='distinct' on container data-mode", () => {
    render(<GroupingModeToggle mode="distinct" onChange={() => {}} />);
    expect(screen.getByTestId("grouping-mode")).toHaveAttribute("data-mode", "distinct");
  });
});
