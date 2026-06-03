import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { AssignmentCart } from "../../src/print/components/AssignmentCart";
import { PhaseBlock } from "../../src/print/components/PhaseBlock";

const block = (phase: string) => (
  <PhaseBlock key={phase} phase={phase} score={0.8} meta="a = 14.2 nm · 5 reflections" />
);

describe("AssignmentCart", () => {
  it("shows the empty state when there are no phase blocks", () => {
    render(<AssignmentCart />);
    expect(screen.getByTestId("assignment-empty")).toBeInTheDocument();
    expect(screen.getByText(/no phase assigned/i)).toBeInTheDocument();
    expect(screen.queryByTestId("coexistence-tag")).not.toBeInTheDocument();
  });
  it("lets the caller override the empty copy", () => {
    render(<AssignmentCart empty="Nothing here yet" />);
    expect(screen.getByText("Nothing here yet")).toBeInTheDocument();
  });
  it("renders one block and NO coexistence tag for a single phase", () => {
    render(<AssignmentCart>{block("Im3m")}</AssignmentCart>);
    expect(screen.getAllByTestId("phase-block")).toHaveLength(1);
    expect(screen.queryByTestId("coexistence-tag")).not.toBeInTheDocument();
    expect(screen.getByTestId("assignment-cart").dataset.phaseCount).toBe("1");
  });
  it("shows a 'Coexistence · N phases' tag when more than one phase is present", () => {
    render(<AssignmentCart>{[block("Ia3d"), block("Im3m")]}</AssignmentCart>);
    expect(screen.getAllByTestId("phase-block")).toHaveLength(2);
    const tag = screen.getByTestId("coexistence-tag");
    expect(tag).toHaveTextContent(/coexistence/i);
    expect(tag).toHaveTextContent(/2 phases/);
    expect(screen.getByTestId("assignment-cart").dataset.phaseCount).toBe("2");
  });
  it("draws a divider on every block after the first (full-bleed .pc-block separator)", () => {
    render(<AssignmentCart>{[block("Ia3d"), block("Im3m"), block("Pn3m")]}</AssignmentCart>);
    const blocks = screen.getAllByTestId("cart-block");
    expect(blocks).toHaveLength(3);
    expect(blocks.map((b) => b.dataset.divider)).toEqual(["false", "true", "true"]);
  });
  it("forwards a placement-only className", () => {
    render(<AssignmentCart className="mt-4">{block("Im3m")}</AssignmentCart>);
    expect(screen.getByTestId("assignment-cart").className).toContain("mt-4");
  });
  it("renders the custom-index footer when onCustomIndex is given, and fires it", () => {
    const onCustomIndex = vi.fn();
    const { getByTestId } = render(
      <AssignmentCart onCustomIndex={onCustomIndex}>
        <div>blockA</div>
      </AssignmentCart>,
    );
    const foot = getByTestId("custom-index-trigger");
    expect(foot.textContent).toContain("custom index");
    fireEvent.click(foot);
    expect(onCustomIndex).toHaveBeenCalledTimes(1);
  });
  it("shows the footer in the empty state too", () => {
    const { getByTestId } = render(<AssignmentCart onCustomIndex={() => {}} />);
    expect(getByTestId("assignment-empty")).toBeTruthy();
    expect(getByTestId("custom-index-trigger")).toBeTruthy();
  });
  it("omits the footer when onCustomIndex is absent", () => {
    const { queryByTestId } = render(<AssignmentCart><div>x</div></AssignmentCart>);
    expect(queryByTestId("custom-index-trigger")).toBeNull();
  });
});
