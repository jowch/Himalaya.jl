import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { AssignmentRail } from "../../src/print/components/AssignmentRail";

describe("AssignmentRail", () => {
  it("renders an aside with Assignment (count) and Candidates sections wrapping the slots", () => {
    const { getByTestId, getByText, getAllByTestId } = render(
      <AssignmentRail
        assignmentCount={2}
        assignment={<div data-testid="slot-cart">CART</div>}
        candidates={<div data-testid="slot-cands">CANDS</div>}
        candidatesNote="ranked by fit"
      />,
    );
    expect(getByTestId("assignment-rail").tagName).toBe("ASIDE");
    const heads = getAllByTestId("rail-section-head").map((h) => h.textContent);
    expect(heads.some((t) => t?.includes("Assignment") && t?.includes("2"))).toBe(true);
    expect(heads.some((t) => t?.includes("Candidates"))).toBe(true);
    expect(getByTestId("slot-cart")).toBeTruthy();
    expect(getByTestId("slot-cands")).toBeTruthy();
    expect(getByText("ranked by fit")).toBeTruthy();
  });
});
