import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { RailSection } from "../../src/print/components/RailSection";

describe("RailSection", () => {
  it("renders the label, optional count, children, and note", () => {
    const { getByTestId, getByText } = render(
      <RailSection label="Assignment" count={3} note="ranked note">
        <div>body</div>
      </RailSection>,
    );
    const head = getByTestId("rail-section-head");
    expect(head.textContent).toContain("Assignment");
    expect(head.textContent).toContain("3");
    expect(getByText("body")).toBeTruthy();
    expect(getByTestId("rail-section-note").textContent).toContain("ranked note");
  });
  it("omits count and note when not given", () => {
    const { queryByTestId, getByTestId } = render(
      <RailSection label="Candidates"><div>x</div></RailSection>,
    );
    expect(getByTestId("rail-section-head").textContent).toContain("Candidates");
    expect(queryByTestId("rail-section-note")).toBeNull();
  });
});
