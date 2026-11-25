import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import { SubtractToggle } from "../../src/ui/SubtractToggle";
import { minusExplanation, RAW } from "../../src/ui/explanations";

describe("SubtractToggle", () => {
  it("explains Raw and Minus on their buttons", () => {
    render(<SubtractToggle name="buffer_blank" on={false} onChange={() => {}} />);
    expect(screen.getByRole("button", { name: "Raw" })).toHaveAccessibleDescription(RAW);
    expect(screen.getByRole("button", { name: "Minus buffer_blank" })).toHaveAccessibleDescription(minusExplanation("buffer_blank"));
  });

  it("adds no Tab stops of its own", () => {
    const { container } = render(<SubtractToggle name="buffer_blank" on={false} onChange={() => {}} />);
    expect(container.querySelectorAll("[tabindex]")).toHaveLength(0);
  });
});
