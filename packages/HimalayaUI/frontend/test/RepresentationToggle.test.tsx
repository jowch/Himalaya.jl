import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { RepresentationToggle, type Representation } from "../src/components/RepresentationToggle";

function setup(value: Representation = "waterfall") {
  const onChange = vi.fn();
  render(<RepresentationToggle value={value} onChange={onChange} />);
  return { onChange };
}

describe("RepresentationToggle", () => {
  it("marks the active representation with aria-pressed", () => {
    setup("waterfall");
    expect(screen.getByTestId("repr-waterfall")).toHaveAttribute("aria-pressed", "true");
  });

  it("calls onChange with 'waterfall' when the waterfall button is clicked", () => {
    const { onChange } = setup("waterfall");
    fireEvent.click(screen.getByTestId("repr-waterfall"));
    expect(onChange).toHaveBeenCalledWith("waterfall");
  });

  it("renders heatmap as a disabled, not-yet-available option (deferred to #208)", () => {
    const { onChange } = setup("waterfall");
    const heatmap = screen.getByTestId("repr-heatmap");
    expect(heatmap).toBeDisabled();
    fireEvent.click(heatmap);
    expect(onChange).not.toHaveBeenCalled();
  });
});
