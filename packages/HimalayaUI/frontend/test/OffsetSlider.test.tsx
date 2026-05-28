import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { OffsetSlider } from "../src/components/OffsetSlider";

describe("OffsetSlider", () => {
  it("renders the current value to two decimals", () => {
    render(<OffsetSlider value={1.2} onChange={vi.fn()} />);
    expect(screen.getByTestId("offset-slider-value")).toHaveTextContent("1.20");
  });

  it("calls onChange with the parsed numeric value on input", () => {
    const onChange = vi.fn();
    render(<OffsetSlider value={1.2} onChange={onChange} />);
    fireEvent.change(screen.getByTestId("offset-slider"), { target: { value: "0.8" } });
    expect(onChange).toHaveBeenCalledWith(0.8);
  });

  it("exposes the slider range and step matching the mockup", () => {
    render(<OffsetSlider value={1.2} onChange={vi.fn()} />);
    const input = screen.getByTestId("offset-slider") as HTMLInputElement;
    expect(input.min).toBe("0.4");
    expect(input.max).toBe("1.4");
    expect(input.step).toBe("0.05");
  });
});
