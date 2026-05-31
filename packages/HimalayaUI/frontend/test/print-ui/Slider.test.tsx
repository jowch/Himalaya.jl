import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Slider } from "../../src/print/ui/Slider";

describe("<Slider>", () => {
  it("renders a slider role input with min/max/value reflected", () => {
    render(<Slider value={0.4} min={0} max={1.4} onChange={() => {}} />);
    const input = screen.getByRole("slider") as HTMLInputElement;
    expect(input).toHaveAttribute("min", "0");
    expect(input).toHaveAttribute("max", "1.4");
    expect(input.value).toBe("0.4");
  });

  it("reflects the step attribute when given", () => {
    render(<Slider value={0.4} min={0} max={1.4} step={0.1} onChange={() => {}} />);
    const input = screen.getByRole("slider");
    expect(input).toHaveAttribute("step", "0.1");
  });

  it("does not set a step attribute when step is omitted", () => {
    render(<Slider value={0.4} min={0} max={1.4} onChange={() => {}} />);
    const input = screen.getByRole("slider");
    expect(input).not.toHaveAttribute("step");
  });

  it("fires onChange with the NUMERIC value", () => {
    const onChange = vi.fn();
    render(<Slider value={0.4} min={0} max={1.4} step={0.1} onChange={onChange} />);
    const input = screen.getByRole("slider");
    fireEvent.change(input, { target: { value: "0.8" } });
    expect(onChange).toHaveBeenCalledWith(0.8);
    expect(typeof onChange.mock.calls[0]?.[0]).toBe("number");
  });

  it("renders the label when given", () => {
    render(<Slider value={0.4} min={0} max={1.4} onChange={() => {}} label="Offset" />);
    expect(screen.getByText("Offset")).toBeInTheDocument();
  });

  it("renders the valueDisplay when given", () => {
    render(
      <Slider value={0.4} min={0} max={1.4} onChange={() => {}} valueDisplay="0.40" />,
    );
    expect(screen.getByText("0.40")).toBeInTheDocument();
  });
});
