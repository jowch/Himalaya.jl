import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ScaleToggle } from "../src/components/ScaleToggle";

describe("ScaleToggle", () => {
  it("marks the active scale via aria-pressed", () => {
    render(<ScaleToggle value="log" onChange={vi.fn()} />);
    expect(screen.getByTestId("scale-log")).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByTestId("scale-linear")).toHaveAttribute("aria-pressed", "false");
  });

  it("calls onChange('linear') when linear is clicked", () => {
    const onChange = vi.fn();
    render(<ScaleToggle value="log" onChange={onChange} />);
    fireEvent.click(screen.getByTestId("scale-linear"));
    expect(onChange).toHaveBeenCalledWith("linear");
  });

  it("calls onChange('log') when log is clicked", () => {
    const onChange = vi.fn();
    render(<ScaleToggle value="linear" onChange={onChange} />);
    fireEvent.click(screen.getByTestId("scale-log"));
    expect(onChange).toHaveBeenCalledWith("log");
  });
});
