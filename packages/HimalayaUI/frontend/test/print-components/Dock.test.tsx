import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Dock } from "../../src/print/components/Dock";
describe("<Dock>", () => {
  it("renders the offset value and an offset slider", () => {
    render(<Dock offset={1.2} onOffsetChange={() => {}} />);
    expect(screen.getByTestId("dock")).toHaveTextContent("1.20×");
    expect(screen.getByRole("slider", { name: /trace offset/i })).toBeInTheDocument();
  });
  it("fires onOffsetChange", () => {
    const onOffsetChange = vi.fn();
    render(<Dock offset={1.2} onOffsetChange={onOffsetChange} />);
    fireEvent.change(screen.getByRole("slider", { name: /trace offset/i }), { target: { value: "0.8" } });
    expect(onOffsetChange).toHaveBeenCalled();
  });
});
