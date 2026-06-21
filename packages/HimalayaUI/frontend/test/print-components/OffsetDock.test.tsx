import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { OffsetDock } from "../../src/print/components/OffsetDock";
describe("<OffsetDock>", () => {
  it("renders the offset value and an offset slider", () => {
    render(<OffsetDock offset={1.2} onOffsetChange={() => {}} />);
    expect(screen.getByTestId("offset-dock")).toHaveTextContent("1.20×");
    expect(screen.getByRole("slider", { name: /trace offset/i })).toBeInTheDocument();
  });
  it("fires onOffsetChange", () => {
    const onOffsetChange = vi.fn();
    render(<OffsetDock offset={1.2} onOffsetChange={onOffsetChange} />);
    fireEvent.change(screen.getByRole("slider", { name: /trace offset/i }), { target: { value: "0.8" } });
    expect(onOffsetChange).toHaveBeenCalled();
  });
});
