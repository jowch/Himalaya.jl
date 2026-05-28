import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { OffsetDock } from "../src/components/OffsetDock";

describe("OffsetDock", () => {
  it("does not render when hidden", () => {
    render(<OffsetDock show={false} value={1.2} onChange={vi.fn()} />);
    expect(screen.queryByTestId("offset-dock")).not.toBeInTheDocument();
  });

  it("renders the dock + value when shown", () => {
    render(<OffsetDock show value={1.2} onChange={vi.fn()} />);
    expect(screen.getByTestId("offset-dock")).toBeInTheDocument();
    expect(screen.getByTestId("offset-dock-value")).toHaveTextContent("1.20");
  });

  it("calls onChange with the parsed value on input", () => {
    const onChange = vi.fn();
    render(<OffsetDock show value={1.2} onChange={onChange} />);
    fireEvent.change(screen.getByTestId("offset-dock-slider"), { target: { value: "0.6" } });
    expect(onChange).toHaveBeenCalledWith(0.6);
  });
});
