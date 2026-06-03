import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { RailBack } from "../../src/print/components/RailBack";
describe("<RailBack>", () => {
  it("renders the default label and fires onClick", () => {
    const onClick = vi.fn();
    render(<RailBack onClick={onClick} />);
    const b = screen.getByTestId("rail-back");
    expect(b.tagName).toBe("BUTTON");
    expect(b).toHaveTextContent("Compose");
    fireEvent.click(b);
    expect(onClick).toHaveBeenCalledOnce();
  });
  it("accepts a custom label", () => {
    render(<RailBack label="Edit" />);
    expect(screen.getByTestId("rail-back")).toHaveTextContent("Edit");
  });
});
