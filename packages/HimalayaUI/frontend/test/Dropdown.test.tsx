import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { Dropdown } from "../src/print/ui/Dropdown";

const OPTIONS = [
  { value: "a", label: "Alpha" },
  { value: "b", label: "Beta" },
] as const;

describe("Dropdown (Phase E1 primitive)", () => {
  it("renders the active option's label on the trigger", () => {
    render(<Dropdown aria-label="Pick" options={OPTIONS} value="b" onChange={() => {}} />);
    expect(screen.getByTestId("dropdown-trigger")).toHaveTextContent("Beta");
  });

  it("opens the menu on trigger click and selects an option", () => {
    const onChange = vi.fn();
    render(<Dropdown aria-label="Pick" options={OPTIONS} value="a" onChange={onChange} />);
    expect(screen.queryByRole("menu")).toBeNull();
    fireEvent.click(screen.getByTestId("dropdown-trigger"));
    expect(screen.getByRole("menu")).toBeInTheDocument();
    fireEvent.click(screen.getByText("Beta"));
    expect(onChange).toHaveBeenCalledWith("b");
    // closes after select
    expect(screen.queryByRole("menu")).toBeNull();
  });

  it("renders a placeholder when value matches no option", () => {
    render(
      <Dropdown aria-label="Pick" options={OPTIONS} value={undefined}
        placeholder="Choose…" onChange={() => {}} />,
    );
    expect(screen.getByTestId("dropdown-trigger")).toHaveTextContent("Choose…");
  });

  it("trigger carries aria-haspopup=menu and aria-expanded", () => {
    render(<Dropdown aria-label="Pick" options={OPTIONS} value="a" onChange={() => {}} />);
    const trigger = screen.getByTestId("dropdown-trigger");
    expect(trigger).toHaveAttribute("aria-haspopup", "menu");
    expect(trigger).toHaveAttribute("aria-expanded", "false");
    fireEvent.click(trigger);
    expect(trigger).toHaveAttribute("aria-expanded", "true");
  });
});
