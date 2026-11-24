import { describe, it, expect, vi } from "vitest";
import { fireEvent, render, screen } from "@testing-library/react";
import { SortMenu } from "../../src/ui/SortMenu";
import { SAMPLE_SORTS } from "../../src/workspace/sampleSort";

function open(onSort = vi.fn()) {
  render(<SortMenu what="samples" current="measured" choices={SAMPLE_SORTS} onSort={onSort} />);
  const button = screen.getByRole("button", { name: "Sort samples by Measured" });
  fireEvent.click(button);
  return { button, onSort };
}

describe("SortMenu", () => {
  it("opens on the current order, ticked and focused", () => {
    open();
    const current = screen.getByRole("menuitemradio", { name: /^Measured/ });
    expect(current).toHaveAttribute("aria-checked", "true");
    expect(current).toHaveFocus();
    expect(screen.getByRole("menuitemradio", { name: /^Name/ })).toHaveAttribute("aria-checked", "false");
  });

  it("moves with the arrow keys and picks with Enter, then closes and gives focus back", () => {
    const { button, onSort } = open();
    fireEvent.keyDown(screen.getByRole("menu"), { key: "ArrowDown" });
    const name = screen.getByRole("menuitemradio", { name: /^Name/ });
    expect(name).toHaveFocus();
    fireEvent.click(name);
    expect(onSort).toHaveBeenCalledWith("name");
    expect(screen.queryByRole("menu")).toBeNull();
    expect(button).toHaveFocus();
  });

  it("closes on Esc without changing the order, and keeps workspace keys out while open", () => {
    const { button, onSort } = open();
    expect(screen.getByRole("menu")).toHaveAttribute("data-keys-trap");
    fireEvent.keyDown(screen.getByRole("menu"), { key: "Escape" });
    expect(screen.queryByRole("menu")).toBeNull();
    expect(onSort).not.toHaveBeenCalled();
    expect(button).toHaveFocus();
  });
});
