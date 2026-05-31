import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Menu } from "../../src/print/ui/Menu";

const opts = [
  { value: "a", label: "Apple" },
  { value: "b", label: "Banana" },
  { value: "c", label: "Cherry" },
] as const;

describe("<Menu>", () => {
  it("renders N menuitems when open", () => {
    render(
      <Menu aria-label="fruit" open options={opts} onSelect={() => {}} onClose={() => {}} />,
    );
    expect(screen.getAllByRole("menuitem")).toHaveLength(3);
  });

  it("renders nothing when closed", () => {
    render(
      <Menu aria-label="fruit" open={false} options={opts} onSelect={() => {}} onClose={() => {}} />,
    );
    expect(screen.queryByRole("menu")).toBeNull();
  });

  it("clicking an option fires onSelect(value) then onClose", () => {
    const onSelect = vi.fn();
    const onClose = vi.fn();
    render(
      <Menu aria-label="fruit" open options={opts} onSelect={onSelect} onClose={onClose} />,
    );
    fireEvent.click(screen.getByRole("menuitem", { name: "Banana" }));
    expect(onSelect).toHaveBeenCalledWith("b");
    expect(onClose).toHaveBeenCalledTimes(1);
    expect(onSelect.mock.invocationCallOrder[0]).toBeLessThan(onClose.mock.invocationCallOrder[0]);
  });

  it("Escape fires onClose", () => {
    const onClose = vi.fn();
    render(
      <Menu aria-label="fruit" open options={opts} onSelect={() => {}} onClose={onClose} />,
    );
    fireEvent.keyDown(screen.getByRole("menu"), { key: "Escape" });
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("ArrowDown moves focus to the next menuitem", () => {
    render(
      <Menu aria-label="fruit" open options={opts} onSelect={() => {}} onClose={() => {}} />,
    );
    const items = screen.getAllByRole("menuitem");
    items[0].focus();
    expect(document.activeElement).toBe(items[0]);
    fireEvent.keyDown(items[0], { key: "ArrowDown" });
    expect(document.activeElement).toBe(items[1]);
  });

  it("a disabled option does not fire onSelect on click", () => {
    const onSelect = vi.fn();
    const onClose = vi.fn();
    render(
      <Menu
        aria-label="fruit"
        open
        options={[
          { value: "a", label: "Apple" },
          { value: "b", label: "Banana", disabled: true },
        ]}
        onSelect={onSelect}
        onClose={onClose}
      />,
    );
    fireEvent.click(screen.getByRole("menuitem", { name: "Banana" }));
    expect(onSelect).not.toHaveBeenCalled();
  });

  it("the active option has data-active=true", () => {
    render(
      <Menu
        aria-label="fruit"
        open
        options={opts}
        activeValue="b"
        onSelect={() => {}}
        onClose={() => {}}
      />,
    );
    expect(screen.getByRole("menuitem", { name: "Banana" }).getAttribute("data-active")).toBe("true");
    expect(screen.getByRole("menuitem", { name: "Apple" }).getAttribute("data-active")).toBeNull();
  });
});
