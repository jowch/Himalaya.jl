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

  // APG menu-button: opening the menu moves DOM focus into it.
  describe("focus on open", () => {
    it("focuses the active item when activeValue matches an enabled option", () => {
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
      expect(document.activeElement).toBe(screen.getByRole("menuitem", { name: "Banana" }));
    });

    it("focuses the first enabled item when there is no activeValue", () => {
      render(
        <Menu aria-label="fruit" open options={opts} onSelect={() => {}} onClose={() => {}} />,
      );
      expect(document.activeElement).toBe(screen.getByRole("menuitem", { name: "Apple" }));
    });

    it("skips a disabled first option", () => {
      render(
        <Menu
          aria-label="fruit"
          open
          options={[
            { value: "a", label: "Apple", disabled: true },
            { value: "b", label: "Banana" },
          ]}
          onSelect={() => {}}
          onClose={() => {}}
        />,
      );
      expect(document.activeElement).toBe(screen.getByRole("menuitem", { name: "Banana" }));
    });

    it('initialFocus="last" focuses the last enabled item (ArrowUp-opens)', () => {
      render(
        <Menu
          aria-label="fruit"
          open
          options={opts}
          initialFocus="last"
          onSelect={() => {}}
          onClose={() => {}}
        />,
      );
      expect(document.activeElement).toBe(screen.getByRole("menuitem", { name: "Cherry" }));
    });

    it('initialFocus="last" skips a disabled last option', () => {
      render(
        <Menu
          aria-label="fruit"
          open
          options={[
            { value: "a", label: "Apple" },
            { value: "b", label: "Banana" },
            { value: "c", label: "Cherry", disabled: true },
          ]}
          initialFocus="last"
          onSelect={() => {}}
          onClose={() => {}}
        />,
      );
      expect(document.activeElement).toBe(screen.getByRole("menuitem", { name: "Banana" }));
    });

    it("does not match a disabled option as the active item", () => {
      render(
        <Menu
          aria-label="fruit"
          open
          options={[
            { value: "a", label: "Apple" },
            { value: "b", label: "Banana", disabled: true },
          ]}
          activeValue="b"
          onSelect={() => {}}
          onClose={() => {}}
        />,
      );
      expect(document.activeElement).toBe(screen.getByRole("menuitem", { name: "Apple" }));
    });
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
