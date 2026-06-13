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
      // activeValue makes this a value-selector → menuitemradio.
      expect(document.activeElement).toBe(screen.getByRole("menuitemradio", { name: "Banana" }));
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
      // Value-selector (activeValue given): the enabled fallback is a radio too.
      expect(document.activeElement).toBe(screen.getByRole("menuitemradio", { name: "Apple" }));
    });
  });

  it("menu items are a single tab stop: each carries tabIndex=-1 (UI-MENUTAB)", () => {
    // APG menu pattern — the popup is reached via the trigger + arrow keys, not
    // by tabbing through items. Keeping items out of the tab order means one Tab
    // from a focused item leaves the whole menu, letting the owner's focus-out
    // handler close it instead of walking item-to-item.
    render(
      <Menu aria-label="fruit" open options={opts} onSelect={() => {}} onClose={() => {}} />,
    );
    for (const item of screen.getAllByRole("menuitem")) {
      expect(item).toHaveAttribute("tabindex", "-1");
    }
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
    // Value-selector menu: items are radios, so query by menuitemradio role.
    expect(screen.getByRole("menuitemradio", { name: "Banana" }).getAttribute("data-active")).toBe("true");
    expect(screen.getByRole("menuitemradio", { name: "Apple" }).getAttribute("data-active")).toBeNull();
  });

  // A value-selector menu (activeValue given) must tell AT which option is the
  // current value — role=menuitemradio + aria-checked, not a visual-only
  // highlight. An action menu (no activeValue) keeps the plain menuitem.
  describe("value-selector vs action discriminator (activeValue)", () => {
    it("with activeValue: options are menuitemradio with aria-checked marking the active one", () => {
      render(
        <Menu
          aria-label="order by"
          open
          options={opts}
          activeValue="b"
          onSelect={() => {}}
          onClose={() => {}}
        />,
      );
      const radios = screen.getAllByRole("menuitemradio");
      expect(radios).toHaveLength(3);
      // No plain menuitems remain.
      expect(screen.queryAllByRole("menuitem")).toHaveLength(0);
      expect(screen.getByRole("menuitemradio", { name: "Banana" }).getAttribute("aria-checked")).toBe("true");
      expect(screen.getByRole("menuitemradio", { name: "Apple" }).getAttribute("aria-checked")).toBe("false");
      expect(screen.getByRole("menuitemradio", { name: "Cherry" }).getAttribute("aria-checked")).toBe("false");
    });

    it("without activeValue: options stay plain menuitems with no aria-checked (action menu)", () => {
      render(
        <Menu aria-label="download formats" open options={opts} onSelect={() => {}} onClose={() => {}} />,
      );
      const items = screen.getAllByRole("menuitem");
      expect(items).toHaveLength(3);
      expect(screen.queryAllByRole("menuitemradio")).toHaveLength(0);
      for (const it of items) expect(it.getAttribute("aria-checked")).toBeNull();
    });
  });
});
