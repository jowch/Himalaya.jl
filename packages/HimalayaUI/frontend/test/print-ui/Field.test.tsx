import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Field } from "../../src/print/ui/Field";

describe("<Field>", () => {
  it("with neither options nor onClick renders a STATIC read-only value (not a button, no chevron)", () => {
    render(<Field value="LL37 : lipid ratio" />);
    const f = screen.getByTestId("field");
    // controls-don't-lie: inert Field is not interactive.
    expect(f.tagName).not.toBe("BUTTON");
    expect(f).toHaveTextContent("LL37 : lipid ratio");
    // No clickable-looking chevron.
    expect(f).not.toHaveTextContent("▾");
  });
  it("renders a clickable trigger button when given onClick", () => {
    const onClick = vi.fn();
    render(<Field value="x" onClick={onClick} />);
    const f = screen.getByTestId("field");
    expect(f.tagName).toBe("BUTTON");
    fireEvent.click(f);
    expect(onClick).toHaveBeenCalledOnce();
  });
  it("shows the placeholder when value is empty", () => {
    render(<Field value="" placeholder="Choose a variable" />);
    expect(screen.getByTestId("field")).toHaveTextContent("Choose a variable");
  });
  it("composes the accessible name from srLabel + value (WCAG 4.1.2)", () => {
    render(<Field value="lipid ratio" srLabel="Ordered by" onClick={() => {}} />);
    // The visually-hidden label prefix joins the value in the computed name —
    // it does NOT override it (which an aria-label attribute would).
    expect(
      screen.getByRole("button", { name: /ordered by lipid ratio/i }),
    ).toBeInTheDocument();
  });
  it("without srLabel the accessible name is just the bare value (additive prop)", () => {
    render(<Field value="lipid ratio" onClick={() => {}} />);
    // No label prefix, the value is never dropped.
    expect(screen.getByRole("button", { name: "lipid ratio" })).toBeInTheDocument();
  });
  it("opens a menu and fires onSelect in dropdown mode", () => {
    const onSelect = vi.fn();
    render(<Field value="Time" options={[{value:"Time",label:"Time"},{value:"Dose",label:"Dose"}]} onSelect={onSelect} />);
    expect(screen.queryByTestId("menu")).toBeNull();          // closed initially
    fireEvent.click(screen.getByTestId("field"));
    expect(screen.getByTestId("menu")).toBeInTheDocument();    // opens
    // Field always passes activeValue → value-selector menu (menuitemradio).
    fireEvent.click(screen.getByRole("menuitemradio", { name: "Dose" }));
    expect(onSelect).toHaveBeenCalledWith("Dose");
    expect(screen.queryByTestId("menu")).toBeNull();          // closes on select
  });

  // APG menu-button keyboard pattern (dropdown mode).
  describe("dropdown keyboard pattern", () => {
    const opts = [
      { value: "Time", label: "Time" },
      { value: "Dose", label: "Dose" },
      { value: "Mass", label: "Mass" },
    ];
    function setup() {
      render(<Field value="Dose" options={opts} onSelect={() => {}} />);
      return screen.getByTestId("field");
    }

    it("ArrowDown on the closed trigger opens and focuses the active item", () => {
      const trigger = setup();
      trigger.focus();
      fireEvent.keyDown(trigger, { key: "ArrowDown" });
      expect(screen.getByTestId("menu")).toBeInTheDocument();
      expect(document.activeElement).toBe(screen.getByRole("menuitemradio", { name: "Dose" }));
    });

    it("ArrowUp on the closed trigger opens and focuses the LAST item", () => {
      const trigger = setup();
      trigger.focus();
      fireEvent.keyDown(trigger, { key: "ArrowUp" });
      expect(screen.getByTestId("menu")).toBeInTheDocument();
      expect(document.activeElement).toBe(screen.getByRole("menuitemradio", { name: "Mass" }));
    });

    it("opening by click also moves focus into the menu (APG)", () => {
      const trigger = setup();
      fireEvent.click(trigger);
      expect(document.activeElement).toBe(screen.getByRole("menuitemradio", { name: "Dose" }));
    });

    it("Escape on the trigger while open closes; focus stays on the trigger", () => {
      const trigger = setup();
      trigger.focus();
      fireEvent.click(trigger);
      trigger.focus(); // user shift-tabbed back / focus never left in this path
      fireEvent.keyDown(trigger, { key: "Escape" });
      expect(screen.queryByTestId("menu")).toBeNull();
      expect(document.activeElement).toBe(trigger);
    });

    it("Escape inside the menu closes and RETURNS focus to the trigger", () => {
      const trigger = setup();
      fireEvent.click(trigger);
      fireEvent.keyDown(screen.getByTestId("menu"), { key: "Escape" });
      expect(screen.queryByTestId("menu")).toBeNull();
      expect(document.activeElement).toBe(trigger);
    });

    it("focus leaving the widget (Tab away) closes the menu (UI-MENUTAB)", () => {
      const trigger = setup();
      fireEvent.click(trigger);
      expect(screen.getByTestId("menu")).toBeInTheDocument();
      // APG: Tab closes the menu. jsdom does not move focus on Tab; model the
      // focusout to an element OUTSIDE the field widget (real Tab is e2e).
      fireEvent.focusOut(screen.getByRole("menuitemradio", { name: "Dose" }), {
        relatedTarget: document.body,
      });
      expect(screen.queryByTestId("menu")).toBeNull();
    });

    it("focus moving to the trigger (inside the widget) keeps the menu open (no toggle-reopen)", () => {
      const trigger = setup();
      fireEvent.click(trigger);
      fireEvent.focusOut(screen.getByRole("menuitemradio", { name: "Dose" }), {
        relatedTarget: trigger,
      });
      expect(screen.getByTestId("menu")).toBeInTheDocument();
    });

    it("selecting an item closes and returns focus to the trigger", () => {
      const trigger = setup();
      fireEvent.click(trigger);
      fireEvent.click(screen.getByRole("menuitemradio", { name: "Mass" }));
      expect(screen.queryByTestId("menu")).toBeNull();
      expect(document.activeElement).toBe(trigger);
    });

    it("Escape that closes the menu does NOT reach document listeners (innermost popup only)", () => {
      // ModalShell dismisses on a document-level Escape listener; a menu
      // inside a future modal must not double-close it (APG: Escape closes
      // only the innermost popup).
      const trigger = setup();
      const docSpy = vi.fn();
      document.addEventListener("keydown", docSpy);
      try {
        fireEvent.click(trigger);
        fireEvent.keyDown(screen.getByTestId("menu"), { key: "Escape" });
        expect(screen.queryByTestId("menu")).toBeNull();
        expect(docSpy).not.toHaveBeenCalled();
        // Escape on the trigger while open is equally contained.
        fireEvent.click(trigger);
        fireEvent.keyDown(trigger, { key: "Escape" });
        expect(screen.queryByTestId("menu")).toBeNull();
        expect(docSpy).not.toHaveBeenCalled();
      } finally {
        document.removeEventListener("keydown", docSpy);
      }
    });

    it("an outside pointerdown closes WITHOUT stealing focus back to the trigger", () => {
      const trigger = setup();
      fireEvent.click(trigger);
      fireEvent.pointerDown(document.body);
      expect(screen.queryByTestId("menu")).toBeNull();
      expect(document.activeElement).not.toBe(trigger);
    });
  });
});
