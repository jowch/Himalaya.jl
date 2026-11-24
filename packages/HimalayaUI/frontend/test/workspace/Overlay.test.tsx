import { describe, it, expect, vi } from "vitest";
import { startTransition, useState } from "react";
import { flushSync } from "react-dom";
import { render, screen, fireEvent } from "@testing-library/react";
import { Overlay } from "../../src/ui/Overlay";

function Harness({ onClose = () => undefined, empty = false }: { onClose?: () => void; empty?: boolean }): JSX.Element {
  const [open, setOpen] = useState(false);
  return (
    <>
      <button type="button" onClick={() => setOpen(true)}>Open</button>
      {open && (
        <Overlay
          label="Test"
          onClose={() => {
            onClose();
            setOpen(false);
          }}
        >
          {!empty && (
            <>
              <input aria-label="First" />
              <button type="button">Middle</button>
              <button type="button">Last</button>
            </>
          )}
        </Overlay>
      )}
    </>
  );
}

function open(props: Parameters<typeof Harness>[0] = {}): HTMLElement {
  render(<Harness {...props} />);
  const opener = screen.getByRole("button", { name: "Open" });
  opener.focus();
  fireEvent.click(opener);
  return opener;
}

describe("Overlay", () => {
  it("moves focus into the dialog and gives it back to the opener on close", () => {
    const opener = open();
    expect(document.activeElement).toBe(screen.getByLabelText("First"));
    fireEvent.keyDown(screen.getByRole("dialog"), { key: "Escape" });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(document.activeElement).toBe(opener);
  });

  it("wraps Tab from the last item to the first, and Shift Tab the other way", () => {
    open();
    screen.getByRole("button", { name: "Last" }).focus();
    fireEvent.keyDown(document.activeElement!, { key: "Tab" });
    expect(document.activeElement).toBe(screen.getByLabelText("First"));
    fireEvent.keyDown(document.activeElement!, { key: "Tab", shiftKey: true });
    expect(document.activeElement).toBe(screen.getByRole("button", { name: "Last" }));
  });

  it("Shift Tab from the panel itself goes to the last item", () => {
    open();
    screen.getByRole("dialog").focus();
    fireEvent.keyDown(screen.getByRole("dialog"), { key: "Tab", shiftKey: true });
    expect(document.activeElement).toBe(screen.getByRole("button", { name: "Last" }));
  });

  it("focuses its body when it holds nothing focusable, so keys scroll it, and Tab keeps it there", () => {
    const opener = open({ empty: true });
    const body = screen.getByRole("region", { name: "Test" });
    expect(document.activeElement).toBe(body);
    expect(body.tabIndex).toBe(0);
    fireEvent.keyDown(body, { key: "Tab" });
    expect(document.activeElement).toBe(body);
    fireEvent.keyDown(body, { key: "Escape" });
    expect(document.activeElement).toBe(opener);
  });

  it("keeps its body out of the Tab order when it holds something focusable", () => {
    open();
    expect(screen.getByTestId("overlay-body").tabIndex).toBe(-1);
    expect(screen.queryByRole("region")).toBeNull();
  });

  it("with a heading, names the dialog by its visible title once and offers a Close button", () => {
    const onClose = vi.fn();
    render(<Overlay label="Things" heading onClose={onClose}><p>x</p></Overlay>);
    const dialog = screen.getByRole("dialog", { name: "Things" });
    expect(dialog.getAttribute("aria-label")).toBeNull();
    expect(screen.getAllByText("Things")).toHaveLength(1);
    expect(screen.getByRole("region", { name: "Things" })).toBe(document.activeElement);
    fireEvent.click(screen.getByRole("button", { name: "Close" }));
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("closes on a press on the backdrop but not on one inside the panel", () => {
    const onClose = vi.fn();
    open({ onClose });
    fireEvent.mouseDown(screen.getByRole("button", { name: "Middle" }));
    fireEvent.mouseDown(screen.getByRole("dialog"));
    expect(onClose).not.toHaveBeenCalled();
    fireEvent.mouseDown(screen.getByRole("dialog").parentElement!);
    expect(onClose).toHaveBeenCalledTimes(1);
  });

  it("survives a scroll update that is still waiting when an urgent one lands, as when Home is pressed mid-scroll", () => {
    render(<Overlay label="Long" onClose={() => undefined}><p>x</p></Overlay>);
    const body = screen.getByTestId("overlay-body");
    Object.defineProperty(body, "scrollHeight", { configurable: true, value: 1000 });
    Object.defineProperty(body, "clientHeight", { configurable: true, value: 100 });
    const errors = vi.spyOn(console, "error").mockImplementation(() => undefined);
    try {
      startTransition(() => { body.dispatchEvent(new Event("scroll")); });
      expect(() => flushSync(() => { window.dispatchEvent(new Event("resize")); })).not.toThrow();
    } finally {
      errors.mockRestore();
    }
    expect(screen.getByRole("dialog", { name: "Long" })).toBeInTheDocument();
    expect(screen.getByTestId("overlay-more")).toBeInTheDocument();
  });
});
