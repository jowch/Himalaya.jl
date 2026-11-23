import { describe, it, expect, vi } from "vitest";
import { useState } from "react";
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
    const body = screen.getByTestId("overlay-body");
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
  });

  it("with a heading, names the dialog by its visible title once and offers a Close button", () => {
    const onClose = vi.fn();
    render(<Overlay label="Things" heading onClose={onClose}><p>x</p></Overlay>);
    const dialog = screen.getByRole("dialog", { name: "Things" });
    expect(dialog.getAttribute("aria-label")).toBeNull();
    expect(screen.getAllByText("Things")).toHaveLength(1);
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
});
