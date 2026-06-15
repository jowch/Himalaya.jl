import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { useRef } from "react";
import { Popover } from "../../src/print/ui/Popover";

describe("<Popover>", () => {
  it("focuses initialFocusRef instead of the panel when provided (search-first)", () => {
    function Harness(): JSX.Element {
      const ref = useRef<HTMLInputElement>(null);
      return (
        <Popover trigger={<button>open</button>} label="Search" initialFocusRef={ref}>
          <input ref={ref} data-testid="inner-input" aria-label="search" />
        </Popover>
      );
    }
    render(<Harness />);
    fireEvent.click(screen.getByRole("button", { name: "open" }));
    expect(screen.getByTestId("inner-input")).toHaveFocus();
  });

  it("is closed by default — no panel, trigger collapsed", () => {
    render(
      <Popover trigger={<button>open</button>} label="Details">
        <button>inner</button>
      </Popover>,
    );
    expect(screen.queryByTestId("popover")).toBeNull();
    const trigger = screen.getByRole("button", { name: "open" });
    expect(trigger.getAttribute("aria-haspopup")).toBe("dialog");
    expect(trigger.getAttribute("aria-expanded")).toBe("false");
  });

  it("opens on trigger click; the labelled dialog panel appears", () => {
    render(
      <Popover trigger={<button>open</button>} label="Details">
        <button>inner</button>
      </Popover>,
    );
    fireEvent.click(screen.getByRole("button", { name: "open" }));
    const panel = screen.getByTestId("popover");
    expect(panel).toBeTruthy();
    expect(panel.getAttribute("aria-label")).toBe("Details");
    expect(screen.getByRole("button", { name: "open" }).getAttribute("aria-expanded")).toBe("true");
  });

  it("toggles closed on a second trigger click", () => {
    render(
      <Popover trigger={<button>open</button>} label="Details">
        <button>inner</button>
      </Popover>,
    );
    const trigger = screen.getByRole("button", { name: "open" });
    fireEvent.click(trigger);
    expect(screen.getByTestId("popover")).toBeTruthy();
    fireEvent.click(trigger);
    expect(screen.queryByTestId("popover")).toBeNull();
    expect(trigger.getAttribute("aria-expanded")).toBe("false");
  });

  it("renders the (interactive) content when open and it is clickable", () => {
    const onInner = vi.fn();
    render(
      <Popover trigger={<button>open</button>} label="Details">
        <button onClick={onInner}>inner</button>
      </Popover>,
    );
    fireEvent.click(screen.getByRole("button", { name: "open" }));
    const inner = screen.getByRole("button", { name: "inner" });
    fireEvent.click(inner);
    expect(onInner).toHaveBeenCalledTimes(1);
  });

  it("Escape closes the panel and returns focus to the trigger", () => {
    render(
      <Popover trigger={<button>open</button>} label="Details">
        <button>inner</button>
      </Popover>,
    );
    const trigger = screen.getByRole("button", { name: "open" });
    fireEvent.click(trigger);
    const panel = screen.getByTestId("popover");
    fireEvent.keyDown(panel, { key: "Escape" });
    expect(screen.queryByTestId("popover")).toBeNull();
    expect(document.activeElement).toBe(trigger);
  });

  it("a pointerdown outside the wrapper closes the panel", () => {
    render(
      <div>
        <Popover trigger={<button>open</button>} label="Details">
          <button>inner</button>
        </Popover>
        <span data-testid="outside">elsewhere</span>
      </div>,
    );
    fireEvent.click(screen.getByRole("button", { name: "open" }));
    expect(screen.getByTestId("popover")).toBeTruthy();
    fireEvent.pointerDown(screen.getByTestId("outside"));
    expect(screen.queryByTestId("popover")).toBeNull();
  });

  it("a pointerdown inside the panel keeps it open", () => {
    render(
      <Popover trigger={<button>open</button>} label="Details">
        <button>inner</button>
      </Popover>,
    );
    fireEvent.click(screen.getByRole("button", { name: "open" }));
    fireEvent.pointerDown(screen.getByRole("button", { name: "inner" }));
    expect(screen.getByTestId("popover")).toBeTruthy();
  });

  it("places the panel below by default and above when side='top'", () => {
    const { rerender } = render(
      <Popover trigger={<button>open</button>} label="Details">
        <button>inner</button>
      </Popover>,
    );
    fireEvent.click(screen.getByRole("button", { name: "open" }));
    expect(screen.getByTestId("popover").getAttribute("data-side")).toBe("bottom");
    rerender(
      <Popover trigger={<button>open</button>} label="Details" side="top">
        <button>inner</button>
      </Popover>,
    );
    expect(screen.getByTestId("popover").getAttribute("data-side")).toBe("top");
  });
});
