import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { useShortcuts } from "../../src/print/shell/useShortcuts";

function Harness({ bindings, enabled }: { bindings: Parameters<typeof useShortcuts>[0]; enabled?: boolean }) {
  useShortcuts(bindings, enabled);
  return <input data-testid="field" />;
}

describe("useShortcuts", () => {
  it("dispatches a window keydown to the matching action handler and preventDefaults", () => {
    const drop = vi.fn();
    render(<Harness bindings={{ drop }} />);
    const e = new KeyboardEvent("keydown", { key: "x", cancelable: true });
    window.dispatchEvent(e);
    expect(drop).toHaveBeenCalledTimes(1);
    expect(e.defaultPrevented).toBe(true);
  });

  it("ignores keys with no binding (no preventDefault)", () => {
    const drop = vi.fn();
    render(<Harness bindings={{ drop }} />);
    const e = new KeyboardEvent("keydown", { key: "k", cancelable: true });
    window.dispatchEvent(e);
    expect(drop).not.toHaveBeenCalled();
    expect(e.defaultPrevented).toBe(false);
  });

  it("respects suppressGlobalKeys: a keystroke from an INPUT does not fire page shortcuts", () => {
    const drop = vi.fn();
    const { getByTestId } = render(<Harness bindings={{ drop }} />);
    fireEvent.keyDown(getByTestId("field"), { key: "x" });
    expect(drop).not.toHaveBeenCalled();
  });

  it("unbinds on unmount", () => {
    const drop = vi.fn();
    const { unmount } = render(<Harness bindings={{ drop }} />);
    unmount();
    window.dispatchEvent(new KeyboardEvent("keydown", { key: "x" }));
    expect(drop).not.toHaveBeenCalled();
  });

  it("does nothing when disabled", () => {
    const drop = vi.fn();
    render(<Harness bindings={{ drop }} enabled={false} />);
    window.dispatchEvent(new KeyboardEvent("keydown", { key: "x" }));
    expect(drop).not.toHaveBeenCalled();
  });

  it("a handler returning false DECLINES the key: no preventDefault, event keeps propagating", () => {
    const decline = vi.fn(() => false as const);
    render(<Harness bindings={{ dismiss: decline }} />);
    const e = new KeyboardEvent("keydown", { key: "Escape", cancelable: true });
    window.dispatchEvent(e);
    expect(decline).toHaveBeenCalledTimes(1);
    expect(e.defaultPrevented).toBe(false); // declined → still un-prevented for the next rung
  });

  it("reads the latest handler without re-binding (stable across renders)", () => {
    const first = vi.fn();
    const second = vi.fn();
    const { rerender } = render(<Harness bindings={{ drop: first }} />);
    rerender(<Harness bindings={{ drop: second }} />);
    window.dispatchEvent(new KeyboardEvent("keydown", { key: "x" }));
    expect(first).not.toHaveBeenCalled();
    expect(second).toHaveBeenCalledTimes(1);
  });
});
