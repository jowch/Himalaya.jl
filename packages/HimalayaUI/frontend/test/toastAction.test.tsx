import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { ToastContainer } from "../src/print/ui/Toast";
import { showToast } from "../src/lib/toast";

describe("imperative toast action slot", () => {
  it("renders an action button and fires onClick", () => {
    render(<ToastContainer />); // singleton — registers showToast on mount
    const onUndo = vi.fn();
    act(() => { showToast("Moved a1.tif", "info", { label: "Undo", onClick: onUndo }); });
    const btn = screen.getByRole("button", { name: "Undo" });
    fireEvent.click(btn);
    expect(onUndo).toHaveBeenCalledTimes(1);
  });

  it("renders no action button when no action is passed", () => {
    render(<ToastContainer />);
    act(() => { showToast("Saved", "success"); });
    expect(screen.queryByRole("button", { name: "Undo" })).toBeNull();
  });
});
