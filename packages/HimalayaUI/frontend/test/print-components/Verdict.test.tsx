import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Verdict } from "../../src/print/components/Verdict";

// SA-SCREENED: the verdict is tri-state. `dropped` = rejected, `kept` =
// explicitly accepted, neither = unscreened. The old binary read "Kept" for a
// null status, which is the lie the Keep verb fixes.

describe("<Verdict> unscreened state (neither dropped nor kept)", () => {
  it('shows "Unscreened" state text, not "Kept"', () => {
    render(<Verdict dropped={false} />);
    expect(screen.getByText("Unscreened")).toBeInTheDocument();
    expect(screen.queryByText("Kept")).not.toBeInTheDocument();
  });

  it("shows the honest unscreened hint", () => {
    render(<Verdict dropped={false} />);
    expect(
      screen.getByText("Keep or drop to screen this exposure.")
    ).toBeInTheDocument();
  });

  it("offers Keep and Drop", () => {
    const onToggle = vi.fn();
    const onToggleKeep = vi.fn();
    render(<Verdict dropped={false} onToggle={onToggle} onToggleKeep={onToggleKeep} />);
    fireEvent.click(screen.getByRole("button", { name: "Keep" }));
    expect(onToggleKeep).toHaveBeenCalledTimes(1);
    fireEvent.click(screen.getByRole("button", { name: "Drop" }));
    expect(onToggle).toHaveBeenCalledTimes(1);
  });

  it("status Dot has data-tone=neutral when unscreened", () => {
    const { container } = render(<Verdict dropped={false} />);
    expect(container.querySelector("[data-tone='neutral']")).toBeInTheDocument();
  });
});

describe("<Verdict> kept state (kept=true)", () => {
  it('shows "Kept" state text', () => {
    render(<Verdict dropped={false} kept />);
    expect(screen.getByText("Kept")).toBeInTheDocument();
  });

  it("offers Restore (the keep toggle) and Drop", () => {
    const onToggle = vi.fn();
    const onToggleKeep = vi.fn();
    render(<Verdict dropped={false} kept onToggle={onToggle} onToggleKeep={onToggleKeep} />);
    // Restore clears the accepted call back to unscreened → keep toggle.
    fireEvent.click(screen.getByRole("button", { name: "Restore" }));
    expect(onToggleKeep).toHaveBeenCalledTimes(1);
    // Drop on a kept frame is direct: last verb wins, no trip through null.
    fireEvent.click(screen.getByRole("button", { name: "Drop" }));
    expect(onToggle).toHaveBeenCalledTimes(1);
  });

  it("status Dot has data-tone=success when kept", () => {
    const { container } = render(<Verdict dropped={false} kept />);
    expect(container.querySelector("[data-tone='success']")).toBeInTheDocument();
  });
});

describe("<Verdict> dropped state (dropped=true)", () => {
  it('shows "Dropped" state text', () => {
    render(<Verdict dropped={true} />);
    expect(screen.getByText("Dropped")).toBeInTheDocument();
  });

  it("offers Keep and Restore (the drop toggle)", () => {
    const onToggle = vi.fn();
    const onToggleKeep = vi.fn();
    render(<Verdict dropped={true} onToggle={onToggle} onToggleKeep={onToggleKeep} />);
    // Restore clears the rejected call back to unscreened → drop toggle.
    fireEvent.click(screen.getByRole("button", { name: "Restore" }));
    expect(onToggle).toHaveBeenCalledTimes(1);
    // Keep on a dropped frame is direct: last verb wins.
    fireEvent.click(screen.getByRole("button", { name: "Keep" }));
    expect(onToggleKeep).toHaveBeenCalledTimes(1);
  });

  it("status Dot has data-tone=accent when dropped", () => {
    const { container } = render(<Verdict dropped={true} />);
    expect(container.querySelector("[data-tone='accent']")).toBeInTheDocument();
  });
});

describe("<Verdict> hint override", () => {
  it("renders a custom hint when provided", () => {
    render(<Verdict dropped={false} hint="Custom hint text." />);
    expect(screen.getByText("Custom hint text.")).toBeInTheDocument();
    expect(
      screen.queryByText("Keep or drop to screen this exposure.")
    ).not.toBeInTheDocument();
  });
});

describe("<Verdict> no handlers", () => {
  it("does not throw when clicked without onToggle/onToggleKeep", () => {
    render(<Verdict dropped={false} />);
    expect(() => {
      fireEvent.click(screen.getByRole("button", { name: "Drop" }));
      fireEvent.click(screen.getByRole("button", { name: "Keep" }));
    }).not.toThrow();
  });
});
