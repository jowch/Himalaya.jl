import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Verdict } from "../../src/print/components/Verdict";

// Binary screening verdict: `dropped` = status rejected, else kept (null).
// One verb only — Drop, a toggle. No Keep, no Restore, no Unscreened.

describe("<Verdict> kept state (dropped=false)", () => {
  it('shows "Kept", not "Unscreened"', () => {
    render(<Verdict dropped={false} />);
    expect(screen.getByText("Kept")).toBeInTheDocument();
    expect(screen.queryByText("Unscreened")).not.toBeInTheDocument();
  });

  it("status Dot has data-tone=success when kept", () => {
    const { container } = render(<Verdict dropped={false} />);
    expect(container.querySelector("[data-tone='success']")).toBeInTheDocument();
  });

  it("the Drop toggle is not pressed when kept", () => {
    render(<Verdict dropped={false} />);
    expect(screen.getByRole("button", { name: "Drop" })).toHaveAttribute("aria-pressed", "false");
  });
});

describe("<Verdict> dropped state (dropped=true)", () => {
  it('shows "Dropped"', () => {
    render(<Verdict dropped={true} />);
    expect(screen.getByText("Dropped")).toBeInTheDocument();
  });

  it("status Dot has data-tone=accent when dropped", () => {
    const { container } = render(<Verdict dropped={true} />);
    expect(container.querySelector("[data-tone='accent']")).toBeInTheDocument();
  });

  it("the Drop toggle is pressed when dropped", () => {
    render(<Verdict dropped={true} />);
    expect(screen.getByRole("button", { name: "Drop" })).toHaveAttribute("aria-pressed", "true");
  });
});

describe("<Verdict> the single Drop toggle", () => {
  it("offers exactly one verb (Drop) — never Keep or Restore", () => {
    render(<Verdict dropped={false} onToggle={vi.fn()} />);
    expect(screen.getByRole("button", { name: "Drop" })).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: /keep/i })).toBeNull();
    expect(screen.queryByRole("button", { name: /restore/i })).toBeNull();
  });

  it("Drop calls onToggle from either state", () => {
    const onToggle = vi.fn();
    const { rerender } = render(<Verdict dropped={false} onToggle={onToggle} />);
    fireEvent.click(screen.getByRole("button", { name: "Drop" }));
    rerender(<Verdict dropped={true} onToggle={onToggle} />);
    fireEvent.click(screen.getByRole("button", { name: "Drop" }));
    expect(onToggle).toHaveBeenCalledTimes(2);
  });

  it("does not throw when clicked without a handler", () => {
    render(<Verdict dropped={false} />);
    expect(() => fireEvent.click(screen.getByRole("button", { name: "Drop" }))).not.toThrow();
  });
});

describe("<Verdict> hint override", () => {
  it("renders a custom hint when provided", () => {
    render(<Verdict dropped={false} hint="Custom hint text." />);
    expect(screen.getByText("Custom hint text.")).toBeInTheDocument();
    expect(screen.queryByText("Drop to cull this frame.")).not.toBeInTheDocument();
  });
});
