import { render, screen } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import userEvent from "@testing-library/user-event";
import { IconButton } from "../../src/print/ui/IconButton";

describe("IconButton", () => {
  it("uses label as the accessible name", () => {
    render(<IconButton label="Close notes" dismiss />);
    expect(screen.getByRole("button", { name: "Close notes" })).toBeInTheDocument();
  });

  it("dismiss renders the canonical × glyph (U+00D7)", () => {
    render(<IconButton label="Dismiss" dismiss />);
    expect(screen.getByRole("button", { name: "Dismiss" })).toHaveTextContent("×");
  });

  it("renders children instead of the glyph when not dismiss", () => {
    render(<IconButton label="Move up">{"▲"}</IconButton>);
    expect(screen.getByRole("button", { name: "Move up" })).toHaveTextContent("▲");
  });

  it("defaults tone to ghost and exposes data-tone", () => {
    render(<IconButton label="x" dismiss />);
    expect(screen.getByRole("button", { name: "x" })).toHaveAttribute("data-tone", "ghost");
  });

  it("exposes data-tone for accent and danger", () => {
    const { rerender } = render(<IconButton label="remove" dismiss tone="accent" />);
    expect(screen.getByRole("button", { name: "remove" })).toHaveAttribute("data-tone", "accent");
    rerender(<IconButton label="delete" dismiss tone="danger" />);
    expect(screen.getByRole("button", { name: "delete" })).toHaveAttribute("data-tone", "danger");
  });

  it("shows a Tooltip caption from the label on hover; a passed title overrides it (item 7)", async () => {
    const user = userEvent.setup();
    const { rerender } = render(<IconButton label="Remove from recipe" dismiss />);
    // The native browser `title` is gone — the caption lives in our Tooltip,
    // which surfaces on hover/focus.
    const btn = screen.getByRole("button", { name: "Remove from recipe" });
    expect(btn).not.toHaveAttribute("title");
    await user.hover(btn);
    expect(screen.getByRole("tooltip")).toHaveTextContent("Remove from recipe");
    await user.unhover(btn);

    rerender(<IconButton label="Remove from recipe" dismiss title="Custom" />);
    await user.hover(screen.getByRole("button", { name: "Remove from recipe" }));
    expect(screen.getByRole("tooltip")).toHaveTextContent("Custom");
  });

  it("defaults type to button (does not submit forms)", () => {
    render(<IconButton label="x" dismiss />);
    expect(screen.getByRole("button", { name: "x" })).toHaveAttribute("type", "button");
  });

  it("is disabled and does not fire onClick when disabled", async () => {
    const onClick = vi.fn();
    render(<IconButton label="Move up" disabled onClick={onClick}>{"▲"}</IconButton>);
    const btn = screen.getByRole("button", { name: "Move up" });
    expect(btn).toBeDisabled();
    await userEvent.click(btn);
    expect(onClick).not.toHaveBeenCalled();
  });

  it("fires onClick when enabled", async () => {
    const onClick = vi.fn();
    render(<IconButton label="Dismiss" dismiss onClick={onClick} />);
    await userEvent.click(screen.getByRole("button", { name: "Dismiss" }));
    expect(onClick).toHaveBeenCalledTimes(1);
  });

  it("forwards data-testid and placement className to the button", () => {
    render(<IconButton label="x" dismiss data-testid="my-id" className="shrink-0" />);
    const btn = screen.getByTestId("my-id");
    expect(btn).toHaveClass("shrink-0");
  });
});
