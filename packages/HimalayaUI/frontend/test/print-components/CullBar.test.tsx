import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { CullBar } from "../../src/print/components/CullBar";

describe("CullBar", () => {
  it("shows the selected-frame count", () => {
    render(<CullBar count={3} show />);
    const bar = screen.getByTestId("cull-bar");
    expect(bar).toHaveAttribute("data-show", "true");
    expect(bar.textContent).toContain("3");
    expect(bar.textContent).toContain("frames selected");
  });

  it("is hidden (data-show=false) when show is omitted", () => {
    render(<CullBar count={0} />);
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");
  });

  it("wires reject, restore, and clear", async () => {
    const onReject = vi.fn(), onRestore = vi.fn(), onClear = vi.fn();
    render(<CullBar count={2} show onReject={onReject} onRestore={onRestore} onClear={onClear} />);
    await userEvent.click(screen.getByRole("button", { name: /drop/i }));
    await userEvent.click(screen.getByRole("button", { name: /restore/i }));
    await userEvent.click(screen.getByRole("button", { name: /clear/i }));
    expect(onReject).toHaveBeenCalledOnce();
    expect(onRestore).toHaveBeenCalledOnce();
    expect(onClear).toHaveBeenCalledOnce();
  });

  it("uses the accent reject button and ghostInverse for restore/clear", () => {
    render(<CullBar count={1} show />);
    expect(screen.getByRole("button", { name: /drop/i }).getAttribute("data-variant")).toBe("accent");
    expect(screen.getByRole("button", { name: /restore/i }).getAttribute("data-variant")).toBe("ghostInverse");
    expect(screen.getByRole("button", { name: /clear/i }).getAttribute("data-variant")).toBe("ghostInverse");
  });

  it("removes the bar from the a11y tree and tab order when hidden", () => {
    render(<CullBar count={0} />);
    const bar = screen.getByTestId("cull-bar");
    expect(bar).toHaveAttribute("aria-hidden", "true");
    screen.getAllByRole("button", { hidden: true }).forEach((b) =>
      expect(b).toHaveAttribute("tabindex", "-1"),
    );
  });
});
