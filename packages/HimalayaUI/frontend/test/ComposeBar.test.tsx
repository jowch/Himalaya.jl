// test/ComposeBar.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { ComposeBar } from "../src/print/ui/ComposeBar";

describe("ComposeBar", () => {
  it("is hidden (inert) when show=false", () => {
    // JSDOM note: assert the inert ATTRIBUTE only — JSDOM does not implement
    // inert focus/AT semantics.
    render(<ComposeBar count={0} show={false} />);
    const bar = screen.getByTestId("compose-bar");
    expect(bar).toHaveAttribute("data-show", "false");
    expect(bar).toHaveAttribute("inert");
    expect(bar).not.toHaveAttribute("aria-hidden");
  });

  it("is visible (not inert) when show=true", () => {
    render(<ComposeBar count={2} show />);
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");
    expect(screen.getByTestId("compose-bar")).not.toHaveAttribute("inert");
  });

  it("displays the sample count", () => {
    render(<ComposeBar count={3} show />);
    expect(screen.getByTestId("compose-bar")).toHaveTextContent("3");
  });

  it("calls onNewSeries when + New series is clicked", async () => {
    const user = userEvent.setup();
    const onNewSeries = vi.fn();
    render(<ComposeBar count={2} show onNewSeries={onNewSeries} />);
    await user.click(screen.getByRole("button", { name: /new series/i }));
    expect(onNewSeries).toHaveBeenCalledTimes(1);
  });

  it("calls onClear when Clear is clicked", async () => {
    const user = userEvent.setup();
    const onClear = vi.fn();
    render(<ComposeBar count={1} show onClear={onClear} />);
    await user.click(screen.getByRole("button", { name: /clear/i }));
    expect(onClear).toHaveBeenCalledTimes(1);
  });

  it("buttons carry no per-button tabindex (inert owns the tab order)", () => {
    render(<ComposeBar count={0} show={false} />);
    screen.getAllByRole("button").forEach((b) =>
      expect(b).not.toHaveAttribute("tabindex"),
    );
  });
});
