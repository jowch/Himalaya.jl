// test/ComposeBar.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { ComposeBar } from "../src/print/ui/ComposeBar";

describe("ComposeBar", () => {
  it("is hidden (aria-hidden) when show=false", () => {
    render(<ComposeBar count={0} show={false} />);
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("aria-hidden", "true");
  });

  it("is visible when show=true", () => {
    render(<ComposeBar count={2} show />);
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");
    expect(screen.getByTestId("compose-bar")).not.toHaveAttribute("aria-hidden", "true");
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

  it("+ New series button has tabIndex=-1 when hidden", () => {
    render(<ComposeBar count={0} show={false} />);
    expect(screen.getByRole("button", { name: /new series/i, hidden: true }))
      .toHaveAttribute("tabindex", "-1");
  });
});
