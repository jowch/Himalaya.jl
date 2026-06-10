import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { FlagButton } from "../../src/print/ui/FlagButton";

describe("<FlagButton>", () => {
  it("exposes purpose + value in the accessible name", () => {
    render(<FlagButton value="1 : 0.25" />);
    // Anchored: the accname is EXACTLY purpose + value, nothing more.
    const btn = screen.getByRole("button", { name: /^skip this read: 1 : 0\.25$/i });
    expect(btn).toBe(screen.getByTestId("flag-button"));
    expect(btn).toHaveTextContent("1 : 0.25");
  });

  it("keeps the same accessible name in the flagged state (state lives in aria-pressed)", () => {
    render(<FlagButton value="1 : 0.25" flagged />);
    // Anchored so this FAILS if the "▸ skipped" caption ever joins the name.
    expect(
      screen.getByRole("button", { name: /^skip this read: 1 : 0\.25$/i }),
    ).toBeInTheDocument();
  });

  it("defaults to the ok state: not pressed, skip tooltip, no caption", () => {
    render(<FlagButton value="1 : 1" />);
    const btn = screen.getByTestId("flag-button");
    expect(btn.getAttribute("data-state")).toBe("ok");
    expect(btn.getAttribute("aria-pressed")).toBe("false");
    expect(btn.getAttribute("title")).toBe("Skip this read");
    expect(screen.queryByText(/skipped/i)).toBeNull();
  });

  it("renders the flagged state: pressed, restore tooltip, skipped caption", () => {
    render(<FlagButton value="1 : 0" flagged />);
    const btn = screen.getByTestId("flag-button");
    expect(btn.getAttribute("data-state")).toBe("flagged");
    expect(btn.getAttribute("aria-pressed")).toBe("true");
    expect(btn.getAttribute("title")).toBe("Restore this read");
    expect(screen.getByText(/▸ skipped/i)).toBeInTheDocument();
  });

  it("never speaks the retired uncertain-parse vocabulary", () => {
    const { unmount } = render(<FlagButton value="1 : 0" flagged />);
    expect(screen.queryByText(/check the read/i)).toBeNull();
    unmount();
    render(<FlagButton value="1 : 0" />);
    expect(screen.getByTestId("flag-button").getAttribute("title")).not.toMatch(/re-open/i);
  });

  it("fires onClick", () => {
    const onClick = vi.fn();
    render(<FlagButton value="1 : 0" flagged onClick={onClick} />);
    fireEvent.click(screen.getByTestId("flag-button"));
    expect(onClick).toHaveBeenCalledOnce();
  });
});
