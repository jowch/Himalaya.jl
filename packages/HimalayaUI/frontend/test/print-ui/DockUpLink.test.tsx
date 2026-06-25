import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { DockUpLink } from "../../src/print/ui/DockUpLink";

describe("<DockUpLink>", () => {
  it("renders a button with the chevron-prefixed label and the testid contract", () => {
    render(<DockUpLink label="Corpus" onClick={() => {}} />);
    const el = screen.getByTestId("dock-up-link");
    expect(el.tagName).toBe("BUTTON");
    expect(el).toHaveTextContent("‹ Corpus");
  });

  it("fires onClick", () => {
    const onClick = vi.fn();
    render(<DockUpLink label="All series" onClick={onClick} />);
    fireEvent.click(screen.getByTestId("dock-up-link"));
    expect(onClick).toHaveBeenCalledTimes(1);
  });

  it("renders an anchor when href is set, and preventDefaults navigation", () => {
    const onClick = vi.fn();
    render(<DockUpLink label="Experiments" href="/experiments" onClick={onClick} />);
    const el = screen.getByTestId("dock-up-link");
    expect(el.tagName).toBe("A");
    expect(el).toHaveAttribute("href", "/experiments");
    const ev = new MouseEvent("click", { bubbles: true, cancelable: true });
    el.dispatchEvent(ev);
    expect(ev.defaultPrevented).toBe(true);
    expect(onClick).toHaveBeenCalledTimes(1);
  });

  it("appends placement-only className", () => {
    render(<DockUpLink label="Corpus" onClick={() => {}} className="mr-1" />);
    expect(screen.getByTestId("dock-up-link")).toHaveClass("mr-1");
  });
});
