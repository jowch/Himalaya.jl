import { describe, it, expect, vi, afterEach } from "vitest";
import { act, fireEvent, render, screen } from "@testing-library/react";
import { NAV_PEEK_DELAY_MS, NAV_PEEK_LEAVE_MS, WorkspaceLayout } from "../../src/ui/WorkspaceLayout";

function layout(navHidden: boolean, onShowNav = vi.fn()) {
  return render(
    <WorkspaceLayout header={<span>head</span>} nav={<input aria-label="Filter samples" />} main={<p>plot</p>} tools={<p>tools</p>} navHidden={navHidden} onShowNav={onShowNav} />,
  );
}

const nav = (): HTMLElement => screen.getByRole("navigation", { name: "Samples" });

afterEach(() => {
  vi.useRealTimers();
});

describe("WorkspaceLayout hidden sample list", () => {
  it("shows the list over the plot only after the pointer rests on the handle, so passing over the edge does nothing", () => {
    vi.useFakeTimers();
    layout(true);
    const handle = screen.getByRole("button", { name: "Show the sample list" });
    fireEvent.mouseEnter(handle);
    act(() => void vi.advanceTimersByTime(NAV_PEEK_DELAY_MS - 50));
    fireEvent.mouseLeave(handle);
    act(() => void vi.advanceTimersByTime(NAV_PEEK_DELAY_MS));
    expect(nav().dataset.peek).toBeUndefined();
    fireEvent.mouseEnter(handle);
    act(() => void vi.advanceTimersByTime(NAV_PEEK_DELAY_MS));
    expect(nav().dataset.peek).toBe("true");
  });

  it("keeps the shown list while the pointer is over it, and hides it shortly after the pointer leaves", () => {
    vi.useFakeTimers();
    layout(true);
    const handle = screen.getByRole("button", { name: "Show the sample list" });
    fireEvent.mouseEnter(handle);
    act(() => void vi.advanceTimersByTime(NAV_PEEK_DELAY_MS));
    fireEvent.mouseLeave(handle);
    fireEvent.mouseEnter(nav());
    act(() => void vi.advanceTimersByTime(NAV_PEEK_LEAVE_MS * 2));
    expect(nav().dataset.peek).toBe("true");
    fireEvent.mouseLeave(nav());
    act(() => void vi.advanceTimersByTime(NAV_PEEK_LEAVE_MS));
    expect(nav().dataset.peek).toBeUndefined();
  });

  it("keeps the list shown while focus is inside it", () => {
    layout(true);
    const field = screen.getByRole("textbox", { name: "Filter samples" });
    fireEvent.focus(field);
    expect(nav().dataset.peek).toBe("true");
    fireEvent.blur(field);
    expect(nav().dataset.peek).toBeUndefined();
  });

  it("clicking the handle asks to keep the list open", () => {
    const onShowNav = vi.fn();
    layout(true, onShowNav);
    fireEvent.click(screen.getByRole("button", { name: "Show the sample list" }));
    expect(onShowNav).toHaveBeenCalledTimes(1);
  });

  it("has no handle while the list is shown", () => {
    layout(false);
    expect(screen.queryByRole("button", { name: "Show the sample list" })).toBeNull();
    expect(nav().dataset.peek).toBeUndefined();
  });
});
