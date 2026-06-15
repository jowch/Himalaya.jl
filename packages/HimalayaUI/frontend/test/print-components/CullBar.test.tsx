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

  it("uses the singular grain for a single frame", () => {
    render(<CullBar count={1} show />);
    const bar = screen.getByTestId("cull-bar");
    expect(bar.textContent).toContain("frame selected");
    expect(bar.textContent).not.toContain("frames selected");
  });

  it("stays count-only when sampleCount is omitted", () => {
    render(<CullBar count={3} show />);
    const bar = screen.getByTestId("cull-bar");
    expect(bar.textContent).toContain("frames selected");
    expect(bar.textContent).not.toContain("across");
  });

  it("stays count-only when the selection sits within one sample (sampleCount=1)", () => {
    // No noise when nothing is hidden: a single-sample selection reads
    // exactly "N frames selected" — and never "across 1 samples".
    render(<CullBar count={3} sampleCount={1} show />);
    const bar = screen.getByTestId("cull-bar");
    expect(bar.textContent).toContain("frames selected");
    expect(bar.textContent).not.toContain("across");
  });

  it("discloses the spread when the selection spans samples (SA-F3)", () => {
    render(<CullBar count={5} sampleCount={2} show />);
    const bar = screen.getByTestId("cull-bar");
    expect(bar.textContent).toContain("5");
    expect(bar.textContent).toContain("frames selected across");
    expect(bar.textContent).toContain("2");
    expect(bar.textContent).toContain("samples");
  });

  it("keeps both grains correct: singular frame, plural samples", () => {
    render(<CullBar count={1} sampleCount={2} show />);
    const bar = screen.getByTestId("cull-bar");
    expect(bar.textContent).toContain("frame selected across");
    expect(bar.textContent).not.toContain("frames selected");
    // The suffix only ever appears for a spread (>1), so "samples" is plural.
    expect(bar.textContent).toContain("samples");
    expect(bar.textContent).not.toContain("1 samples");
  });

  it("is hidden (data-show=false) when show is omitted", () => {
    render(<CullBar count={0} />);
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");
  });

  it("wires reject, keep, restore, and clear", async () => {
    const onReject = vi.fn(), onKeep = vi.fn(), onRestore = vi.fn(), onClear = vi.fn();
    render(<CullBar count={2} show onReject={onReject} onKeep={onKeep} onRestore={onRestore} onClear={onClear} />);
    await userEvent.click(screen.getByRole("button", { name: /drop/i }));
    await userEvent.click(screen.getByRole("button", { name: /keep/i }));
    await userEvent.click(screen.getByRole("button", { name: /restore/i }));
    await userEvent.click(screen.getByRole("button", { name: /clear/i }));
    expect(onReject).toHaveBeenCalledOnce();
    expect(onKeep).toHaveBeenCalledOnce();
    expect(onRestore).toHaveBeenCalledOnce();
    expect(onClear).toHaveBeenCalledOnce();
  });

  it("Keep carries the K key hint (matching Drop's X hint)", () => {
    render(<CullBar count={2} show />);
    const keep = screen.getByRole("button", { name: /keep/i });
    expect(keep.textContent).toContain("K");
    const drop = screen.getByRole("button", { name: /drop/i });
    expect(drop.textContent).toContain("X");
  });

  it("uses the accent reject button, a constructive (non-accent) Keep, and ghostInverse for restore/clear", () => {
    render(<CullBar count={1} show />);
    expect(screen.getByRole("button", { name: /drop/i }).getAttribute("data-variant")).toBe("accent");
    expect(screen.getByRole("button", { name: /keep/i }).getAttribute("data-variant")).toBe("success");
    expect(screen.getByRole("button", { name: /restore/i }).getAttribute("data-variant")).toBe("ghostInverse");
    expect(screen.getByRole("button", { name: /clear/i }).getAttribute("data-variant")).toBe("ghostInverse");
  });

  it("is inert when hidden (removes the subtree from the a11y tree and tab order)", () => {
    // JSDOM note: assert the inert ATTRIBUTE only — JSDOM does not implement
    // inert focus/AT semantics. Real-browser behavior is pinned by the
    // corpus-culling e2e (focus leaves the bar, no aria-hidden warning).
    render(<CullBar count={0} />);
    const bar = screen.getByTestId("cull-bar");
    expect(bar).toHaveAttribute("inert");
    expect(bar).not.toHaveAttribute("aria-hidden");
    screen.getAllByRole("button").forEach((b) =>
      expect(b).not.toHaveAttribute("tabindex"),
    );
  });

  it("is not inert when shown", () => {
    render(<CullBar count={2} show />);
    expect(screen.getByTestId("cull-bar")).not.toHaveAttribute("inert");
  });

  it("becomes inert when show flips false after a click empties the selection", () => {
    const { rerender } = render(<CullBar count={2} show />);
    rerender(<CullBar count={0} show={false} />);
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("inert");
  });
});
