import { describe, it, expect, afterEach, vi } from "vitest";
import { act, render, screen } from "@testing-library/react";
import { LiveRegion } from "../src/print/ui/LiveRegion";
import { announce, setAnnounceImpl } from "../src/lib/announce";

afterEach(() => {
  setAnnounceImpl(null);
});

describe("announce indirection", () => {
  it("falls back to a no-op (console.debug) before any region mounts", () => {
    const spy = vi.spyOn(console, "debug").mockImplementation(() => {});
    // Must not throw with no impl registered.
    expect(() => announce("untethered")).not.toThrow();
    spy.mockRestore();
  });

  it("restores the fallback when setAnnounceImpl(null) is called", () => {
    const impl = vi.fn();
    setAnnounceImpl(impl);
    announce("one");
    expect(impl).toHaveBeenCalledWith("one", undefined);
    setAnnounceImpl(null);
    // After clearing, the registered impl must not be called again.
    announce("two");
    expect(impl).toHaveBeenCalledTimes(1);
  });
});

describe("LiveRegion", () => {
  it("mounts a visually-hidden polite and assertive region", () => {
    render(<LiveRegion />);
    const polite = screen.getByTestId("live-region-polite");
    const assertive = screen.getByTestId("live-region-assertive");
    expect(polite).toHaveAttribute("aria-live", "polite");
    expect(assertive).toHaveAttribute("aria-live", "assertive");
    // Visually-hidden: uses the sr-only utility (semantics, not a class string
    // assertion of appearance — sr-only IS the project hiding mechanism).
    expect(polite.className).toContain("sr-only");
    expect(assertive.className).toContain("sr-only");
  });

  it("writes a polite message into the polite node", () => {
    render(<LiveRegion />);
    act(() => {
      announce("Peak added");
    });
    expect(screen.getByTestId("live-region-polite")).toHaveTextContent("Peak added");
    expect(screen.getByTestId("live-region-assertive")).toHaveTextContent("");
  });

  it("writes an assertive message into the assertive node", () => {
    render(<LiveRegion />);
    act(() => {
      announce("Something failed", { assertive: true });
    });
    expect(screen.getByTestId("live-region-assertive")).toHaveTextContent("Something failed");
    expect(screen.getByTestId("live-region-polite")).toHaveTextContent("");
  });

  it("registers the announce impl on mount and clears it on unmount", () => {
    const { unmount } = render(<LiveRegion />);
    act(() => {
      announce("live");
    });
    expect(screen.getByTestId("live-region-polite")).toHaveTextContent("live");
    unmount();
    // After unmount, announce must not throw (fallback restored).
    const spy = vi.spyOn(console, "debug").mockImplementation(() => {});
    expect(() => announce("after unmount")).not.toThrow();
    spy.mockRestore();
  });
});

describe("LiveRegion interleaved re-speak", () => {
  it("re-speaks an identical polite message after an interleaved assertive one", () => {
    render(<LiveRegion />);
    const polite = () => screen.getByTestId("live-region-polite");
    act(() => {
      announce("A");
    });
    const firstText = polite().textContent;
    expect(firstText).toContain("A");
    // Interleaved assertive write — must NOT consume the polite region's flip.
    act(() => {
      announce("X", { assertive: true });
    });
    // Same polite string again: the polite node's textContent MUST change so the
    // SR re-speaks. A shared flip would land back on the exact prior text → miss.
    act(() => {
      announce("A");
    });
    expect(polite().textContent).not.toBe(firstText);
    expect(polite().textContent).toContain("A");
  });
});

describe("LiveRegion barrel export", () => {
  it("LiveRegion is exported from the ui barrel", async () => {
    const mod = await import("../src/print/ui");
    expect(typeof mod.LiveRegion).toBe("function");
  });
});
