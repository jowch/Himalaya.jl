import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { act, render, screen, fireEvent } from "@testing-library/react";
import { ToastContainer } from "../src/print/ui/Toast";
import { showToast, setToastImpl } from "../src/lib/toast";

describe("Toast", () => {
  beforeEach(() => {
    vi.useFakeTimers();
  });
  afterEach(() => {
    vi.useRealTimers();
    setToastImpl(null);
  });

  it("renders a toast when showToast is called", () => {
    render(<ToastContainer />);
    act(() => {
      showToast("hello", "error");
    });
    expect(screen.getByTestId("toast")).toHaveTextContent("hello");
    expect(screen.getByTestId("toast")).toHaveAttribute("data-toast-kind", "error");
  });

  it("shows a severity word label per kind", () => {
    render(<ToastContainer />);
    const cases: Array<[import("../src/lib/toast").ToastKind, string]> = [
      ["error", "Error"],
      ["warning", "Warning"],
      ["success", "Success"],
      ["info", "Info"],
    ];
    for (const [kind, word] of cases) {
      act(() => {
        showToast(`msg-${kind}`, kind);
      });
      const toast = screen.getByTestId("toast");
      expect(toast).toHaveAttribute("data-toast-kind", kind);
      // visible severity word
      expect(toast).toHaveTextContent(word);
      // accessible status icon naming the severity (second channel, not hue)
      expect(toast.querySelector(`[aria-label="${word}"]`)).not.toBeNull();
      act(() => {
        fireEvent.click(screen.getByLabelText("Dismiss"));
      });
    }
  });

  it("announces error and warning toasts assertively (role=alert, aria-live=assertive)", () => {
    render(<ToastContainer />);
    for (const kind of ["error", "warning"] as const) {
      act(() => {
        showToast(`msg-${kind}`, kind);
      });
      const toast = screen.getByTestId("toast");
      expect(toast).toHaveAttribute("role", "alert");
      expect(toast).toHaveAttribute("aria-live", "assertive");
      act(() => {
        fireEvent.click(screen.getByLabelText("Dismiss"));
      });
    }
  });

  it("announces info and success toasts politely (role=status)", () => {
    render(<ToastContainer />);
    for (const kind of ["info", "success"] as const) {
      act(() => {
        showToast(`msg-${kind}`, kind);
      });
      const toast = screen.getByTestId("toast");
      expect(toast).toHaveAttribute("role", "status");
      // Positively lock the polite contract — an absent aria-live would also
      // satisfy `not assertive`, so assert the attribute is present + "polite".
      expect(toast).toHaveAttribute("aria-live", "polite");
      act(() => {
        fireEvent.click(screen.getByLabelText("Dismiss"));
      });
    }
  });

  it("uses a full hairline border, not a left-edge severity stripe", () => {
    render(<ToastContainer />);
    act(() => {
      showToast("hello", "error");
    });
    const toast = screen.getByTestId("toast");
    expect(toast.className).toContain("border-hair");
    expect(toast.className).not.toContain("border-l-4");
    expect(toast.className).not.toContain("border-error");
  });

  it("auto-dismisses error toast after 5000ms", () => {
    render(<ToastContainer />);
    act(() => {
      showToast("oops", "error");
    });
    expect(screen.getByTestId("toast")).toBeInTheDocument();
    act(() => {
      vi.advanceTimersByTime(4999);
    });
    expect(screen.queryByTestId("toast")).toBeInTheDocument();
    act(() => {
      vi.advanceTimersByTime(2);
    });
    expect(screen.queryByTestId("toast")).not.toBeInTheDocument();
  });

  it("auto-dismisses non-error toast after 3000ms", () => {
    render(<ToastContainer />);
    act(() => {
      showToast("ok", "info");
    });
    act(() => {
      vi.advanceTimersByTime(2999);
    });
    expect(screen.queryByTestId("toast")).toBeInTheDocument();
    act(() => {
      vi.advanceTimersByTime(2);
    });
    expect(screen.queryByTestId("toast")).not.toBeInTheDocument();
  });

  it("stacks multiple toasts", () => {
    render(<ToastContainer />);
    act(() => {
      showToast("first", "info");
      showToast("second", "warning");
      showToast("third", "success");
    });
    const toasts = screen.getAllByTestId("toast");
    expect(toasts).toHaveLength(3);
    expect(toasts[0]).toHaveTextContent("first");
    expect(toasts[1]).toHaveTextContent("second");
    expect(toasts[2]).toHaveTextContent("third");
  });

  it("close button dismisses immediately", () => {
    render(<ToastContainer />);
    act(() => {
      showToast("dismiss me", "info");
    });
    const btn = screen.getByLabelText("Dismiss");
    act(() => {
      fireEvent.click(btn);
    });
    expect(screen.queryByTestId("toast")).not.toBeInTheDocument();
  });

  it("clears auto-dismiss timer on manual close (no double dismiss)", () => {
    render(<ToastContainer />);
    act(() => {
      showToast("dismiss me", "info");
    });
    const btn = screen.getByLabelText("Dismiss");
    act(() => {
      fireEvent.click(btn);
    });
    expect(screen.queryByTestId("toast")).not.toBeInTheDocument();
    // Advancing past the 3000ms auto-dismiss must be a no-op — the timer
    // should have been cleared on manual close. If it weren't, the stale
    // setItems call would still be a no-op functionally, but vi.getTimerCount
    // gives us a stronger assertion that the handle was actually cleared.
    expect(vi.getTimerCount()).toBe(0);
    act(() => {
      vi.advanceTimersByTime(5000);
    });
    expect(screen.queryByTestId("toast")).not.toBeInTheDocument();
  });

  it("falls back to console.warn after unmount", () => {
    const { unmount } = render(<ToastContainer />);
    unmount();
    const spy = vi.spyOn(console, "warn").mockImplementation(() => {});
    showToast("after unmount", "error");
    expect(spy).toHaveBeenCalledWith("[toast:error] after unmount");
    spy.mockRestore();
  });
});

describe("Toast barrel export", () => {
  it("ToastContainer is exported from the ui barrel", async () => {
    const mod = await import("../src/print/ui");
    expect(typeof mod.ToastContainer).toBe("function");
  });
});
