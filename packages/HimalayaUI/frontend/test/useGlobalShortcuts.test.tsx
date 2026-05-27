import { describe, it, expect, beforeEach } from "vitest";
import { renderHook, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import type { ReactNode } from "react";
import { useGlobalShortcuts } from "../src/hooks/useGlobalShortcuts";
import { useAppState } from "../src/state";

function wrapperAt(path: string) {
  return function Wrapper({ children }: { children: ReactNode }): JSX.Element {
    return <MemoryRouter initialEntries={[path]}>{children}</MemoryRouter>;
  };
}

// Regression guard (#161): the loupe (a corpus surface) binds the arrow keys
// to exposure-flipping. `useGlobalShortcuts` is hoisted app-wide and also
// binds the arrows to a legacy page-tab step (`setActivePage`). Both listeners
// sit on `window`, so the loupe cannot defend itself with stopPropagation —
// the page-tab step must be gated out of corpus routes here instead.
describe("useGlobalShortcuts — arrow-key page-tab step gating", () => {
  beforeEach(() => {
    useAppState.setState({ activePage: "compare" });
  });

  it("does not step activePage on a corpus route (the loupe owns the arrows there)", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/samples/loupe/7"),
    });
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    expect(useAppState.getState().activePage).toBe("compare");
  });

  it("arrow keys are a no-op on a legacy route now that TAB_ORDER is single-tab (#181)", () => {
    // I4.4 (#181): Index retired → TAB_ORDER = ["compare"]. With one tab,
    // ArrowRight/ArrowLeft clamp to the current tab — nothing to step to.
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/compare/all"),
    });
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    expect(useAppState.getState().activePage).toBe("compare");
    fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    expect(useAppState.getState().activePage).toBe("compare");
  });
});
