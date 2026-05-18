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
    useAppState.setState({ activePage: "index" });
  });

  it("does not step activePage on a corpus route (the loupe owns the arrows there)", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/samples/loupe/7"),
    });
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    expect(useAppState.getState().activePage).toBe("index");
  });

  it("still steps activePage on a legacy route", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/index"),
    });
    // TAB_ORDER is ["inspect","index","compare"] → ArrowRight: index → compare.
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    expect(useAppState.getState().activePage).toBe("compare");
  });
});
