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
    useAppState.setState({ activePage: "none" });
  });

  it("does not step activePage on a corpus route (the loupe owns the arrows there)", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/samples/loupe/7"),
    });
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    expect(useAppState.getState().activePage).toBe("none");
  });

  it("arrow keys are a no-op on a legacy route now that every surface is retired (#177)", () => {
    // I3.6 (#177): Compare retired (Index #181 / Inspect #163 before it). The
    // only PageId left is the inert "none" sentinel → TAB_ORDER = ["none"], so
    // ArrowRight/ArrowLeft clamp to it and never step. (The page-tab step + the
    // whole dual-nav model retire in I5.1.) A stale `/compare/all` URL also
    // redirects to /series at the router, so it never mounts AppShell anyway.
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/compare/all"),
    });
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    expect(useAppState.getState().activePage).toBe("none");
    fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    expect(useAppState.getState().activePage).toBe("none");
  });
});
