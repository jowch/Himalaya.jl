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

// I5.1 (#182): the dual-nav `activePage` model + its ArrowLeft/Right page-tab
// step are deleted. The arrow keys are no longer bound by useGlobalShortcuts,
// so corpus surfaces (e.g. the loupe at /samples/loupe/:id) own them outright.
// These tests guard against a regression that re-binds them globally.
describe("useGlobalShortcuts — arrow keys are unbound (no page-tab step)", () => {
  beforeEach(() => {
    useAppState.setState({ activeSampleId: undefined });
  });

  it("does not mutate any nav state on a corpus route (loupe owns the arrows)", () => {
    const before = useAppState.getState();
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/samples/loupe/7"),
    });
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    const after = useAppState.getState();
    // The hook must not touch the active ids on an arrow press.
    expect(after.activeSampleId).toBe(before.activeSampleId);
    expect(after.activeExperimentId).toBe(before.activeExperimentId);
  });

  it("is a no-op on a non-corpus route too (no global arrow binding remains)", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/totally/unknown"),
    });
    // Should neither throw nor mutate state.
    expect(() => {
      fireEvent.keyDown(document.body, { key: "ArrowRight" });
      fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    }).not.toThrow();
  });

  // R0a (#221): the `T` theme-toggle shortcut is retired with the dark theme.
  // Pressing T must no longer mutate any store state (no `theme` field exists).
  it("does not bind T (theme toggle retired with the dark theme)", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/samples"),
    });
    const before = useAppState.getState();
    expect("theme" in (before as Record<string, unknown>)).toBe(false);
    expect(() => {
      fireEvent.keyDown(document.body, { key: "T" });
      fireEvent.keyDown(document.body, { key: "t" });
    }).not.toThrow();
    // No store mutation from T.
    expect(useAppState.getState().activeSampleId).toBe(before.activeSampleId);
  });
});
