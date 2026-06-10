import { describe, it, expect, beforeEach } from "vitest";
import { renderHook, render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, useLocation } from "react-router-dom";
import type { ReactNode } from "react";
import { useGlobalShortcuts } from "../src/hooks/useGlobalShortcuts";
import { useAppState } from "../src/state";
import type { Sample } from "../src/api";

function wrapperAt(path: string) {
  return function Wrapper({ children }: { children: ReactNode }): JSX.Element {
    return <MemoryRouter initialEntries={[path]}>{children}</MemoryRouter>;
  };
}

function sample(id: number): Sample {
  return { id, experiment_id: 1, name: `S${id}`, display_name: null,
    notes: null, tags: [], q_units: "A-1" } as Sample;
}

// Renders the hook alongside a pathname probe so a test can observe whether the
// ,/. shortcut navigated the URL (focus route) vs. only set store state.
function Harness({ samples }: { samples: Sample[] }): JSX.Element {
  useGlobalShortcuts(samples);
  const loc = useLocation();
  return <div data-testid="pathname">{loc.pathname}</div>;
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

// The ,/. sample step was dead on the focus route: setActiveSample alone is
// reverted by the one-way URL->store sync, so the step must navigate the URL
// there (mirroring the topbar stepper + NavModal's M1 fix).
describe("useGlobalShortcuts — ,/. sample step", () => {
  const SAMPLES = [sample(10), sample(11), sample(12)];

  beforeEach(() => {
    useAppState.setState({ activeSampleId: 11 });
  });

  it("navigates to the next sibling URL on '.' from /sample/:id", () => {
    render(<Harness samples={SAMPLES} />, { wrapper: wrapperAt("/sample/11") });
    fireEvent.keyDown(document.body, { key: "." });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/12");
  });

  it("navigates to the previous sibling URL on ',' from /sample/:id", () => {
    render(<Harness samples={SAMPLES} />, { wrapper: wrapperAt("/sample/11") });
    fireEvent.keyDown(document.body, { key: "," });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/10");
  });

  it("does not wrap past the last sample on the focus route", () => {
    useAppState.setState({ activeSampleId: 12 });
    render(<Harness samples={SAMPLES} />, { wrapper: wrapperAt("/sample/12") });
    fireEvent.keyDown(document.body, { key: "." });
    expect(screen.getByTestId("pathname")).toHaveTextContent("/sample/12");
  });

  it("on a corpus route, '.' sets store state without navigating", () => {
    render(<Harness samples={SAMPLES} />, { wrapper: wrapperAt("/samples") });
    fireEvent.keyDown(document.body, { key: "." });
    // Store advanced, URL unchanged (corpus surfaces are store-driven).
    expect(useAppState.getState().activeSampleId).toBe(12);
    expect(screen.getByTestId("pathname")).toHaveTextContent("/samples");
  });
});

// LO-KEYD: the shared suppressGlobalKeys guard must NOT swallow modifier
// chords — ⌘K is the one binding here that REQUIRES metaKey to pass through.
describe("useGlobalShortcuts — ⌘K passes the suppression guard", () => {
  beforeEach(() => {
    useAppState.setState({ navModalOpen: false, activeExperimentId: undefined });
  });

  it("opens the nav modal on ⌘K from a non-editing target", () => {
    renderHook(() => useGlobalShortcuts(undefined), {
      wrapper: wrapperAt("/samples"),
    });
    fireEvent.keyDown(document.body, { key: "k", metaKey: true });
    expect(useAppState.getState().navModalOpen).toBe(true);
  });
});
