/**
 * Task 6.2: shellActions slot — scope-exemption regression.
 *
 * Strategy: mount a minimal TestShell (just useKeyboardLayer, like keyboardLayer.test.tsx),
 * seed shellActions via the registry directly (same actions AppShell sets in its useEffect),
 * and fire key events from document.body to prove the bare-key scope guard is NOT applied.
 *
 * The critical case is the global-scope regression: `/` is a bare key (length === 1),
 * so page actions would require [data-interaction-scope] to fire.  Shell actions bypass
 * the guard and fire from any target — including document.body which has no ancestor scope.
 */
import { describe, it, expect, vi, afterEach } from "vitest";
import { render, fireEvent, cleanup } from "@testing-library/react";
import { useKeyboardLayer } from "../../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../../src/print/interaction/registry";
import { core } from "../../src/print/interaction/core";
import { useAppState } from "../../src/state";

/** Mirrors AppShell's minimal role: mounts the keyboard layer once. */
function TestShell(): JSX.Element {
  useKeyboardLayer();
  return <div />;
}

afterEach(() => {
  useInteraction.getState().clearPage();
  useInteraction.getState().setShellActions([]);
  cleanup();
  vi.restoreAllMocks();
});

/** Seed the registry with the same shell actions AppShell's useEffect would set. */
function seedShellActions() {
  useInteraction.getState().setShellActions([
    core("find", {
      run: () => {
        const s = useAppState.getState();
        s.openNavModal(s.activeExperimentId === undefined ? "experiment" : "sample");
      },
    }),
    core("help", { run: () => useAppState.getState().openHelpOverlay() }),
  ]);
}

describe("shellActions — find/help bypass the bare-key scope guard", () => {
  it("GLOBAL SCOPE REGRESSION: `/` from document.body (no scope ancestor) calls openNavModal", () => {
    // This is the defining regression: a bare key (`/`, length === 1) that would be
    // blocked by `!inPageScope(e.target)` for page actions — but shell actions skip
    // that guard entirely.  A page action registered with key `/` would NOT fire here.
    const spy = vi.spyOn(useAppState.getState(), "openNavModal");
    seedShellActions();
    render(<TestShell />);
    fireEvent.keyDown(document.body, { key: "/" });
    expect(spy).toHaveBeenCalledTimes(1);
  });

  it("`/` passes the right step when no experiment is selected (activeExperimentId is undefined)", () => {
    // Default store state: activeExperimentId === undefined → step "experiment".
    const spy = vi.spyOn(useAppState.getState(), "openNavModal");
    seedShellActions();
    render(<TestShell />);
    fireEvent.keyDown(document.body, { key: "/" });
    expect(spy).toHaveBeenCalledWith("experiment");
  });

  it("⌘K (metaKey+k) calls openNavModal from document.body", () => {
    const spy = vi.spyOn(useAppState.getState(), "openNavModal");
    seedShellActions();
    render(<TestShell />);
    fireEvent.keyDown(document.body, { key: "k", metaKey: true });
    expect(spy).toHaveBeenCalledTimes(1);
  });

  it("`?` from document.body calls openHelpOverlay", () => {
    const spy = vi.spyOn(useAppState.getState(), "openHelpOverlay");
    seedShellActions();
    render(<TestShell />);
    fireEvent.keyDown(document.body, { key: "?", shiftKey: true }); // real `?` is Shift+/
    expect(spy).toHaveBeenCalledTimes(1);
  });

  it("`/` typed inside an <input> does NOT fire (isTyping guard)", () => {
    // The isTyping guard fires BEFORE the shellActions loop, so even scope-exempt
    // shell actions do not fire when the user is typing.
    const spy = vi.spyOn(useAppState.getState(), "openNavModal");
    seedShellActions();
    render(<><TestShell /><input data-testid="field" /></>);
    const input = document.querySelector<HTMLInputElement>('[data-testid="field"]')!;
    input.focus();
    fireEvent.keyDown(input, { key: "/" });
    expect(spy).not.toHaveBeenCalled();
  });

  it("page action with key `/` does NOT fire when focus is on an out-of-scope element (scope guard is page-only)", () => {
    // Proves the scope-exemption is selective: a PAGE action registered with a bare key
    // is blocked when focus is on an HTMLElement with no [data-interaction-scope] ancestor
    // and is not document.body.  (inPageScope returns false for such elements.)
    // Shell actions bypass this guard — so the same bare key would fire as a shellAction.
    const pageRun = vi.fn();
    useInteraction.getState().setPage(null, [
      { id: "page-slash", label: "Test", keys: ["/"], group: "Act", run: pageRun },
    ]);
    // No shellActions seeded — only the page action is registered.
    render(
      <>
        <TestShell />
        {/* A focusable div with NO data-interaction-scope ancestor. inPageScope() → false. */}
        <div data-testid="out-of-scope" tabIndex={0} />
      </>,
    );
    const el = document.querySelector<HTMLElement>('[data-testid="out-of-scope"]')!;
    el.focus();
    fireEvent.keyDown(el, { key: "/" });
    expect(pageRun).not.toHaveBeenCalled();
  });
});
