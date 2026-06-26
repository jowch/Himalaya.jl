// test/KbdOverlay.test.tsx
//
// T2.5: KbdOverlay — the `?` help overlay.
//
// Pins:
//   - Renders the KbdLegend (keyboard shortcut list) when open.
//   - Esc closes it (via ModalShell's built-in closeOnEsc).
//   - ? key opens it (wired through the shell-action registry + useKeyboardLayer).
//
// The overlay is controlled by useAppState.helpOverlayOpen / openHelpOverlay /
// closeHelpOverlay. We drive it through the real Zustand store so the
// shellActions ↔ KbdOverlay integration is tested end-to-end.
import { describe, it, expect, vi, afterEach, beforeEach } from "vitest";
import { render, screen, fireEvent, cleanup } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { useAppState } from "../src/state";
import { KbdOverlay } from "../src/print/shell/KbdOverlay";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../src/print/interaction/registry";
import { core } from "../src/print/interaction/core";

// useFocusTrap uses MutationObserver which JSDOM stubs but doesn't implement
// fully. Stub it to avoid noise.
vi.mock("../src/hooks/useFocusTrap", () => ({
  useFocusTrap: () => {},
}));

/** Mirrors AppShell's role: mounts the keyboard layer and seeds the help shell action. */
function HarnessWithOverlay(): JSX.Element {
  useKeyboardLayer();
  return <KbdOverlay />;
}

afterEach(() => {
  useInteraction.getState().setShellActions([]);
  cleanup();
});

function renderHarness() {
  return render(
    <MemoryRouter>
      <HarnessWithOverlay />
    </MemoryRouter>,
  );
}

beforeEach(() => {
  // Reset the overlay to closed before each test.
  useAppState.setState({ helpOverlayOpen: false });
});

describe("KbdOverlay — T2.5 ? help overlay", () => {
  it("is not rendered when helpOverlayOpen is false", () => {
    renderHarness();
    expect(screen.queryByTestId("kbd-overlay")).toBeNull();
  });

  it("renders when helpOverlayOpen is true (contains the legend)", () => {
    useAppState.setState({ helpOverlayOpen: true });
    renderHarness();
    expect(screen.getByTestId("kbd-overlay")).toBeInTheDocument();
    // The legend inside must be present (registry-driven shortcut list).
    expect(screen.getByTestId("kbd-overlay-legend")).toBeInTheDocument();
  });

  it("Esc closes the overlay", () => {
    useAppState.setState({ helpOverlayOpen: true });
    renderHarness();
    expect(screen.getByTestId("kbd-overlay")).toBeInTheDocument();
    fireEvent.keyDown(document, { key: "Escape" });
    expect(screen.queryByTestId("kbd-overlay")).toBeNull();
    expect(useAppState.getState().helpOverlayOpen).toBe(false);
  });

  it("? key opens the overlay via the shell-action registry + useKeyboardLayer", () => {
    // Seed the help shell action — this is what AppShell's useEffect does on mount.
    useInteraction.getState().setShellActions([
      core("help", { run: () => useAppState.getState().openHelpOverlay() }),
    ]);
    renderHarness();
    expect(screen.queryByTestId("kbd-overlay")).toBeNull();
    fireEvent.keyDown(window, { key: "?", shiftKey: true }); // real `?` is Shift+/
    expect(useAppState.getState().helpOverlayOpen).toBe(true);
    // Re-render to reflect state change (not using act here as state is sync).
    // The overlay should now render.
    expect(screen.getByTestId("kbd-overlay")).toBeInTheDocument();
  });
});
