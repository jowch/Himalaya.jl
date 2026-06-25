// test/KbdOverlay.test.tsx
//
// T2.5: KbdOverlay — the `?` help overlay.
//
// Pins:
//   - Renders the KbdLegend (keyboard shortcut list) when open.
//   - Esc closes it (via ModalShell's built-in closeOnEsc).
//   - ? key opens it (wired in useGlobalShortcuts).
//
// The overlay is controlled by useAppState.helpOverlayOpen / openHelpOverlay /
// closeHelpOverlay. We drive it through the real Zustand store so the
// useGlobalShortcuts ↔ KbdOverlay integration is tested end-to-end.
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { useAppState } from "../src/state";
import { KbdOverlay } from "../src/print/shell/KbdOverlay";
import { useGlobalShortcuts } from "../src/hooks/useGlobalShortcuts";

// useFocusTrap uses MutationObserver which JSDOM stubs but doesn't implement
// fully. Stub it to avoid noise.
vi.mock("../src/hooks/useFocusTrap", () => ({
  useFocusTrap: () => {},
}));

function HarnessWithOverlay(): JSX.Element {
  useGlobalShortcuts();
  return <KbdOverlay />;
}

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

  it("? key opens the overlay via useGlobalShortcuts", () => {
    renderHarness();
    expect(screen.queryByTestId("kbd-overlay")).toBeNull();
    fireEvent.keyDown(window, { key: "?" });
    expect(useAppState.getState().helpOverlayOpen).toBe(true);
    // Re-render to reflect state change (not using act here as state is sync).
    // The overlay should now render.
    expect(screen.getByTestId("kbd-overlay")).toBeInTheDocument();
  });
});
