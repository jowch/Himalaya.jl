/**
 * Task 0.7: proves that AppShell mounts both InteractionDock and useKeyboardLayer.
 *
 * Strategy (not the false-green from the brief's Step 1):
 *   - Render AppRoutes at /series (lightweight; mocking listSeries suffices).
 *   - After the shell mounts, seed the interaction store.
 *   - dock-primary appearing proves InteractionDock is rendered inside the shell.
 *   - Enter firing the spy proves useKeyboardLayer's window listener is live.
 *
 * Scope note: Enter is NOT a bare key (`"Enter".length === 5 !== 1`), so
 * isBareKey() returns false and the page-scope guard never applies. Firing
 * keyDown on document.body is therefore unambiguous — no data-interaction-scope
 * element needed.
 */
import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../../src/print/shell/AppRoutes";
import { useInteraction } from "../../src/print/interaction/registry";
import { core } from "../../src/print/interaction/core";
import * as api from "../../src/api";

function renderAtSeries() {
  const qc = new QueryClient({
    defaultOptions: { queries: { retry: false } },
  });
  vi.spyOn(api, "listSeries").mockResolvedValue([]);
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/series"]}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

afterEach(() => {
  useInteraction.getState().clearPage();
  vi.restoreAllMocks();
});

describe("AppShell mounts the interaction layer (Task 0.7)", () => {
  it("InteractionDock is mounted in the shell (dock-primary appears after setPage)", async () => {
    renderAtSeries();
    // Wait for the shell to be in the DOM.
    await screen.findByTestId("app-shell");

    const spy = vi.fn();
    // Seed the store after mount — InteractionDock must be subscribed inside the
    // shell and re-render when the store changes. Wrap in act so the Zustand
    // subscriber's re-render is batched correctly.
    act(() => {
      useInteraction.getState().setPage(null, [core("openFocus", { run: spy, dock: "primary" })]);
    });

    // If InteractionDock is mounted, this resolves; if not, it times out → RED.
    expect(await screen.findByTestId("dock-primary")).toBeInTheDocument();
  });

  it("useKeyboardLayer window listener is live (Enter fires registered action)", async () => {
    renderAtSeries();
    await screen.findByTestId("app-shell");

    const spy = vi.fn();
    act(() => {
      useInteraction.getState().setPage(null, [core("openFocus", { run: spy, dock: "primary" })]);
    });

    // Wait for the dock to confirm the store is seeded and the component settled.
    await screen.findByTestId("dock-primary");

    // Enter is NOT bare (isBareKey checks e.key.length === 1; "Enter".length === 5),
    // so the page-scope guard is skipped and the action fires from any target.
    fireEvent.keyDown(document.body, { key: "Enter" });

    // If useKeyboardLayer is mounted, the spy fires; if not → RED.
    expect(spy).toHaveBeenCalledTimes(1);
  });
});
