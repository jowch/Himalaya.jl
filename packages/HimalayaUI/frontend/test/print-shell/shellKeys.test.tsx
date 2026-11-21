import { describe, it, expect, beforeEach, vi, afterEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../../src/state";
import { AppRoutes } from "../../src/print/shell/AppRoutes";
import * as api from "../../src/api";

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("shell keys on the pages outside the workspace", () => {
  beforeEach(() => {
    useAppState.setState({ activeExperimentId: undefined, staleUrlContext: null, navModalOpen: false });
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  });
  afterEach(() => vi.restoreAllMocks());

  it("? opens the shortcuts list for the experiments home and Esc closes it", async () => {
    renderAt("/experiments");
    await screen.findByTestId("app-shell");
    fireEvent.keyDown(document.body, { key: "?" });
    const help = await screen.findByRole("dialog", { name: "Keyboard shortcuts" });
    expect(help).toHaveTextContent("Move between experiments");
    fireEvent.keyDown(help, { key: "Escape" });
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
  });

  it("⌘K opens the finder from inside a text field, where ? and / keep typing", async () => {
    renderAt("/experiments/new");
    await screen.findByTestId("app-shell");
    const field = screen.getByTestId("dirpicker-input").querySelector("input")!;
    expect(document.activeElement).toBe(field);
    fireEvent.keyDown(field, { key: "?" });
    fireEvent.keyDown(field, { key: "/" });
    expect(useAppState.getState().navModalOpen).toBe(false);
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
    fireEvent.keyDown(field, { key: "k", metaKey: true });
    expect(useAppState.getState().navModalOpen).toBe(true);
  });

  it("⌘K does nothing inside an open dialog", async () => {
    renderAt("/experiments");
    await screen.findByTestId("app-shell");
    fireEvent.keyDown(document.body, { key: "?" });
    const help = await screen.findByRole("dialog", { name: "Keyboard shortcuts" });
    fireEvent.keyDown(help, { key: "k", ctrlKey: true });
    expect(useAppState.getState().navModalOpen).toBe(false);
  });

  it.each([
    [{ key: "k", metaKey: true }],
    [{ key: "k", ctrlKey: true }],
    [{ key: "/" }],
  ])("%o opens the finder", async (init) => {
    renderAt("/experiments");
    await screen.findByTestId("app-shell");
    fireEvent.keyDown(document.body, init);
    expect(useAppState.getState().navModalOpen).toBe(true);
  });

  it("ignores ? and / while typing in a field", async () => {
    renderAt("/experiments");
    await screen.findByTestId("app-shell");
    const input = document.createElement("input");
    document.body.appendChild(input);
    input.focus();
    fireEvent.keyDown(input, { key: "?" });
    fireEvent.keyDown(input, { key: "/" });
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
    expect(useAppState.getState().navModalOpen).toBe(false);
    input.remove();
  });
});
