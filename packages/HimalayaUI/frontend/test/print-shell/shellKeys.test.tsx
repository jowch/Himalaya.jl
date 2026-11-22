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
    useAppState.setState({ activeExperimentId: undefined, staleUrlContext: null, paletteOpen: false });
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

  it("⌘K opens the palette from inside a text field, where ? and / type once it has text", async () => {
    renderAt("/experiments/new");
    await screen.findByTestId("app-shell");
    const field = screen.getByTestId("dirpicker-input").querySelector("input")!;
    expect(document.activeElement).toBe(field);
    fireEvent.change(field, { target: { value: "/lab" } });
    fireEvent.keyDown(field, { key: "?" });
    fireEvent.keyDown(field, { key: "/" });
    expect(useAppState.getState().paletteOpen).toBe(false);
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
    fireEvent.keyDown(field, { key: "k", metaKey: true });
    expect(useAppState.getState().paletteOpen).toBe(true);
  });

  it("? from the empty folder field opens the shortcuts list, and types once the field has text", async () => {
    renderAt("/experiments/new");
    await screen.findByTestId("app-shell");
    const field = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.change(field, { target: { value: "/lab" } });
    fireEvent.keyDown(field, { key: "?" });
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
    fireEvent.change(field, { target: { value: "" } });
    fireEvent.keyDown(field, { key: "?" });
    expect(await screen.findByRole("dialog", { name: "Keyboard shortcuts" })).toBeInTheDocument();
  });

  it("⌘? and Ctrl+? in the empty folder field do not open the shortcuts", async () => {
    renderAt("/experiments/new");
    await screen.findByTestId("app-shell");
    const field = screen.getByTestId("dirpicker-input").querySelector("input")!;
    fireEvent.keyDown(field, { key: "?", metaKey: true });
    fireEvent.keyDown(field, { key: "?", ctrlKey: true });
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
  });

  it("? types into an empty field inside a dialog instead of opening the shortcuts over it", async () => {
    renderAt("/experiments");
    await screen.findByTestId("app-shell");
    const dialog = document.createElement("div");
    dialog.setAttribute("role", "dialog");
    dialog.setAttribute("aria-modal", "true");
    const input = document.createElement("input");
    dialog.appendChild(input);
    document.body.appendChild(dialog);
    input.focus();
    try {
      fireEvent.keyDown(input, { key: "?" });
      expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
    } finally {
      dialog.remove();
    }
  });

  it("? types into the empty palette field on the experiments home instead of opening the shortcuts", async () => {
    renderAt("/experiments");
    await screen.findByTestId("app-shell");
    fireEvent.keyDown(document.body, { key: "k", metaKey: true });
    const palette = await screen.findByRole("dialog");
    const input = palette.querySelector("input")!;
    expect(input).toHaveValue("");
    fireEvent.keyDown(input, { key: "?" });
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
  });

  it("⌘K does nothing inside an open dialog", async () => {
    renderAt("/experiments");
    await screen.findByTestId("app-shell");
    fireEvent.keyDown(document.body, { key: "?" });
    const help = await screen.findByRole("dialog", { name: "Keyboard shortcuts" });
    fireEvent.keyDown(help, { key: "k", ctrlKey: true });
    expect(useAppState.getState().paletteOpen).toBe(false);
  });

  it.each([
    [{ key: "k", metaKey: true }],
    [{ key: "k", ctrlKey: true }],
    [{ key: "/" }],
  ])("%o opens the palette", async (init) => {
    renderAt("/experiments");
    await screen.findByTestId("app-shell");
    fireEvent.keyDown(document.body, init);
    expect(useAppState.getState().paletteOpen).toBe(true);
  });

  it("ignores ? and / while typing in a field", async () => {
    renderAt("/experiments");
    await screen.findByTestId("app-shell");
    const input = document.createElement("input");
    input.value = "a";
    document.body.appendChild(input);
    input.focus();
    fireEvent.keyDown(input, { key: "?" });
    fireEvent.keyDown(input, { key: "/" });
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
    expect(useAppState.getState().paletteOpen).toBe(false);
    input.remove();
  });
});
