import { describe, it, expect, beforeEach, vi, afterEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { StrictMode } from "react";
import { MemoryRouter, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../../src/state";
import { AppRoutes } from "../../src/shell/AppRoutes";
import * as api from "../../src/api";
import type { CorpusSample, Experiment } from "../../src/api";

let location = "";
function Probe(): null {
  location = useLocation().pathname;
  return null;
}

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <StrictMode>
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={[path]}>
          <AppRoutes />
          <Probe />
        </MemoryRouter>
      </QueryClientProvider>
    </StrictMode>,
  );
}

const experiment = (id: number, name: string): Experiment => ({
  id, name, description: null, path: `/data/${name}`, created_at: "2026-11-01T00:00:00Z",
  last_scanned_at: null, ingest_status: "idle",
});

const sample = (id: number, experiment_id: number, name: string): CorpusSample => ({
  id, experiment_id, name, notes: null, q_units: "A^-1",
});

async function openPalette(path: string): Promise<HTMLElement> {
  renderAt(path);
  await screen.findByTestId("app-shell");
  fireEvent.keyDown(document.body, { key: "k", metaKey: true });
  return screen.findByTestId("workspace-palette");
}

function input(): HTMLElement {
  return screen.getByRole("combobox");
}

describe("the palette on the pages outside the workspace", () => {
  beforeEach(() => {
    localStorage.clear();
    location = "";
    useAppState.setState({ username: "alice", paletteOpen: false });
    vi.spyOn(api, "listExperiments").mockResolvedValue([experiment(1, "Lipids"), experiment(2, "Peptides")]);
    vi.spyOn(api, "listCorpusSamples").mockResolvedValue([sample(11, 1, "MO_PEG5"), sample(12, 2, "AMP_weak")]);
  });
  afterEach(() => vi.restoreAllMocks());

  it("lists the experiments and the app commands on the home page", async () => {
    const palette = await openPalette("/experiments");
    await within(palette).findByRole("option", { name: /Open Lipids/ });
    expect(within(palette).getByRole("option", { name: /Open Peptides/ })).toBeInTheDocument();
    expect(within(palette).getByRole("option", { name: /New experiment/ })).toBeInTheDocument();
    expect(within(palette).queryByRole("option", { name: /All experiments/ })).toBeNull();
    expect(within(palette).getByRole("option", { name: /Theme/ })).toBeInTheDocument();
    expect(within(palette).getByRole("option", { name: /Switch user/ })).toBeInTheDocument();
    expect(within(palette).getByRole("option", { name: /Keyboard shortcuts/ })).toBeInTheDocument();
  });

  it("offers All experiments, not New experiment, in the New-experiment funnel", async () => {
    const palette = await openPalette("/experiments/new");
    await within(palette).findByRole("option", { name: /Open Lipids/ });
    expect(within(palette).queryByRole("option", { name: /New experiment/ })).toBeNull();
    fireEvent.click(within(palette).getByRole("option", { name: /All experiments/ }));
    expect(location).toBe("/experiments");
    expect(screen.queryByTestId("workspace-palette")).toBeNull();
    expect(useAppState.getState().paletteOpen).toBe(false);
  });

  it("opens an experiment", async () => {
    const palette = await openPalette("/experiments");
    fireEvent.click(await within(palette).findByRole("option", { name: /Open Peptides/ }));
    expect(location).toBe("/experiments/2");
  });

  it("finds a sample in any experiment by name and opens it in the workspace", async () => {
    const palette = await openPalette("/experiments");
    await within(palette).findByRole("option", { name: /Open Lipids/ });
    expect(within(palette).queryByRole("option", { name: /Go to MO_PEG5/ })).toBeNull();
    expect(api.listCorpusSamples).not.toHaveBeenCalled();
    fireEvent.change(input(), { target: { value: "peg5" } });
    const hit = await within(palette).findByRole("option", { name: /Go to MO_PEG5/ });
    expect(hit).toHaveTextContent("Lipids");
    fireEvent.click(hit);
    expect(location).toBe("/experiments/1/workspace/11");
  });

  it("sets the theme", async () => {
    const palette = await openPalette("/experiments");
    fireEvent.click(await within(palette).findByRole("option", { name: /Theme/ }));
    fireEvent.click(within(screen.getByTestId("workspace-palette")).getByRole("option", { name: /Dark/ }));
    expect(localStorage.getItem("himalaya-ui:theme")).toBe("dark");
  });

  it("switches user", async () => {
    const palette = await openPalette("/experiments");
    fireEvent.click(await within(palette).findByRole("option", { name: /Switch user/ }));
    expect(useAppState.getState().username).toBeUndefined();
  });

  it("opens the shortcuts list", async () => {
    const palette = await openPalette("/experiments");
    fireEvent.click(await within(palette).findByRole("option", { name: /Keyboard shortcuts/ }));
    expect(await screen.findByRole("dialog", { name: "Keyboard shortcuts" })).toBeInTheDocument();
  });

  it("⌘K again and Esc close it", async () => {
    await openPalette("/experiments");
    fireEvent.keyDown(input(), { key: "k", metaKey: true });
    expect(screen.queryByTestId("workspace-palette")).toBeNull();
    fireEvent.keyDown(document.body, { key: "k", metaKey: true });
    fireEvent.keyDown(await screen.findByRole("combobox"), { key: "Escape" });
    expect(screen.queryByTestId("workspace-palette")).toBeNull();
    expect(useAppState.getState().paletteOpen).toBe(false);
  });

});
