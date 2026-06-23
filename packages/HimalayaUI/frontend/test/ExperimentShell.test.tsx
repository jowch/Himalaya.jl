// test/ExperimentShell.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentShell } from "../src/print/shell/ExperimentShell";
import * as api from "../src/api";
import { useAppState } from "../src/state";

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}>
        <Routes>
          <Route path="/experiments/:id" element={<ExperimentShell />}>
            <Route path="corpus" element={<div>CORPUS BODY</div>} />
            <Route path="config" element={<div>CONFIG BODY</div>} />
          </Route>
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

const EXP: api.Experiment = {
  id: 7, name: "SSRL · 1p7m", description: null, path: "/d", data_dir: "/d/data", analysis_dir: "/d/an",
  manifest_path: null, created_at: "2026-04-12", q_units: "A^-1",
  beam_center_x: 1, beam_center_y: 1, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8,
  energy_kev_source: "prp", flight_path_m_source: "setup", beam_center_x_source: "setup",
  beam_center_y_source: "setup", pixel_size_um_source: "prp", q_units_source: "prp",
  last_scanned_at: "2026-04-12T10:00:00", scan_signature: "sig", ingest_status: "complete",
  image_pattern: null, metadata_pattern: null, integration_pattern: null,
  stats: { loads: 4, samples: 12, exposures: 48, sessions: 2, span_hours: 6.5, started_at: "2026-04-12T10:00:00" },
};

describe("ExperimentShell (Phase E1)", () => {
  beforeEach(() => { vi.spyOn(api, "getExperiment").mockResolvedValue(EXP); });
  afterEach(() => vi.restoreAllMocks());

  it("renders the experiment header name (T3.2: TopNav lives in outer AppShell)", async () => {
    // T3.2: ExperimentShell is pure page content; TopNav comes from the outer
    // AppShell (not rendered in this isolated unit test).
    renderAt("/experiments/7/corpus");
    // The header name is an Input variant='title' (edit-in-place). Its wrapper
    // div carries data-testid; the inner <input> carries the value attribute.
    const nameWrapper = await screen.findByTestId("experiment-header-name");
    expect(nameWrapper.querySelector("input")?.value).toBe("SSRL · 1p7m");
  });

  it("retires the tab bar; Configuration is reached via the ⚙ gear (M3)", async () => {
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    expect(screen.queryByTestId("experiment-tab-bar")).toBeNull();
    expect(screen.getByTestId("experiment-config-gear")).toBeInTheDocument();
  });

  it("renders the child route via Outlet", async () => {
    renderAt("/experiments/7/config");
    expect(await screen.findByText("CONFIG BODY")).toBeInTheDocument();
  });

  it("no corpus-topbar in isolation (T3.2: TopNav lives in outer AppShell)", async () => {
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    // T3.2: corpus-topbar is gone; ExperimentShell is pure page content.
    // TopNav is provided by the outer AppShell (not in this isolated unit test).
    expect(screen.queryByTestId("corpus-topbar")).toBeNull();
    expect(screen.queryByTestId("topnav")).toBeNull();
  });

  it("StatBar shows real counts from stats when available (not placeholder dashes)", async () => {
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    // The stats fixture has loads=4, samples=12, exposures=48, sessions=2.
    // None of these are "—" — the StatBar must surface the real numbers.
    expect(screen.getByText("4")).toBeInTheDocument();
    expect(screen.getByText("12")).toBeInTheDocument();
    expect(screen.getByText("48")).toBeInTheDocument();
    expect(screen.getByText("2")).toBeInTheDocument();
    // Placeholder dashes must not appear in the stat bar.
    const statBar = screen.getByRole("group", { name: /experiment stats/i });
    expect(statBar).not.toHaveTextContent("—");
  });

  it("commits name on blur: calls updateExperiment with trimmed value", async () => {
    // Spy on updateExperiment before rendering so we capture all calls.
    const updateSpy = vi.spyOn(api, "updateExperiment").mockResolvedValue(EXP);

    renderAt("/experiments/7/corpus");

    // Wait for the Input to appear (it is suppressed during loading).
    const wrapper = await screen.findByTestId("experiment-header-name");
    const input = wrapper.querySelector("input") as HTMLInputElement;
    expect(input).not.toBeNull();

    // Simulate the user editing the field: fire a change event with surrounding
    // whitespace (proves the component trims before calling the API), then fire
    // blur to trigger commitName → updateExperiment.
    // fireEvent.change directly sets input.value and fires onChange, which wires
    // through onValueChange → setLocalEdit(true) + setPendingDraft(v).
    fireEvent.change(input, { target: { value: "  BL 7.3.3 SMB  " } });
    fireEvent.blur(input); // triggers onBlur → commitName

    expect(updateSpy).toHaveBeenCalledTimes(1);
    expect(updateSpy).toHaveBeenCalledWith(
      7,                         // expId parsed from the URL param
      { name: "BL 7.3.3 SMB" }, // trimmed
      {},                        // authOpts(undefined, undefined) with no username set
    );
  });

  // --- Rescan button ---

  it("renders a Rescan button when not processing", async () => {
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    expect(screen.getByTestId("experiment-rescan-button")).toBeInTheDocument();
  });

  it("Rescan button calls api.triggerScan on click", async () => {
    const scanSpy = vi.spyOn(api, "triggerScan").mockResolvedValue(EXP);
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    fireEvent.click(screen.getByTestId("experiment-rescan-button"));
    // triggerScan(id, authOpts, force) — the manual Rescan button is a routine
    // (additive) rescan, so force=false (a pattern edit passes force=true).
    await waitFor(() => expect(scanSpy).toHaveBeenCalledWith(7, expect.anything(), false));
  });

  it("Rescan button is disabled when isProcessing (status=scanning)", async () => {
    // A real in-flight INITIAL scan has the persisted row at "scanning" too (the
    // create route sets it). effectiveIngestStatus reads a scanning overlay as
    // stale unless the persisted row agrees — so mirror reality here, otherwise
    // this would exercise the 8c stale-overlay path (overlay scanning + row
    // complete → resolves to complete), not an active scan.
    vi.spyOn(api, "getExperiment").mockResolvedValue({ ...EXP, ingest_status: "scanning" });
    useAppState.setState({
      ingestInFlight: { 7: { status: "scanning", processed: 5, total: 20 } },
    });
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    expect(screen.getByTestId("experiment-rescan-button")).toBeDisabled();
    // Reset state for subsequent tests.
    useAppState.setState({ ingestInFlight: null });
  });

  it("Rescan button is disabled when isProcessing (status=analyzing)", async () => {
    useAppState.setState({
      ingestInFlight: { 7: { status: "analyzing", processed: 20, total: 20 } },
    });
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    expect(screen.getByTestId("experiment-rescan-button")).toBeDisabled();
    useAppState.setState({ ingestInFlight: null });
  });
});
