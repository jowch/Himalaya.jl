// test/ExperimentShell.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentShell } from "../src/print/shell/ExperimentShell";
import * as api from "../src/api";

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
};

describe("ExperimentShell (Phase E1)", () => {
  beforeEach(() => { vi.spyOn(api, "getExperiment").mockResolvedValue(EXP); });
  afterEach(() => vi.restoreAllMocks());

  it("renders own chrome (top nav) + the experiment header name", async () => {
    renderAt("/experiments/7/corpus");
    expect(screen.getByTestId("experiment-top-nav")).toBeInTheDocument();
    // The header name is an Input variant='title' (edit-in-place). Its wrapper
    // div carries data-testid; the inner <input> carries the value attribute.
    const nameWrapper = await screen.findByTestId("experiment-header-name");
    expect(nameWrapper.querySelector("input")?.value).toBe("SSRL · 1p7m");
  });

  it("renders the Corpus | Configuration tab bar with the active tab", async () => {
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    expect(screen.getByTestId("exp-tab-corpus")).toHaveAttribute("aria-current", "page");
    expect(screen.getByTestId("exp-tab-config")).not.toHaveAttribute("aria-current");
  });

  it("renders the child route via Outlet", async () => {
    renderAt("/experiments/7/config");
    expect(await screen.findByText("CONFIG BODY")).toBeInTheDocument();
    expect(screen.getByTestId("exp-tab-config")).toHaveAttribute("aria-current", "page");
  });

  it("does NOT render the corpus topbar (it lives outside CorpusShell)", async () => {
    renderAt("/experiments/7/corpus");
    await screen.findByTestId("experiment-header-name");
    expect(screen.queryByTestId("corpus-topbar")).toBeNull();
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
});
