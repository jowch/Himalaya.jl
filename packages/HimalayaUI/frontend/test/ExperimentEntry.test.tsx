import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentEntry } from "../src/print/shell/ExperimentEntry";
import { useAppState } from "../src/state";
import * as api from "../src/api";
import type { Experiment } from "../src/api";

const EXPERIMENT: Experiment = {
  id: 7, name: "TestExp", description: null, path: "/data/exp1",
  data_dir: "/data/exp1/raw", analysis_dir: "/data/exp1/analysis",
  manifest_path: null, created_at: "2026-01-01T00:00:00Z",
  q_units: null, beam_center_x: null, beam_center_y: null,
  pixel_size_um: null, energy_kev: null, flight_path_m: null,
  energy_kev_source: "computed", flight_path_m_source: "computed",
  beam_center_x_source: "computed", beam_center_y_source: "computed",
  pixel_size_um_source: "computed", q_units_source: "computed",
  last_scanned_at: "2026-01-01T10:00:00Z", ingest_status: "complete",
  image_pattern: "*.tif", metadata_pattern: null, integration_pattern: null,
};

function Where(): JSX.Element {
  return <div data-testid="where">{useLocation().pathname}</div>;
}

function renderAt(experiment: Experiment, inFlight: { processed: number; total: number; status: "scanning" | "analyzing" } | null = null) {
  useAppState.setState({ ingestInFlight: inFlight ? { 7: inFlight } : null });
  vi.spyOn(api, "getExperiment").mockResolvedValue(experiment);
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7"]}>
        <Routes>
          <Route path="/experiments/:id" element={<ExperimentEntry />} />
          <Route path="*" element={<Where />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

afterEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ ingestInFlight: null });
});

describe("ExperimentEntry", () => {
  it("opens a scanned experiment in the workspace", async () => {
    renderAt(EXPERIMENT);
    expect(await screen.findByTestId("where")).toHaveTextContent("/experiments/7/workspace");
  });

  it("opens the workspace during a rescan", async () => {
    renderAt({ ...EXPERIMENT, ingest_status: "scanning" }, { processed: 1, total: 4, status: "analyzing" });
    expect(await screen.findByTestId("where")).toHaveTextContent("/experiments/7/workspace");
  });

  it("opens the workspace during the first ingest", async () => {
    renderAt({ ...EXPERIMENT, last_scanned_at: null, ingest_status: "scanning" }, { processed: 1, total: 4, status: "scanning" });
    expect(await screen.findByTestId("where")).toHaveTextContent("/experiments/7/workspace");
  });

  it("opens the workspace after a failed ingest", async () => {
    renderAt({ ...EXPERIMENT, ingest_status: "failed" });
    expect(await screen.findByTestId("where")).toHaveTextContent("/experiments/7/workspace");
  });
});
