// test/ConfigurationPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ConfigurationPage } from "../src/print/pages/ConfigurationPage";
import * as api from "../src/api";

const EXP: api.Experiment = {
  id: 7, name: "SSRL", description: null, path: "/d", data_dir: "/d/data", analysis_dir: "/d/an",
  manifest_path: null, created_at: "2026-04-12", q_units: "A^-1", beam_center_x: 421,
  beam_center_y: 836, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8095,
  energy_kev_source: "prp", flight_path_m_source: "setup", beam_center_x_source: "setup",
  beam_center_y_source: "setup", pixel_size_um_source: "prp", q_units_source: "prp",
  last_scanned_at: null, scan_signature: null, ingest_status: "complete",
  image_pattern: null, metadata_pattern: null, integration_pattern: null,
};

function renderPage() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7/config"]}>
        <Routes>
          <Route path="/experiments/:id/config" element={<ConfigurationPage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("ConfigurationPage (Phase E1 shell)", () => {
  beforeEach(() => { vi.spyOn(api, "getExperiment").mockResolvedValue(EXP); });
  afterEach(() => vi.restoreAllMocks());

  it("renders the Geometry and Sources regions", async () => {
    renderPage();
    expect(await screen.findByTestId("config-geometry-region")).toBeInTheDocument();
    expect(screen.getByTestId("config-sources-region")).toBeInTheDocument();
  });
});
