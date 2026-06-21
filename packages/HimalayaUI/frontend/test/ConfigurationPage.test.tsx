// test/ConfigurationPage.test.tsx
import { describe, it, test, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ConfigurationPage } from "../src/print/pages/ConfigurationPage";
import * as api from "../src/api";
import { useDraftExperiment } from "../src/lib/draftExperiment";

// Module-level navigate mock — canonical pattern (mirrors NewExperimentPage.test.tsx).
const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

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

function renderConfiguration({ route }: { route: string }) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[route]}>
        <Routes>
          <Route path="/experiments/new/config" element={<ConfigurationPage />} />
          <Route path="/experiments/:id/config" element={<ConfigurationPage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

function mockFetchManifest(result: api.ManifestResponse) {
  vi.spyOn(api, "fetchManifest").mockResolvedValue(result);
}

describe("ConfigurationPage (Phase E1 shell)", () => {
  beforeEach(() => {
    navigate.mockClear();
    vi.spyOn(api, "getExperiment").mockResolvedValue(EXP);
  });
  afterEach(() => vi.restoreAllMocks());

  it("renders the Geometry and Sources regions", async () => {
    renderPage();
    expect(await screen.findByTestId("config-geometry-region")).toBeInTheDocument();
    expect(screen.getByTestId("config-sources-region")).toBeInTheDocument();
  });
});

describe("ConfigurationPage (first-run mode)", () => {
  beforeEach(() => {
    navigate.mockClear();
  });
  afterEach(() => {
    vi.restoreAllMocks();
    useDraftExperiment.setState({ path: "", patterns: {} });
  });

  test("Configuration first-run runs the manifest, hides geometry, gates Approve, then creates", async () => {
    useDraftExperiment.setState({ path: "/data/run42", patterns: {} });
    mockFetchManifest({ total: 4, matched: { image: 2, metadata: 1, integration: 0 }, unmatched: [{ file: "s2", miss: "metadata" }] });
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 9 } as any);
    renderConfiguration({ route: "/experiments/new/config" });   // first-run = draft route, no :id
    expect(screen.getByRole("button", { name: /approve/i })).toBeDisabled();        // while indexing
    expect(await screen.findByText(/2 image/i)).toBeInTheDocument();
    expect(screen.queryByText(/geometry/i)).toBeNull();                              // hidden first-run
    expect(screen.getByRole("button", { name: /approve/i })).toBeEnabled();
    fireEvent.click(screen.getByRole("button", { name: /approve/i }));
    await waitFor(() => expect(create).toHaveBeenCalledWith(expect.objectContaining({ path: "/data/run42" })));
    await waitFor(() => expect(navigate).toHaveBeenCalledWith("/experiments/9/corpus"));
  });
});
