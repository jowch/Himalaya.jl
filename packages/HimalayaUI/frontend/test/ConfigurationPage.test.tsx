// test/ConfigurationPage.test.tsx
import { describe, it, test, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ConfigurationPage } from "../src/print/pages/ConfigurationPage";
import * as api from "../src/api";
import { useDraftExperiment } from "../src/lib/draftExperiment";

// Module-level navigate mock -- canonical pattern (mirrors NewExperimentPage.test.tsx).
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
    useDraftExperiment.getState().clear();
  });
  afterEach(() => {
    vi.restoreAllMocks();
    useDraftExperiment.getState().clear();
  });

  test("first-run resolves the root, shows geometry, gates Approve, then creates with confirmed values", async () => {
    // The picker committed a root; the page resolves it structurally.
    useDraftExperiment.getState().setRoot("/data/run42");
    vi.spyOn(api, "resolveLayout").mockResolvedValue({
      name: "run42",
      data_dir: "/data/run42/data",
      analysis_dir: "/data/run42/analysis",
      setup_file: "/data/run42/analysis/setup_info_x.txt",
      setup_ambiguous: false,
    });
    mockFetchManifest({
      total: 4,
      matched: { image: 2, metadata: 1, integration: 0 },
      unmatched: [{ file: "s2", miss: "metadata" }],
      geometry: {
        beam_center_x: 421.3, beam_center_x_source: "setup",
        beam_center_y: 836.7, beam_center_y_source: "setup",
        flight_path_m: 1.8095, flight_path_m_source: "setup",
        pixel_size_um: 172.0, pixel_size_um_source: "prp",
        energy_kev: 9.0, energy_kev_source: "prp",
      },
      matched_files: ["JC_001.tif", "JC_002.tif"],
    });
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 9 } as any);
    renderConfiguration({ route: "/experiments/new/config" });   // first-run = draft route, no :id

    // Resolve → manifest → body. The matched count + geometry are now visible.
    expect(await screen.findByText(/2 matched/i)).toBeInTheDocument();
    expect(screen.getByText(/auto-derived/i)).toBeInTheDocument();   // the Geometry card heading
    expect(screen.getByText(/421\.3/)).toBeInTheDocument();
    // The Name field is prefilled from the resolver (editable).
    expect(
      (screen.getByTestId("config-name").querySelector("input") as HTMLInputElement).value,
    ).toBe("run42");

    // Approve enabled once indexing resolved; creates with the CONFIRMED values.
    await waitFor(() => expect(screen.getByRole("button", { name: /approve/i })).toBeEnabled());
    fireEvent.click(screen.getByRole("button", { name: /approve/i }));
    await waitFor(() => expect(create).toHaveBeenCalledWith(expect.objectContaining({
      name: "run42", data_dir: "/data/run42/data",
    })));
    // Approve lands on the combined scan + grouping-review surface.
    await waitFor(() => expect(navigate).toHaveBeenCalledWith("/experiments/9/grouping"));
  });
});
