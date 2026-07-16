// test/ConfigurationPage.test.tsx
import { describe, it, test, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor, within } from "@testing-library/react";
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
      image_pattern: null,
      metadata_pattern: null,
      integration_pattern: null,
    });
    mockFetchManifest({
      total: 6,
      matched: { image: 2, metadata: 2, integration: 2 },   // clean triple: Approve enables
      unmatched: [],
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
    // The title reads as text at rest (prefilled from the resolver); the pencil
    // opens the rename field, prefilled with the same value.
    expect(screen.getByTestId("config-name")).toHaveTextContent("run42");
    fireEvent.click(screen.getByRole("button", { name: /rename experiment/i }));
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

  test("pattern fields use the Edit mechanism: read-only value, Edit reveals the input", async () => {
    useDraftExperiment.getState().setRoot("/data/run42");
    vi.spyOn(api, "resolveLayout").mockResolvedValue({
      name: "run42", data_dir: "/data/run42/data", analysis_dir: "/data/run42/analysis",
      setup_file: "/data/run42/analysis/setup_info_x.txt", setup_ambiguous: false,
      image_pattern: null, metadata_pattern: null, integration_pattern: null,
    });
    mockFetchManifest({
      total: 2, matched: { image: 2, metadata: 2, integration: 2 }, unmatched: [],
      geometry: {
        beam_center_x: 421.3, beam_center_x_source: "setup",
        beam_center_y: 836.7, beam_center_y_source: "setup",
        flight_path_m: 1.8095, flight_path_m_source: "setup",
        pixel_size_um: 172.0, pixel_size_um_source: "prp",
        energy_kev: 9.0, energy_kev_source: "prp",
      },
      matched_files: ["JC_001.tif"],
    });
    renderConfiguration({ route: "/experiments/new/config" });
    await screen.findByText(/2 matched/i);

    // At rest the pattern is read-only text + an Edit affordance (no input yet).
    expect(screen.getByText("{name}.tif")).toBeInTheDocument();
    expect(screen.queryByTestId("config-image-pattern")).toBeNull();
    // The Exposure-pattern row's Edit reveals the prefilled input.
    const row = screen.getByText("Exposure pattern").closest("div")!;
    fireEvent.click(within(row).getByRole("button", { name: /edit/i }));
    expect(
      (screen.getByTestId("config-image-pattern").querySelector("input") as HTMLInputElement).value,
    ).toBe("{name}.tif");
  });

  test("geometry fields use the Edit mechanism and an edit is sent as a create override", async () => {
    useDraftExperiment.getState().setRoot("/data/run42");
    vi.spyOn(api, "resolveLayout").mockResolvedValue({
      name: "run42", data_dir: "/data/run42/data", analysis_dir: "/data/run42/analysis",
      setup_file: "/data/run42/analysis/setup_info_x.txt", setup_ambiguous: false,
      image_pattern: null, metadata_pattern: null, integration_pattern: null,
    });
    mockFetchManifest({
      total: 2, matched: { image: 2, metadata: 2, integration: 2 }, unmatched: [],
      geometry: {
        beam_center_x: 421.3, beam_center_x_source: "setup",
        beam_center_y: 836.7, beam_center_y_source: "setup",
        flight_path_m: 1.8095, flight_path_m_source: "setup",
        pixel_size_um: 172.0, pixel_size_um_source: "prp",
        energy_kev: 9.0, energy_kev_source: "prp",
      },
      matched_files: ["JC_001.tif"],
    });
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 9 } as any);
    renderConfiguration({ route: "/experiments/new/config" });
    await screen.findByText(/2 matched/i);

    // Energy reads as text at rest (9.0 keV); its Edit opens a numeric input.
    expect(screen.getByText("9.0 keV")).toBeInTheDocument();
    const energyRow = screen.getByText("Energy").closest("div")!;
    fireEvent.click(within(energyRow).getByRole("button", { name: /edit/i }));
    const input = within(energyRow).getByLabelText("Energy") as HTMLInputElement;
    expect(input.value).toBe("9");
    fireEvent.change(input, { target: { value: "11.2" } });
    fireEvent.keyDown(input, { key: "Enter" });
    // Now reads the edited value + "edited" chip.
    expect(screen.getByText("11.2 keV")).toBeInTheDocument();

    await waitFor(() => expect(screen.getByRole("button", { name: /approve/i })).toBeEnabled());
    fireEvent.click(screen.getByRole("button", { name: /approve/i }));
    // The whole preview geometry is committed; the edited field is source='user',
    // the rest keep their derived sources.
    await waitFor(() => expect(create).toHaveBeenCalledWith(expect.objectContaining({
      geometry: expect.objectContaining({
        energy_kev: 11.2, energy_kev_source: "user",
        beam_center_x: 421.3, beam_center_x_source: "setup",
        flight_path_m: 1.8095, flight_path_m_source: "setup",
      }),
    })));
  });

  test("beam center edits both x and y in one override (Done commits the pair)", async () => {
    useDraftExperiment.getState().setRoot("/data/run42");
    vi.spyOn(api, "resolveLayout").mockResolvedValue({
      name: "run42", data_dir: "/data/run42/data", analysis_dir: "/data/run42/analysis",
      setup_file: "/data/run42/analysis/setup_info_x.txt", setup_ambiguous: false,
      image_pattern: null, metadata_pattern: null, integration_pattern: null,
    });
    mockFetchManifest({
      total: 2, matched: { image: 2, metadata: 2, integration: 2 }, unmatched: [],
      geometry: {
        beam_center_x: 421.3, beam_center_x_source: "setup",
        beam_center_y: 836.7, beam_center_y_source: "setup",
        flight_path_m: 1.8095, flight_path_m_source: "setup",
        pixel_size_um: 172.0, pixel_size_um_source: "prp",
        energy_kev: 9.0, energy_kev_source: "prp",
      },
      matched_files: ["JC_001.tif"],
    });
    const create = vi.spyOn(api, "createExperiment").mockResolvedValue({ id: 9 } as any);
    renderConfiguration({ route: "/experiments/new/config" });
    await screen.findByText(/2 matched/i);

    const bcRow = screen.getByText("Beam center").closest("div")!;
    fireEvent.click(within(bcRow).getByRole("button", { name: /edit/i }));
    fireEvent.change(within(bcRow).getByLabelText("Beam center X"), { target: { value: "400" } });
    fireEvent.change(within(bcRow).getByLabelText("Beam center Y"), { target: { value: "800" } });
    fireEvent.click(within(bcRow).getByRole("button", { name: /done/i }));
    expect(screen.getByText("400.0, 800.0 px")).toBeInTheDocument();

    await waitFor(() => expect(screen.getByRole("button", { name: /approve/i })).toBeEnabled());
    fireEvent.click(screen.getByRole("button", { name: /approve/i }));
    await waitFor(() => expect(create).toHaveBeenCalledWith(expect.objectContaining({
      geometry: expect.objectContaining({
        beam_center_x: 400, beam_center_x_source: "user",
        beam_center_y: 800, beam_center_y_source: "user",
      }),
    })));
  });

  test("D4: a whole type matching zero everywhere hard-blocks Approve", async () => {
    useDraftExperiment.getState().setRoot("/data/run42");
    vi.spyOn(api, "resolveLayout").mockResolvedValue({
      name: "run42",
      data_dir: "/data/run42/data",
      analysis_dir: "/data/run42/analysis",
      setup_file: "/data/run42/analysis/setup_info_x.txt",
      setup_ambiguous: false,
      image_pattern: null, metadata_pattern: null, integration_pattern: null,
    });
    mockFetchManifest({
      total: 3,
      matched: { image: 3, metadata: 3, integration: 0 },   // integration matched nowhere
      unmatched: [],
      geometry: {
        beam_center_x: 421.3, beam_center_x_source: "setup",
        beam_center_y: 836.7, beam_center_y_source: "setup",
        flight_path_m: 1.8095, flight_path_m_source: "setup",
        pixel_size_um: 172.0, pixel_size_um_source: "prp",
        energy_kev: 9.0, energy_kev_source: "prp",
      },
      matched_files: ["JC_001.tif"],
    });
    renderConfiguration({ route: "/experiments/new/config" });

    // Headline swaps to the block message; Approve stays disabled; no `/0` shown.
    expect(await screen.findByTestId("manifest-block")).toHaveTextContent(/No integration matched/i);
    expect(screen.getByRole("button", { name: /approve/i })).toBeDisabled();
  });
});
