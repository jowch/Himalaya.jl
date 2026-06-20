// test/ExperimentsHomePage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentsHomePage } from "../src/print/pages/ExperimentsHomePage";
import * as api from "../src/api";

const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

function renderPage() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter><ExperimentsHomePage /></MemoryRouter>
    </QueryClientProvider>,
  );
}

const EXP = (over: Partial<api.Experiment>): api.Experiment => ({
  id: 7, name: "SSRL · 1p7m", description: null, path: "/d", data_dir: "/d/data", analysis_dir: "/d/an",
  manifest_path: null, created_at: "2026-04-12", q_units: "A^-1", beam_center_x: 1,
  beam_center_y: 1, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8,
  energy_kev_source: "prp", flight_path_m_source: "setup", beam_center_x_source: "setup",
  beam_center_y_source: "setup", pixel_size_um_source: "prp", q_units_source: "prp",
  last_scanned_at: null, scan_signature: null, ingest_status: "complete",
  image_pattern: null, metadata_pattern: null, integration_pattern: null, ...over,
});

describe("ExperimentsHomePage (Phase E1)", () => {
  beforeEach(() => { navigate.mockClear(); });
  afterEach(() => vi.restoreAllMocks());

  it("lists experiment cards", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({}), EXP({ id: 8, name: "JC plate" })]);
    renderPage();
    expect(await screen.findByText("SSRL · 1p7m")).toBeInTheDocument();
    expect(screen.getByText("JC plate")).toBeInTheDocument();
  });

  it("routes to the corpus tab on card click", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({})]);
    renderPage();
    fireEvent.click(await screen.findByTestId("experiment-card-7"));
    expect(navigate).toHaveBeenCalledWith("/experiments/7/corpus");
  });

  it("routes to /experiments/new on the New CTA", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([EXP({})]);
    renderPage();
    fireEvent.click(await screen.findByTestId("new-experiment-cta"));
    expect(navigate).toHaveBeenCalledWith("/experiments/new");
  });

  it("shows an empty state when there are no experiments", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    renderPage();
    expect(await screen.findByText(/no experiments yet/i)).toBeInTheDocument();
  });
});
