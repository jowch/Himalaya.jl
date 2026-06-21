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
  manifest_path: null, created_at: "2026-04-12T00:00:00", q_units: "A^-1", beam_center_x: 1,
  beam_center_y: 1, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8,
  energy_kev_source: "prp", flight_path_m_source: "setup", beam_center_x_source: "setup",
  beam_center_y_source: "setup", pixel_size_um_source: "prp", q_units_source: "prp",
  last_scanned_at: null, scan_signature: null, ingest_status: "complete",
  image_pattern: null, metadata_pattern: null, integration_pattern: null,
  stats: { loads: 13, samples: 170, exposures: 682, sessions: 3, span_hours: 6.5, started_at: "2026-04-01T10:00:00" },
  review_count: 0,
  ...over,
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

  it("renders the counts from stats", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ stats: { loads: 13, samples: 170, exposures: 682, sessions: 3, span_hours: 6.5, started_at: "2026-04-01T10:00:00" } }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("170")).toBeInTheDocument();
    expect(screen.getByText("682")).toBeInTheDocument();
    expect(screen.getByText("13")).toBeInTheDocument();
  });

  it("shows 'up to date' chip when complete and no review_count", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ ingest_status: "complete", review_count: 0 }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("up to date")).toBeInTheDocument();
  });

  it("shows 'N to review' warning chip when review_count > 0", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ ingest_status: "complete", review_count: 4 }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("4 to review")).toBeInTheDocument();
  });

  it("shows 'scanning…' chip when ingest_status is scanning", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ ingest_status: "scanning", review_count: 0 }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("scanning…")).toBeInTheDocument();
  });

  it("groups experiments under year headings", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({
        id: 1, name: "Exp 2026",
        stats: { loads: 1, samples: 1, exposures: 1, sessions: 1, span_hours: 0, started_at: "2026-01-01T00:00:00" },
      }),
      EXP({
        id: 2, name: "Exp 2025",
        stats: { loads: 1, samples: 1, exposures: 1, sessions: 1, span_hours: 0, started_at: "2025-06-01T00:00:00" },
      }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-1");
    // Year headings are h2 elements; both 2026 and 2025 appear in the rail AND
    // as headings. Confirm the heading role specifically.
    expect(screen.getByRole("heading", { level: 2, name: "2026" })).toBeInTheDocument();
    expect(screen.getByRole("heading", { level: 2, name: "2025" })).toBeInTheDocument();
  });

  it("falls back to created_at year when stats.started_at is null", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({
        id: 5, name: "FallbackYear",
        created_at: "2024-03-15T00:00:00",
        stats: { loads: 0, samples: 0, exposures: 0, sessions: 0, span_hours: 0, started_at: null },
      }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-5");
    expect(screen.getByRole("heading", { level: 2, name: "2024" })).toBeInTheDocument();
  });

  it("shows 'indexing…' instead of counts when scanning", async () => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([
      EXP({ ingest_status: "scanning" }),
    ]);
    renderPage();
    await screen.findByTestId("experiment-card-7");
    expect(screen.getByText("indexing…")).toBeInTheDocument();
  });
});
