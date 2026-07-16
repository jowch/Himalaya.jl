// test/ExperimentCorpusPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, waitFor, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentCorpusPage } from "../src/print/pages/ExperimentCorpusPage";
import { useAppState } from "../src/state";
import * as api from "../src/api";
import type { Load, Experiment } from "../src/api";

// Minimal Experiment stub for the failed state.
const FAILED_EXPERIMENT: Experiment = {
  id: 7, name: "TestExp", description: null, path: "/data/exp1",
  data_dir: "/data/exp1/raw", analysis_dir: "/data/exp1/analysis",
  manifest_path: null, created_at: "2026-01-01T00:00:00Z",
  q_units: null, beam_center_x: null, beam_center_y: null,
  pixel_size_um: null, energy_kev: null, flight_path_m: null,
  energy_kev_source: "computed", flight_path_m_source: "computed",
  beam_center_x_source: "computed", beam_center_y_source: "computed",
  pixel_size_um_source: "computed", q_units_source: "computed",
  last_scanned_at: null, scan_signature: null, ingest_status: "failed",
  image_pattern: "*.tif", metadata_pattern: null, integration_pattern: null,
};

// Minimal load stub with a flagged merge discrepancy.
const LOAD_WITH_FLAG: Load = {
  load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null,
  frame_count: 4, note: null,
  samples: [
    { sample_id: 9, name: "HA85", slot_index: 1, grouping_source: "computed",
      name_source: "computed", merged_into_id: null,
      flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "HA86" },
      exposures: [] },
    { sample_id: 10, name: "HA86", slot_index: 2, grouping_source: "computed",
      name_source: "computed", merged_into_id: null, flag: null, exposures: [] },
    { sample_id: 11, name: "HA87", slot_index: 3, grouping_source: "computed",
      name_source: "computed", merged_into_id: null,
      flag: { kind: "split", split_at_index: 2, jump_from: 12.4, jump_to: 48.1 },
      exposures: [] },
    { sample_id: 12, name: "HA88", slot_index: 4, grouping_source: "computed",
      name_source: "computed", merged_into_id: null, flag: null, exposures: [] },
  ],
};

function renderAt(loads: Load[], processing = false, experiment?: Experiment) {
  if (processing) useAppState.setState({ ingestInFlight: { 7: { processed: 10, total: 100, status: "scanning" } } });
  else useAppState.setState({ ingestInFlight: null });
  vi.spyOn(api, "listLoads").mockResolvedValue(loads);
  // The corpus rows now read useCorpusSamples() + useCorpusExposures(); mock
  // both empty so the sheet renders with no rows (these tests assert the banner
  // + state machine, not the rows).
  vi.spyOn(api, "listCorpusSamples").mockResolvedValue([]);
  vi.spyOn(api, "listExposures").mockResolvedValue([]);
  if (experiment) vi.spyOn(api, "getExperiment").mockResolvedValue(experiment);
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7/corpus"]}>
        <Routes>
          <Route path="/experiments/:id/corpus" element={<ExperimentCorpusPage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("ExperimentCorpusPage (Phase E1)", () => {
  beforeEach(() => vi.restoreAllMocks());
  afterEach(() => { useAppState.setState({ ingestInFlight: null }); vi.restoreAllMocks(); });

  it("shows the grouping-review banner with the count + link when reviews are pending", async () => {
    renderAt([LOAD_WITH_FLAG]);
    // 2 flagged samples: LOAD_WITH_FLAG.samples[0] (merge) + [2] (split)
    const banner = await screen.findByTestId("grouping-review-banner");
    expect(banner).toHaveTextContent("2");
    expect(screen.getByTestId("grouping-review-link")).toHaveAttribute("href", "/experiments/7/grouping");
  });

  it("keeps the banner when nothing needs review, in the calm 'clear' state + same link", async () => {
    // Load with no flagged samples. The banner stays (data-state=clear) with
    // reassuring copy and the SAME entry point, so review is always one click away.
    const cleanLoad: Load = { ...LOAD_WITH_FLAG,
      samples: LOAD_WITH_FLAG.samples.map((s) => ({ ...s, flag: null })) };
    renderAt([cleanLoad]);
    const banner = await screen.findByTestId("grouping-review-banner");
    expect(banner).toHaveAttribute("data-state", "clear");
    expect(banner).toHaveTextContent(/settled/i);
    expect(screen.getByTestId("grouping-review-link")).toHaveAttribute("href", "/experiments/7/grouping");
  });

  it("renders the live-ingest placeholder while re-analyzing (analyzing status)", () => {
    // T4.3 state machine: `analyzing` (rescanning) → inline progress slot.
    // `scanning` (initial ingest) → GroupingReviewPage, NOT the slot.
    useAppState.setState({ ingestInFlight: { 7: { processed: 10, total: 100, status: "analyzing" } } });
    vi.spyOn(api, "listLoads").mockResolvedValue([]);
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments/7/corpus"]}>
          <Routes>
            <Route path="/experiments/:id/corpus" element={<ExperimentCorpusPage />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    expect(screen.getByTestId("live-ingest-slot")).toBeInTheDocument();
  });

  it("failed state: passes manifest unmatched + parsedCount from API to ScanFailedPage", async () => {
    // fetchManifest returns a real unmatched file + 3 parsed images.
    const manifestSpy = vi.spyOn(api, "fetchManifest").mockResolvedValue({
      total: 4,
      matched: { image: 3, metadata: 0, integration: 0 },
      unmatched: [{ file: "orphan.prp", miss: "metadata" }],
    });
    renderAt([], false, FAILED_EXPERIMENT);

    // The unmatched file should appear in the ScanFailedPage.
    await screen.findByText("orphan.prp");
    // The per-type pattern test input for the affected miss type renders.
    expect(screen.getByRole("textbox", { name: /metadata pattern/i })).toBeInTheDocument();
    // B1 (§5.5): the fetch must pass the stored leaf analysis_dir (4th arg), or
    // integration (.dat) is matched against data_dir and always reports unmatched.
    expect(manifestSpy).toHaveBeenCalledWith(
      "/data/exp1/raw", expect.anything(), undefined, "/data/exp1/analysis",
    );
  });

  it("failed state: 'Ingest N that parsed' → Confirm calls triggerScan", async () => {
    vi.spyOn(api, "fetchManifest").mockResolvedValue({
      total: 5,
      matched: { image: 3, metadata: 2, integration: 0 },
      unmatched: [{ file: "x.raw", miss: "integration" }],
    });
    const triggerScanSpy = vi.spyOn(api, "triggerScan").mockResolvedValue(undefined as never);
    renderAt([], false, FAILED_EXPERIMENT);

    // Wait for manifest to resolve so the "Ingest 3 that parsed" button appears.
    const ingestBtn = await screen.findByRole("button", { name: /ingest 3 that parsed/i });
    fireEvent.click(ingestBtn);
    // Confirm stage — mutation is async, wait for the spy to be called.
    fireEvent.click(screen.getByTestId("ingest-confirm-yes"));
    await waitFor(() => expect(triggerScanSpy).toHaveBeenCalled());
  });
});
