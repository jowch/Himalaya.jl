// test/ExperimentCorpusPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentCorpusPage } from "../src/print/pages/ExperimentCorpusPage";
import { useAppState } from "../src/state";
import * as api from "../src/api";
import type { Load } from "../src/api";

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

function renderAt(loads: Load[], processing = false) {
  if (processing) useAppState.setState({ ingestInFlight: { 7: { processed: 10, total: 100, status: "scanning" } } });
  else useAppState.setState({ ingestInFlight: null });
  vi.spyOn(api, "listLoads").mockResolvedValue(loads);
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

  it("hides the banner when nothing needs review", async () => {
    // Load with no flagged samples.
    const cleanLoad: Load = { ...LOAD_WITH_FLAG,
      samples: LOAD_WITH_FLAG.samples.map((s) => ({ ...s, flag: null })) };
    renderAt([cleanLoad]);
    await waitFor(() => expect(api.listLoads).toHaveBeenCalled());
    expect(screen.queryByTestId("grouping-review-banner")).toBeNull();
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
});
