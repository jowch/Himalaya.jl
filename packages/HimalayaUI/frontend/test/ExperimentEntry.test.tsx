import { describe, it, expect, vi, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ExperimentEntry, EXPERIMENT_GONE } from "../src/shell/ExperimentEntry";
import * as toast from "../src/lib/toast";
import { useAppState } from "../src/state";
import * as api from "../src/api";
import type { Experiment } from "../src/api";

const EXPERIMENT: Experiment = {
  id: 7, name: "TestExp", description: null, path: "/data/exp1",
  created_at: "2026-01-01T00:00:00Z",
  last_scanned_at: "2026-01-01T10:00:00Z", ingest_status: "complete",
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

  it("goes home with a toast for an experiment that no longer exists", async () => {
    const shown = vi.spyOn(toast, "showToast");
    vi.spyOn(api, "getExperiment").mockRejectedValue(new api.ApiError(404, "experiment not found", null));
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
    expect(await screen.findByTestId("where")).toHaveTextContent("/experiments");
    expect(shown).toHaveBeenCalledWith(EXPERIMENT_GONE, "warning");
  });
});
