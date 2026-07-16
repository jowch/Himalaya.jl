// test/AppRoutes.ingestion.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/print/shell/AppRoutes";
import * as api from "../src/api";

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}><AppRoutes /></MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("AppRoutes ingestion tree (Phase E1)", () => {
  beforeEach(() => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
    vi.spyOn(api, "getExperiment").mockResolvedValue({ id: 7, name: "X", ingest_status: "complete" } as api.Experiment);
  });
  afterEach(() => vi.restoreAllMocks());

  it("/experiments renders the home gallery", async () => {
    renderAt("/experiments");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
  });

  it("/experiments/new renders the directory picker", () => {
    renderAt("/experiments/new");
    expect(screen.getByTestId("dirpicker-input")).toBeInTheDocument();
  });

  it("/experiments/:id/corpus mounts ExperimentShell chrome + corpus body", async () => {
    renderAt("/experiments/7/corpus");
    expect(await screen.findByTestId("experiment-shell")).toBeInTheDocument();
    // T3.2: AppShell renders TopNav; ExperimentShell is pure page content.
    // There is exactly ONE topnav in the DOM (from AppShell).
    expect(screen.getByTestId("topnav")).toBeInTheDocument();
    expect(screen.queryByTestId("corpus-topbar")).toBeNull();
  });

  it("/experiments/:id/grouping mounts GroupingReviewPage (E2)", async () => {
    vi.spyOn(api, "listLoads").mockResolvedValue([]);
    renderAt("/experiments/7/grouping");
    // GroupingReviewPage renders the "Check the grouping" heading + filter bar.
    expect(await screen.findByText("Check the grouping")).toBeInTheDocument();
    expect(screen.getByRole("textbox", { name: "Filter samples" })).toBeInTheDocument();
  });

  it("/ redirects to /experiments", async () => {
    renderAt("/");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
  });
});
