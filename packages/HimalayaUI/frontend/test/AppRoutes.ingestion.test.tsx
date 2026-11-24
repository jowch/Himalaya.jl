import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/shell/AppRoutes";
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

  it("/experiments/new renders the folder browser", () => {
    vi.spyOn(api, "listFolder").mockResolvedValue({ path: "/lab", parent: "/", folders: [], inside_setup: false });
    renderAt("/experiments/new");
    expect(screen.getByTestId("folder-browser-list")).toBeInTheDocument();
  });

  it("/ redirects to /experiments", async () => {
    renderAt("/");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
  });
});
