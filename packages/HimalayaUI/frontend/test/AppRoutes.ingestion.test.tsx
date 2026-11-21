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

  it.each(["corpus", "config", "grouping"])("/experiments/:id/%s goes to the workspace", async (retired) => {
    renderAt(`/experiments/7/${retired}`);
    expect(await screen.findByTestId("workspace")).toBeInTheDocument();
  });

  it("/experiments/new/config goes back to the new-experiment screen", () => {
    renderAt("/experiments/new/config");
    expect(screen.getByTestId("dirpicker-input")).toBeInTheDocument();
  });

  it("/ redirects to /experiments", async () => {
    renderAt("/");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
  });
});
