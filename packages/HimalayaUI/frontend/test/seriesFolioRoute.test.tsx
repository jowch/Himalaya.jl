import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/components/AppRoutes";
import { useAppState } from "../src/state";

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeSampleId: undefined,
    activeExposureId: undefined,
    activePage: "index",
    activeExperimentId: undefined,
    username: "tester",
    theme: "dark",
  });
  vi.stubGlobal(
    "fetch",
    // Every endpoint (incl. /api/series) returns an empty-but-valid array.
    vi.fn(async () =>
      new Response("[]", {
        status: 200,
        headers: { "content-type": "application/json" },
      }),
    ),
  );
});

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("/series routing", () => {
  it("mounts SeriesFolioPage under the corpus shell", async () => {
    renderAt("/series");
    expect(await screen.findByTestId("series-folio-page")).toBeInTheDocument();
    // proves it is under CorpusShell (the corpus topbar), not the legacy AppShell
    expect(screen.getByTestId("corpus-topbar")).toBeInTheDocument();
  });
});
