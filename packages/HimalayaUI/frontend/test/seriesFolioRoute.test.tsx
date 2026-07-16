import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/print/shell/AppRoutes";
import { useAppState } from "../src/state";

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeSampleId: undefined,
    activeExposureId: undefined,
    activeExperimentId: undefined,
    username: "tester",
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
  it("mounts SeriesFolioPage under the unified app shell (T3.2)", async () => {
    renderAt("/series");
    expect(await screen.findByTestId("folio-header")).toBeInTheDocument();
    // T3.2: proves it is under the unified AppShell (TopNav), not the legacy CorpusShell
    expect(screen.getByTestId("topnav")).toBeInTheDocument();
  });
});
