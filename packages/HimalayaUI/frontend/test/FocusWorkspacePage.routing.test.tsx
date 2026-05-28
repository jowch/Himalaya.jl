import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/components/AppRoutes";
import { useAppState } from "../src/state";

const CORPUS = [
  { id: 33, experiment_id: 1, name: "JC033", display_name: "Sample 33",
    notes: null, tags: [], q_units: "A-1" },
];

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
    vi.fn(async (url: string) => {
      const u = String(url);
      if (u.includes("/api/samples")) {
        return new Response(JSON.stringify(CORPUS), {
          status: 200, headers: { "content-type": "application/json" },
        });
      }
      // A trace is an object ({q, I, sigma}), never a bare array — return a
      // valid empty-but-shaped trace so PlotCard/TraceViewer never read fields
      // off an array.
      if (u.includes("/trace")) {
        return new Response(JSON.stringify({ q: [], I: [], sigma: [] }), {
          status: 200, headers: { "content-type": "application/json" },
        });
      }
      // exposures / experiments / everything else: empty but valid
      return new Response("[]", {
        status: 200, headers: { "content-type": "application/json" },
      });
    }),
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

describe("/sample/:sampleId routing", () => {
  it("mounts FocusWorkspacePage under the corpus shell", async () => {
    renderAt("/sample/33");
    expect(await screen.findByTestId("focus-workspace-page")).toBeInTheDocument();
    // proves it is under CorpusShell (the corpus topbar), not the legacy AppShell
    expect(screen.getByTestId("corpus-topbar")).toBeInTheDocument();
  });

  it("drives activeSampleId from the route param end-to-end", async () => {
    renderAt("/sample/33");
    await waitFor(() =>
      expect(useAppState.getState().activeSampleId).toBe(33),
    );
  });
});
