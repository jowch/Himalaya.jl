import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/components/AppRoutes";
import { useAppState } from "../src/state";

const SERIES = {
  id: 5, title: "LL37 titration", description: null, content_hash: "h",
  created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
  forked_at_hash: null, forked_from_title: null, view_grouping_mode: null,
  view_show_peak_ticks: null, view_show_peak_labels: null,
  ordering_variable: null, order_rule: "ascending", state: "committed",
  members: [], samples: [],
};

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({
    activeSampleId: undefined, activeExposureId: undefined,
    activePage: "index", activeExperimentId: undefined,
    username: "tester", theme: "dark",
  });
  vi.stubGlobal("fetch", vi.fn(async (url: string) => {
    const u = String(url);
    const body = u.includes("/api/series/5") ? SERIES : [];
    return new Response(JSON.stringify(body), {
      status: 200, headers: { "content-type": "application/json" },
    });
  }));
});

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}><AppRoutes /></MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("/series/:id routing", () => {
  it("mounts SeriesBuilderPage under the corpus shell", async () => {
    renderAt("/series/5");
    expect(await screen.findByTestId("series-builder-page")).toBeInTheDocument();
    expect(screen.getByTestId("corpus-topbar")).toBeInTheDocument();
  });
});
