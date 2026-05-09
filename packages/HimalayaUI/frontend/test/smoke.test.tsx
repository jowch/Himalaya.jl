import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { App } from "../src/App";
import { useAppState } from "../src/state";
import { renderWithProviders } from "./test-utils";

// Plan §Phase 4 introduced URL routing for Compare. Wrap App renders in a
// MemoryRouter so route hooks resolve. We start at "/" so the existing
// Zustand-driven IndexPage smoke assertions still pass.
const renderApp = () =>
  renderWithProviders(
    <MemoryRouter initialEntries={["/"]}>
      <App />
    </MemoryRouter>,
  );

function mockFetch(map: Record<string, unknown>): void {
  vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    // Strip query string so map keys can match by pathname alone.
    const path = url.split("?")[0] ?? url;
    const key = Object.keys(map).find((k) => path.endsWith(k));
    if (!key) return new Response("not found", { status: 404 });
    return new Response(JSON.stringify(map[key]), {
      status: 200, headers: { "Content-Type": "application/json" },
    });
  });
}

describe("App smoke", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    localStorage.clear();
    useAppState.setState({
      username: "alice",
      activeExperimentId: 1,
      activeSampleId: 10,
      activeExposureId: undefined,
      activePage: "index",
      tutorialSeen: true,
      theme: "dark",
      hoveredIndexId: undefined,
      navModalOpen: false,
    });
    mockFetch({
      "/api/experiments": [{
        id: 1, name: "demo", path: "/x", data_dir: "/x/data",
        analysis_dir: "/x/analysis", manifest_path: null,
        created_at: "2026-04-22T00:00:00Z",
      }],
      "/api/experiments/1": {
        id: 1, name: "demo", path: "/x", data_dir: "/x/data",
        analysis_dir: "/x/analysis", manifest_path: null,
        created_at: "2026-04-22T00:00:00Z",
      },
      "/api/experiments/1/samples": [
        { id: 10, experiment_id: 1, display_name: null, name: "s1", notes: null, tags: [] },
      ],
      "/api/samples/10/exposures": [],
      "/api/samples/10/messages": [],
      // Cold-mount root redirect path: useStateFromUrl falls through to
      // resolve-by-id when the TanStack cache is cold (test env never
      // hydrates synchronously). Without this, the redirect 404s and
      // /index wipes the seeded activeExperimentId / activeSampleId.
      "/api/resolve": {
        experiment_id: 1, experiment_name: "demo",
        sample_id: 10, sample_name: "s1",
        exposure_id: undefined, exposure_filename: undefined,
      },
    });
  });

  it("renders the three-card index page when user + scope are set", async () => {
    renderApp();
    // Three-card grid + title button should all appear once the URL→state
    // resolve finishes. Cold-mount path goes through ResolvingFallback
    // briefly because the TanStack cache is empty when useStateFromUrl
    // first fires; we wait for the resolve to land before asserting.
    await waitFor(() =>
      expect(screen.getByTestId("workspace-grid")).toBeInTheDocument(),
      { timeout: 3000 },
    );
    expect(screen.getByTestId("plot-title")).toBeInTheDocument();
    expect(screen.getByTestId("tab-rocker")).toBeInTheDocument();
    // Title should include the experiment and sample name once the queries resolve
    await waitFor(() => expect(screen.getByText(/demo/)).toBeInTheDocument(),
      { timeout: 3000 });
    await waitFor(() => expect(screen.getByText("s1")).toBeInTheDocument(),
      { timeout: 3000 });
  });

  it("shows the onboarding overlay when no user is set", () => {
    useAppState.setState({ username: undefined });
    renderApp();
    expect(screen.getByTestId("onboarding-overlay")).toBeInTheDocument();
  });

  it("switching the tab rocker changes the active page to Compare", async () => {
    renderApp();
    const cmpTab = await screen.findByTestId("tab-compare");
    cmpTab.click();
    await waitFor(() => expect(screen.getByTestId("compare-page")).toBeInTheDocument());
  });
});
