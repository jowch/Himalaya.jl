import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import { App } from "../src/App";
import { useAppState } from "../src/state";
import { renderWithProviders } from "./test-utils";

function mockFetch(map: Record<string, unknown>): void {
  vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    const key = Object.keys(map).find((k) => url.endsWith(k));
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
        { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [] },
      ],
      "/api/samples/10/exposures": [],
      "/api/samples/10/messages": [],
    });
  });

  it("renders the three-card index page when user + scope are set", async () => {
    renderWithProviders(<App />);
    // Three-card grid + title button should all appear
    expect(await screen.findByTestId("three-card-grid")).toBeInTheDocument();
    expect(screen.getByTestId("plot-title")).toBeInTheDocument();
    expect(screen.getByTestId("tab-rocker")).toBeInTheDocument();
    // Title should include the experiment and sample name once the queries resolve
    await waitFor(() => expect(screen.getByText(/demo/)).toBeInTheDocument());
    await waitFor(() => expect(screen.getByText("s1")).toBeInTheDocument());
  });

  it("shows the onboarding overlay when no user is set", () => {
    useAppState.setState({ username: undefined });
    renderWithProviders(<App />);
    expect(screen.getByTestId("onboarding-overlay")).toBeInTheDocument();
  });

  it("switching the tab rocker changes the active page to Compare", async () => {
    renderWithProviders(<App />);
    const cmpTab = await screen.findByTestId("tab-compare");
    cmpTab.click();
    await waitFor(() => expect(screen.getByTestId("compare-page")).toBeInTheDocument());
  });
});
