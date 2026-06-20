import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { PrintApp as App } from "../src/print/App";
import { useAppState } from "../src/state";
import { renderWithProviders } from "./test-utils";

// Plan §Phase 4 introduced URL routing for Compare. Wrap App renders in a
// MemoryRouter so route hooks resolve. I4.4 (#181): the Index surface at "/"
// is retired — "/" redirects to the corpus contact sheet (/samples). Tests
// pass an explicit entry so they land where they assert.
const renderAppAt = (path: string) =>
  renderWithProviders(
    <MemoryRouter initialEntries={[path]}>
      <App />
    </MemoryRouter>,
  );
const renderApp = () => renderAppAt("/");

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
      tutorialSeen: true,
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
      // Corpus list — "/" now lands on the /samples contact sheet (#181).
      "/api/samples": [
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
      // Series folio listing — a /compare* URL now redirects here (#177).
      "/api/series": [],
    });
  });

  it("redirects '/' to the Experiments home (ingestion-redesign IA)", async () => {
    // Ingestion redesign (E1): the new IA makes "/experiments" the home; a cold
    // "/" redirects there (was "/samples" under #181/#182). ExperimentsHomePage
    // sits OUTSIDE CorpusShell (so the experiment chrome never stacks on
    // CorpusTopbar), so corpus-shell is NOT mounted at "/". The "New experiment"
    // CTA renders unconditionally, so it is the stable landing assertion.
    renderApp();
    await waitFor(() =>
      expect(screen.getByTestId("new-experiment-cta")).toBeInTheDocument(),
      { timeout: 3000 },
    );
    expect(screen.queryByTestId("corpus-shell")).toBeNull();
    expect(screen.queryByTestId("app-shell")).toBeNull();
  });

  it("shows the onboarding overlay when no user is set", () => {
    useAppState.setState({ username: undefined });
    renderApp();
    expect(screen.getByTestId("onboarding-overlay")).toBeInTheDocument();
  });

  it("a /compare* URL redirects to the series folio (Compare retired, #177)", async () => {
    renderAppAt("/experiments/1/compare");
    await waitFor(() => expect(screen.getByTestId("folio-header")).toBeInTheDocument(),
      { timeout: 3000 });
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });
});
