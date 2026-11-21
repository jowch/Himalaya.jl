/**
 * AppRoutes — the single hoisted route table (T3.2: AppShell unified).
 * TopNav is the single shell; CorpusShell/CorpusTopbar/ExperimentTopNav deleted.
 */
import { describe, it, expect, beforeEach, vi, afterEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppRoutes } from "../src/print/shell/AppRoutes";
import * as api from "../src/api";

function makeQc() {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function renderRoutes(initialPath: string, initialIndex?: number) {
  const qc = makeQc();
  const entries = initialIndex !== undefined
    ? { initialEntries: [initialPath, "/"], initialIndex }
    : { initialEntries: [initialPath] };
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter {...entries}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("AppRoutes — single unified shell (TopNav, T3.2)", () => {
  beforeEach(() => {
    useAppState.setState({
      activeExperimentId: undefined,
      staleUrlContext: null,
    });
  });
  afterEach(() => vi.restoreAllMocks());

  it("mounts the single app shell at the stale catch-all (T3.2: unified shell)", async () => {
    renderRoutes("/totally/unknown");
    expect(await screen.findByTestId("app-shell")).toBeInTheDocument();
    expect(await screen.findByTestId("stale-url-page")).toBeInTheDocument();
    expect(screen.queryByTestId("corpus-shell")).toBeNull();
  });

  it.each(["/series", "/series/4", "/compare/all", "/experiments/7/compare/123"])(
    "treats the deleted %s as a stale URL",
    async (path) => {
      renderRoutes(path);
      expect(await screen.findByTestId("stale-url-page")).toBeInTheDocument();
    },
  );

  it("mounts the workspace outside the app shell", async () => {
    vi.spyOn(api, "listCorpusSamples").mockResolvedValue([]);
    vi.spyOn(api, "getExperiment").mockResolvedValue({ id: 7, name: "BL-19" } as Awaited<ReturnType<typeof api.getExperiment>>);
    renderRoutes("/experiments/7/workspace");
    expect(await screen.findByTestId("workspace")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();
    expect(await screen.findByText("No samples in this experiment.")).toBeInTheDocument();
  });
});

describe("AppRoutes — I4.4 index cutover redirects", () => {
  beforeEach(() => {
    useAppState.setState({
      activeExperimentId: undefined,
      activeSampleId: undefined,
      staleUrlContext: null,
    });
    // / now redirects to /experiments; mock the list so the home page renders.
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  });
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("bare / redirects to the experiments home (/experiments)", async () => {
    renderRoutes("/");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).not.toBeNull();
  });

  it("/index redirects to /experiments (T3.2)", async () => {
    renderRoutes("/index");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
  });

  it("/index/:experiment (no sample) redirects to /experiments (T3.2)", async () => {
    renderRoutes("/index/lipid");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
  });

  it("/index/:experiment/:sample falls back to /experiments when resolve 404s (T3.2)", async () => {
    vi.spyOn(global, "fetch").mockImplementation(() =>
      Promise.resolve({
        ok: false, status: 404,
        json: () => Promise.resolve({
          error: "not_found", missing: "sample", missing_value: "JC404",
          experiment_resolved: { id: 1, name: "lipid" }, sample_resolved: undefined,
        }),
      } as Response),
    );
    renderRoutes("/index/lipid/JC404");
    // T3.2: /samples redirects → /experiments home
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
  });

  it("an unknown path renders StaleUrlPage (#181 regression, #177)", async () => {
    renderRoutes("/foo/bar");
    expect(await screen.findByTestId("stale-url-page")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });
});

describe("AppRoutes — bare / always lands on experiments home (#77 / I4.4 / E1)", () => {
  beforeEach(() => {
    useAppState.setState({
      activeExperimentId: undefined,
      staleUrlContext: null,
    });
    // / redirects to /experiments; mock the list so the home page renders.
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  });
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("bare / lands on the experiments home gallery", async () => {
    useAppState.setState({ activeExperimentId: undefined });
    renderRoutes("/");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("bare / lands on the experiments home even with an active experiment set", async () => {
    useAppState.setState({ activeExperimentId: 7 });
    renderRoutes("/");
    expect(await screen.findByText("All experiments")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

});
