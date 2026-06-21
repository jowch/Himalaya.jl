/**
 * AppRoutes — the single hoisted route table (T3.2: AppShell unified).
 * TopNav is the single shell; CorpusShell/CorpusTopbar/ExperimentTopNav deleted.
 * Includes the relocated #77 compare-sync tests (formerly AppShell.test.tsx).
 */
import { describe, it, expect, beforeEach, vi, afterEach } from "vitest";
import { render, screen, waitFor, fireEvent } from "@testing-library/react";
import { MemoryRouter, useNavigate } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppRoutes } from "../src/print/shell/AppRoutes";
import * as api from "../src/api";
import type { ResolveSuccess } from "../src/api";

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
      resolving: false,
    });
    vi.spyOn(api, "listSeries").mockResolvedValue([]);
  });
  afterEach(() => vi.restoreAllMocks());

  it("mounts the single app shell at /series", async () => {
    renderRoutes("/series");
    expect(await screen.findByTestId("app-shell")).toBeInTheDocument();
    expect(screen.getByTestId("folio-header")).toBeInTheDocument();
    expect(screen.queryByTestId("corpus-shell")).toBeNull();
  });

  it("mounts the single app shell at the stale catch-all (T3.2: unified shell)", async () => {
    renderRoutes("/totally/unknown");
    expect(await screen.findByTestId("app-shell")).toBeInTheDocument();
    expect(await screen.findByTestId("stale-url-page")).toBeInTheDocument();
    expect(screen.queryByTestId("corpus-shell")).toBeNull();
  });

  it("redirects a /compare* URL to the series folio (Compare retired, #177)", async () => {
    renderRoutes("/compare/all");
    expect(await screen.findByTestId("folio-header")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("redirects an experiment-scoped /compare URL to the series folio (#177)", async () => {
    renderRoutes("/experiments/7/compare/123");
    expect(await screen.findByTestId("folio-header")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("does not flag /series as a stale path", async () => {
    renderRoutes("/series");
    await screen.findByTestId("folio-header");
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  it("keeps the single app shell mounted across the series↔stale boundary (T3.2)", async () => {
    // T3.2: one AppShell. Stays mounted whether the route is a matched surface
    // or the `*` stale catch-all; only the body swaps.
    function NavButtons(): JSX.Element {
      const navigate = useNavigate();
      return (
        <>
          <button data-testid="go-series" onClick={() => navigate("/series")}>series</button>
          <button data-testid="go-stale" onClick={() => navigate("/totally/unknown")}>stale</button>
        </>
      );
    }

    const qc = makeQc();
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/series"]}>
          <NavButtons />
          <AppRoutes />
        </MemoryRouter>
      </QueryClientProvider>,
    );

    // Start: app shell up, folio rendered, no legacy corpus-shell.
    expect(await screen.findByTestId("app-shell")).toBeInTheDocument();
    expect(screen.getByTestId("folio-header")).toBeInTheDocument();
    expect(screen.queryByTestId("corpus-shell")).toBeNull();

    // Cross to the stale catch-all — app shell stays; body becomes StaleUrlPage.
    fireEvent.click(screen.getByTestId("go-stale"));
    await waitFor(() => {
      expect(screen.getByTestId("app-shell")).toBeInTheDocument();
      expect(screen.getByTestId("stale-url-page")).toBeInTheDocument();
    });
    expect(screen.queryByTestId("corpus-shell")).toBeNull();

    // Cross back to series route — folio returns under the same shell.
    fireEvent.click(screen.getByTestId("go-series"));
    await waitFor(() => {
      expect(screen.getByTestId("app-shell")).toBeInTheDocument();
      expect(screen.getByTestId("folio-header")).toBeInTheDocument();
    });
    expect(screen.queryByTestId("corpus-shell")).toBeNull();
  });
});

describe("AppRoutes — I4.4 index cutover redirects", () => {
  beforeEach(() => {
    useAppState.setState({
      activeExperimentId: undefined,
      activeSampleId: undefined,
      staleUrlContext: null,
      resolving: false,
    });
    // / now redirects to /experiments; mock the list so the home page renders.
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  });
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("bare / redirects to the experiments home (/experiments)", async () => {
    renderRoutes("/");
    expect(await screen.findByText("Your experiments")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).not.toBeNull();
  });

  it("/index redirects to /experiments (T3.2)", async () => {
    renderRoutes("/index");
    expect(await screen.findByText("Your experiments")).toBeInTheDocument();
  });

  it("/index/:experiment (no sample) redirects to /experiments (T3.2)", async () => {
    renderRoutes("/index/lipid");
    expect(await screen.findByText("Your experiments")).toBeInTheDocument();
  });

  it("/index/:experiment/:sample resolves the slug then redirects to /sample/:id", async () => {
    const body: ResolveSuccess = {
      experiment_id: 1, experiment_name: "lipid",
      sample_id: 10, sample_name: "JC001",
      exposure_id: undefined, exposure_filename: undefined,
    };
    vi.spyOn(global, "fetch").mockImplementation(() =>
      Promise.resolve({
        ok: true, status: 200, json: () => Promise.resolve(body),
      } as Response),
    );
    renderRoutes("/index/lipid/JC001");
    // Phase-4 cutover: /sample/:id now serves the greenfield FocusPage. The
    // single ResolveSuccess fetch mock doesn't satisfy the page's corpus-sample
    // query, so it renders its `focus-not-found` body — which only mounts on
    // the Focus route. That the Focus route (not /experiments) was reached is
    // the assertion this slug-redirect test makes.
    expect(await screen.findByTestId("focus-not-found")).toBeInTheDocument();
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
    expect(await screen.findByText("Your experiments")).toBeInTheDocument();
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
      resolving: false,
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
    expect(await screen.findByText("Your experiments")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("bare / lands on the experiments home even with an active experiment set", async () => {
    useAppState.setState({ activeExperimentId: 7 });
    renderRoutes("/");
    expect(await screen.findByText("Your experiments")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("a compare URL redirects to the series folio (Compare retired, #177)", async () => {
    useAppState.setState({ activeExperimentId: 7 });
    renderRoutes("/experiments/7/compare/123");
    await waitFor(() => {
      expect(screen.getByTestId("folio-header")).toBeInTheDocument();
    });
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });
});
