/**
 * AppRoutes — the single hoisted route table. After I5.1 (#182) there is a
 * SINGLE shell (CorpusShell): every route — including the `*` stale catch-all —
 * mounts under it. The legacy AppShell + dual-nav model are retired. Includes
 * the relocated #77 compare-sync tests (formerly AppShell.test.tsx).
 */
import { describe, it, expect, beforeEach, vi, afterEach } from "vitest";
import { render, screen, waitFor, fireEvent } from "@testing-library/react";
import { MemoryRouter, useNavigate } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppRoutes } from "../src/print/shell/AppRoutes";
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

describe("AppRoutes — single-shell route table", () => {
  beforeEach(() => {
    // Reset the ephemeral URL-resolution fields too — a prior test that
    // parked the store on a stale path must not leak into the next.
    useAppState.setState({
      activeExperimentId: undefined,
      staleUrlContext: null,
      resolving: false,
    });
  });

  it("mounts the corpus shell at /samples", async () => {
    renderRoutes("/samples");
    expect(await screen.findByTestId("corpus-shell")).toBeInTheDocument();
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();
  });

  it("mounts the corpus shell at the stale catch-all (I5.1: single shell)", async () => {
    // I5.1 (#182): AppShell is retired. The `*` catch-all now mounts under
    // CorpusShell and renders the StaleUrlPage body — there is no second shell.
    renderRoutes("/totally/unknown");
    expect(await screen.findByTestId("corpus-shell")).toBeInTheDocument();
    expect(await screen.findByTestId("stale-url-page")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();
  });

  it("redirects a /compare* URL to the series folio (Compare retired, #177)", async () => {
    renderRoutes("/compare/all");
    expect(await screen.findByTestId("folio-header")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("redirects an experiment-scoped /compare URL to the series folio (#177)", async () => {
    renderRoutes("/experiments/7/compare/123");
    expect(await screen.findByTestId("folio-header")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("does not flag /samples as a stale path", async () => {
    // The stale classifier (useStateFromUrl) runs only inside the `*`
    // catch-all body (PageBody) — a matched corpus route never mounts it, so
    // /samples cannot be parsed as a stale path.
    renderRoutes("/samples");
    await screen.findByTestId("corpus-shell");
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  it("keeps the single corpus shell mounted across the samples↔stale boundary", async () => {
    // I5.1 (#182): one shell. CorpusShell stays mounted whether the route is a
    // matched corpus surface or the `*` stale catch-all; only the body swaps.
    function NavButtons(): JSX.Element {
      const navigate = useNavigate();
      return (
        <>
          <button data-testid="go-samples" onClick={() => navigate("/samples")}>samples</button>
          <button data-testid="go-stale" onClick={() => navigate("/totally/unknown")}>stale</button>
        </>
      );
    }

    const qc = makeQc();
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/samples"]}>
          <NavButtons />
          <AppRoutes />
        </MemoryRouter>
      </QueryClientProvider>,
    );

    // Start: corpus shell up, samples body rendered, no legacy shell.
    expect(await screen.findByTestId("corpus-shell")).toBeInTheDocument();
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();

    // Cross to the stale catch-all — corpus shell stays; body becomes StaleUrlPage.
    fireEvent.click(screen.getByTestId("go-stale"));
    await waitFor(() => {
      expect(screen.getByTestId("corpus-shell")).toBeInTheDocument();
      expect(screen.getByTestId("stale-url-page")).toBeInTheDocument();
    });
    expect(screen.queryByTestId("app-shell")).toBeNull();

    // Cross back to corpus route — samples body returns under the same shell.
    fireEvent.click(screen.getByTestId("go-samples"));
    await waitFor(() => {
      expect(screen.getByTestId("corpus-shell")).toBeInTheDocument();
      expect(screen.getByTestId("samples-page")).toBeInTheDocument();
    });
    expect(screen.queryByTestId("app-shell")).toBeNull();
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
  });
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("bare / redirects to the corpus contact sheet (/samples)", async () => {
    renderRoutes("/");
    expect(await screen.findByTestId("samples-page")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();
  });

  it("/index redirects to /samples", async () => {
    renderRoutes("/index");
    expect(await screen.findByTestId("samples-page")).toBeInTheDocument();
  });

  it("/index/:experiment (no sample) redirects to /samples", async () => {
    renderRoutes("/index/lipid");
    expect(await screen.findByTestId("samples-page")).toBeInTheDocument();
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
    // the Focus route. That the Focus route (not /samples) was reached is the
    // assertion this slug-redirect test makes.
    expect(await screen.findByTestId("focus-not-found")).toBeInTheDocument();
  });

  it("/index/:experiment/:sample falls back to /samples when resolve 404s", async () => {
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
    expect(await screen.findByTestId("samples-page")).toBeInTheDocument();
  });

  it("an unknown path renders StaleUrlPage (#181 regression, #177)", async () => {
    // A stale path cleanly renders "Page not found" via the `*` catch-all
    // PageBody, now mounted under the single CorpusShell (I5.1, #182).
    renderRoutes("/foo/bar");
    expect(await screen.findByTestId("stale-url-page")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });
});

describe("AppRoutes — bare / always lands on the corpus (#77 / I4.4)", () => {
  beforeEach(() => {
    useAppState.setState({
      activeExperimentId: undefined,
      staleUrlContext: null,
      resolving: false,
    });
  });

  // I4.4 (#181): the #77 "empty PageBody at /" risk is eliminated differently
  // now. Bare `/` is a standalone redirect to /samples, so a cold `/` can never
  // strand the user on an empty body. (I5.1, #182: the dual-nav `activePage`
  // model that the old "/ bounces to a compare URL" bridge relied on is gone.)

  it("bare / lands on the corpus contact sheet", async () => {
    useAppState.setState({ activeExperimentId: undefined });
    renderRoutes("/");
    expect(await screen.findByTestId("samples-page")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("bare / lands on the corpus even with an active experiment set", async () => {
    useAppState.setState({ activeExperimentId: 7 });
    renderRoutes("/");
    expect(await screen.findByTestId("samples-page")).toBeInTheDocument();
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
