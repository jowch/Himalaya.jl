/**
 * AppRoutes — the single hoisted route table. Tests the nav-bridge: new
 * routes mount the corpus shell, legacy routes mount AppShell, and the two
 * nav models do not desync. Includes the relocated #77 compare-sync tests
 * (formerly AppShell.test.tsx — AppShell is no longer a router).
 */
import { describe, it, expect, beforeEach, vi, afterEach } from "vitest";
import { render, screen, waitFor, fireEvent } from "@testing-library/react";
import { MemoryRouter, useNavigate } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppRoutes } from "../src/components/AppRoutes";
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

describe("AppRoutes — nav-bridge shell selection", () => {
  beforeEach(() => {
    // Reset the ephemeral URL-resolution fields too — a prior test that
    // parked the store on a stale path must not leak into the next.
    useAppState.setState({
      activePage: "none",
      activeExperimentId: undefined,
      staleUrlContext: null,
      resolving: false,
    });
  });

  it("mounts the corpus shell (not AppShell) at /samples", async () => {
    renderRoutes("/samples");
    expect(await screen.findByTestId("corpus-shell")).toBeInTheDocument();
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();
  });

  it("mounts AppShell (not the corpus shell) at the stale catch-all", async () => {
    // I3.6 (#177): Compare is retired, so the only AppShell route left is the
    // `*` catch-all (stale/unknown paths). An unknown path mounts AppShell and
    // renders the StaleUrlPage body; corpus surfaces stay on CorpusShell.
    renderRoutes("/totally/unknown");
    expect(await screen.findByTestId("app-shell")).toBeInTheDocument();
    expect(screen.queryByTestId("corpus-shell")).toBeNull();
  });

  it("redirects a /compare* URL to the series folio (Compare retired, #177)", async () => {
    renderRoutes("/compare/all");
    expect(await screen.findByTestId("series-folio-page")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("redirects an experiment-scoped /compare URL to the series folio (#177)", async () => {
    renderRoutes("/experiments/7/compare/123");
    expect(await screen.findByTestId("series-folio-page")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("does not flag /samples as a stale path", async () => {
    // The legacy URL-sync hooks live in AppShell, which is not mounted on a
    // corpus route — so /samples cannot be parsed as a stale path.
    renderRoutes("/samples");
    await screen.findByTestId("corpus-shell");
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });

  it("remounts the correct shell when navigating across the corpus/legacy boundary", async () => {
    // Headline structural invariant: AppShell unmounts when you cross to a
    // corpus route, and remounts when you cross back to a legacy route.
    function NavButtons(): JSX.Element {
      const navigate = useNavigate();
      return (
        <>
          <button data-testid="go-samples" onClick={() => navigate("/samples")}>samples</button>
          {/* I3.6 (#177): Compare is retired (redirects out), so the only
              surviving AppShell surface is the `*` catch-all (stale/unknown
              paths). Cross to a stale path to exercise the AppShell boundary. */}
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

    // Start: corpus shell is up, legacy shell is absent.
    expect(await screen.findByTestId("corpus-shell")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();

    // Cross to the legacy catch-all — AppShell mounts, CorpusShell unmounts.
    fireEvent.click(screen.getByTestId("go-stale"));
    await waitFor(() => {
      expect(screen.getByTestId("app-shell")).toBeInTheDocument();
      expect(screen.queryByTestId("corpus-shell")).toBeNull();
    });

    // Cross back to corpus route — CorpusShell should remount, AppShell should unmount.
    fireEvent.click(screen.getByTestId("go-samples"));
    await waitFor(() => {
      expect(screen.getByTestId("corpus-shell")).toBeInTheDocument();
      expect(screen.queryByTestId("app-shell")).toBeNull();
    });
  });
});

describe("AppRoutes — I4.4 index cutover redirects", () => {
  beforeEach(() => {
    useAppState.setState({
      activePage: "none",
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
    expect(await screen.findByTestId("focus-workspace-page")).toBeInTheDocument();
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
    // The Compare nav-bridge that used to bounce typo'd URLs is gone with the
    // Compare page (#177); a stale path now cleanly renders "Page not found"
    // via the AppShell catch-all PageBody.
    renderRoutes("/foo/bar");
    expect(await screen.findByTestId("stale-url-page")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });
});

describe("AppRoutes — bare / always lands on the corpus (#77 / I4.4)", () => {
  beforeEach(() => {
    useAppState.setState({
      activePage: "none",
      activeExperimentId: undefined,
      staleUrlContext: null,
      resolving: false,
    });
  });

  // I4.4 (#181): the #77 "empty PageBody at /" risk is eliminated differently
  // now. Bare `/` is a standalone redirect to /samples (outside AppShell), so
  // a cold `/` can never strand the user on an empty body — regardless of the
  // persisted `activePage`. The old "activePage='compare' + / bounces to a
  // compare URL" bridge is retired with the Index surface.

  it("bare / lands on the corpus contact sheet even when activePage='compare'", async () => {
    useAppState.setState({ activePage: "none", activeExperimentId: undefined });
    renderRoutes("/");
    expect(await screen.findByTestId("samples-page")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("bare / lands on the corpus even with an active experiment set", async () => {
    useAppState.setState({ activePage: "none", activeExperimentId: 7 });
    renderRoutes("/");
    expect(await screen.findByTestId("samples-page")).toBeInTheDocument();
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("a compare URL redirects to the series folio (Compare retired, #177)", async () => {
    useAppState.setState({ activePage: "none", activeExperimentId: 7 });
    renderRoutes("/experiments/7/compare/123");
    await waitFor(() => {
      expect(screen.getByTestId("series-folio-page")).toBeInTheDocument();
    });
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });
});
