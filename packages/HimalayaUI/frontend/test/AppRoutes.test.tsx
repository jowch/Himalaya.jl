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
      activePage: "index",
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

  it("mounts AppShell (not the corpus shell) at /index", async () => {
    renderRoutes("/index");
    expect(await screen.findByTestId("app-shell")).toBeInTheDocument();
    expect(screen.queryByTestId("corpus-shell")).toBeNull();
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
          <button data-testid="go-index" onClick={() => navigate("/index")}>index</button>
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

    // Cross to a legacy route — AppShell should mount, CorpusShell should unmount.
    fireEvent.click(screen.getByTestId("go-index"));
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
      activePage: "compare",
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
});

describe("AppRoutes — Zustand → URL compare-sync (#77)", () => {
  beforeEach(() => {
    // Reset the ephemeral URL-resolution fields too — a prior test that
    // parked the store on a stale path must not leak into the next.
    useAppState.setState({
      activePage: "index",
      activeExperimentId: undefined,
      staleUrlContext: null,
      resolving: false,
    });
  });

  it("activePage='compare' + URL '/' navigates to /compare/all", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: undefined });
    renderRoutes("/");
    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
  });

  it("activePage='compare' + URL '/' navigates to /experiments/:eid/compare when an experiment is set", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: 7 });
    renderRoutes("/");
    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
  });

  it("activePage='index' + URL '/' does NOT navigate to a compare route", async () => {
    useAppState.setState({ activePage: "index", activeExperimentId: undefined });
    renderRoutes("/");
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("back-nav from /experiments/:eid/compare/:id to '/' bounces to a compare URL (intentional)", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: 7 });
    renderRoutes("/experiments/7/compare/123", 1);
    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
  });
});
