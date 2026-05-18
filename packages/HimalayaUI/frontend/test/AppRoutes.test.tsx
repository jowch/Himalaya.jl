/**
 * AppRoutes — the single hoisted route table. Tests the nav-bridge: new
 * routes mount the corpus shell, legacy routes mount AppShell, and the two
 * nav models do not desync. Includes the relocated #77 compare-sync tests
 * (formerly AppShell.test.tsx — AppShell is no longer a router).
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppRoutes } from "../src/components/AppRoutes";

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
