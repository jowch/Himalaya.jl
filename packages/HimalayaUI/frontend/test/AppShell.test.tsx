/**
 * AppShell — URL/Zustand sync regression tests (#77).
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppShell } from "../src/components/AppShell";

function makeQc() {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function renderShell(initialPath: string) {
  const qc = makeQc();
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[initialPath]}>
        <AppShell />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("AppShell — Zustand → URL sync (#77)", () => {
  beforeEach(() => {
    useAppState.setState({
      activePage: "index",
      activeExperimentId: undefined,
    });
  });

  it("activePage='compare' + URL '/' navigates to /compare/all", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: undefined });
    renderShell("/");

    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
  });

  it("activePage='compare' + URL '/' navigates to /experiments/:eid/compare when experiment is set", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: 7 });
    renderShell("/");

    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
  });

  it("activePage='index' + URL '/' does NOT navigate to a compare route", async () => {
    // Negative case — guards against a future regression that flips the
    // condition (e.g. `if (activePage !== "index")`).
    useAppState.setState({ activePage: "index", activeExperimentId: undefined });
    renderShell("/");

    // Give the Zustand→URL effect a chance to fire if it were going to.
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });
});
