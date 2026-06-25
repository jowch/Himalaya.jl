// test/AppRoutes.loupe-route.test.tsx
// Critical regression: /sample/:id/loupe must match the registered Loupe route
// (NOT fall through to the `*` catch-all → StaleUrlPage).
// Unit tests that define their own isolated route tables can't catch this;
// only a test that renders the real AppRoutes can.
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/print/shell/AppRoutes";
import * as api from "../src/api";

function makeQc() {
  return new QueryClient({ defaultOptions: { queries: { retry: false } } });
}

function renderApp(route: string) {
  render(
    <QueryClientProvider client={makeQc()}>
      <MemoryRouter initialEntries={[route]}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("AppRoutes — flat Loupe route /sample/:id/loupe", () => {
  beforeEach(() => {
    // LoupePage calls useCorpusSamples() → listCorpusSamples and
    // useExposures() → listExposures. Return empty arrays so the page renders
    // the loupe-not-found body (sample not in corpus) rather than loading forever.
    vi.spyOn(api, "listCorpusSamples").mockResolvedValue([]);
    vi.spyOn(api, "listExposures").mockResolvedValue([]);
  });
  afterEach(() => vi.restoreAllMocks());

  it("renders LoupePage (not StaleUrlPage) at /sample/42/loupe", async () => {
    renderApp("/sample/42/loupe");

    // LoupePage must render — either its body or its not-found branch.
    // Both carry data-testid="loupe-page". That they appear proves the flat
    // route matched and StaleUrlPage was NOT rendered.
    await waitFor(() =>
      expect(screen.getByTestId("loupe-page")).toBeInTheDocument(),
    );
    // Explicitly assert the catch-all did NOT fire.
    expect(screen.queryByTestId("stale-url-page")).toBeNull();
  });

  it("/sample/42/loupe does NOT render StaleUrlPage", async () => {
    renderApp("/sample/42/loupe");

    // Wait long enough for any redirect/stale classification to settle.
    await waitFor(() =>
      expect(screen.queryByTestId("stale-url-page")).toBeNull(),
    );
  });
});
