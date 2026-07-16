// test/AppRoutes.samples-redirect.test.tsx
// T3.2: /samples must redirect to /experiments
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/print/shell/AppRoutes";
import * as api from "../src/api";

function makeQc() {
  return new QueryClient({ defaultOptions: { queries: { retry: false } } });
}

function renderApp({ route }: { route: string }) {
  render(
    <QueryClientProvider client={makeQc()}>
      <MemoryRouter initialEntries={[route]}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("AppRoutes T3.2: /samples redirect", () => {
  beforeEach(() => {
    vi.spyOn(api, "listExperiments").mockResolvedValue([]);
  });
  afterEach(() => vi.restoreAllMocks());

  it("/samples redirects to /experiments", async () => {
    renderApp({ route: "/samples" });
    // After redirect, the ExperimentsHomePage renders "Your experiments"
    await waitFor(() =>
      expect(screen.getByText("All experiments")).toBeInTheDocument(),
    );
    // No SamplesPage (it's been removed from /samples)
    expect(screen.queryByTestId("samples-page")).toBeNull();
  });
});
