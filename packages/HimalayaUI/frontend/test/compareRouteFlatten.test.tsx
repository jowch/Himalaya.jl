/**
 * Compare /edit redirect — Compare UX Phase B-1.
 *
 * The `/edit` URL segment is being dropped: `/compare/:id` is the only
 * route shape. Old `/edit` deep-links must still resolve, so AppShell
 * registers a redirect Route that 301-equivalents them to the bare URL.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
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

function LocationProbe(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="probe">{loc.pathname}{loc.search}</div>;
}

function renderShell(initialPath: string) {
  const qc = makeQc();
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[initialPath]}>
        <Routes>
          <Route
            path="*"
            element={
              <>
                <LocationProbe />
                <AppShell />
              </>
            }
          />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("Compare /edit redirect — Compare UX B-1", () => {
  beforeEach(() => {
    useAppState.setState({ activePage: "compare", activeExperimentId: undefined });
  });

  it("redirects /experiments/:eid/compare/:id/edit → /experiments/:eid/compare/:id", async () => {
    const { findByTestId } = renderShell("/experiments/10/compare/5/edit");
    const probe = await findByTestId("probe");
    expect(probe.textContent).toBe("/experiments/10/compare/5");
  });

  it("redirects /compare/all/:id/edit → /compare/all/:id", async () => {
    const { findByTestId } = renderShell("/compare/all/7/edit");
    const probe = await findByTestId("probe");
    expect(probe.textContent).toBe("/compare/all/7");
  });

  it("preserves the query string through the redirect", async () => {
    const { findByTestId } = renderShell("/experiments/10/compare/5/edit?tab=peaks");
    const probe = await findByTestId("probe");
    expect(probe.textContent).toBe("/experiments/10/compare/5?tab=peaks");
  });
});
