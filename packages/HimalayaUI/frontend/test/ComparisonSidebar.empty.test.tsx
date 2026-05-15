/**
 * ComparisonSidebar — empty states (Compare UX Phase F, F-3).
 *
 * Two distinct empty branches: no comparisons at all (offers "+ New
 * comparison") vs. a search that matched nothing (offers "Clear search").
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ComparisonSidebar } from "../src/components/ComparisonSidebar";
import { queryKeys } from "../src/queries";
import { useAppState } from "../src/state";
import type { ComparisonSummary } from "../src/api";

function makeQc(): QueryClient {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function makeRow(over: Partial<ComparisonSummary> = {}): ComparisonSummary {
  return {
    id: 1, title: "Cubic vs Hex", description: null, content_hash: "h",
    created_by: 5, created_at: null, updated_at: "2026-05-14T11:00:00Z",
    forked_from_id: null, forked_at_hash: null,
    view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, last_event_at: "2026-05-14T11:00:00Z",
    author_username: "alice", member_count: 4,
    member_phases: ["Pn3m"], has_stale_members: false,
    ...over,
  };
}

function PathProbe(): JSX.Element {
  return <div data-testid="path-probe">{useLocation().pathname}</div>;
}

function renderSidebar(qc: QueryClient): ReturnType<typeof render> {
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7/compare"]}>
        <Routes>
          {["/experiments/:eid/compare", "/experiments/:eid/compare/new"].map(
            (path) => (
              <Route
                key={path}
                path={path}
                element={
                  <>
                    <ComparisonSidebar
                      experimentId={7}
                      scope="experiment"
                      activeComparisonId={undefined}
                    />
                    <PathProbe />
                  </>
                }
              />
            ),
          )}
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("Sidebar empty states — Compare UX F-3", () => {
  beforeEach(() => {
    vi.spyOn(Date, "now").mockReturnValue(Date.parse("2026-05-14T12:00:00Z"));
    useAppState.setState({ activeDraft: null, username: undefined });
  });

  it("shows the no-comparisons-yet state + a + New button", () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparisons(7), []);
    renderSidebar(qc);
    expect(
      screen.getByText(/No comparisons in this experiment yet/),
    ).toBeInTheDocument();
    expect(screen.getByTestId("sidebar-empty-new")).toBeInTheDocument();
  });

  it("the no-comparisons + New button navigates to the new route", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparisons(7), []);
    renderSidebar(qc);
    await user.click(screen.getByTestId("sidebar-empty-new"));
    expect(screen.getByTestId("path-probe")).toHaveTextContent(
      "/experiments/7/compare/new",
    );
  });

  it("shows the search-empty state + a Clear search button", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparisons(7), [makeRow({ title: "Cubic" })]);
    renderSidebar(qc);
    await user.type(
      screen.getByTestId("comparison-sidebar-search"),
      "zzz-no-match",
    );
    expect(screen.getByText(/No matches/)).toBeInTheDocument();
    expect(screen.getByTestId("sidebar-empty-clear")).toBeInTheDocument();
  });

  it("the Clear search button restores the unfiltered list", async () => {
    const user = userEvent.setup();
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparisons(7), [makeRow({ title: "Cubic" })]);
    renderSidebar(qc);
    await user.type(
      screen.getByTestId("comparison-sidebar-search"),
      "zzz-no-match",
    );
    await user.click(screen.getByTestId("sidebar-empty-clear"));
    expect(screen.queryByTestId("sidebar-empty-clear")).not.toBeInTheDocument();
    expect(screen.getByTestId("comparison-list-item")).toBeInTheDocument();
  });
});
