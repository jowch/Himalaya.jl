/**
 * ComparisonSidebar — draft indicator (Compare UX Phase F, F-2).
 *
 * When an `activeDraft` is open for a comparison, its sidebar row gets a •
 * dot prefix, a "(draft)" title suffix, and a "by you · just now" byline.
 */
import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { ComparisonSidebar } from "../src/components/ComparisonSidebar";
import { queryKeys } from "../src/queries";
import { useAppState } from "../src/state";
import type { ComparisonSummary } from "../src/api";
import type { ActiveDraft } from "../src/lib/comparison/draft";

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
    member_phases: ["Pn3m", "Hex", "Lam"], has_stale_members: false,
    ...over,
  };
}

function makeDraft(id: number | undefined): ActiveDraft {
  return {
    id, baseHash: undefined, title: "draft title", description: "",
    members: [], forkedFromId: undefined, forkedAtHash: undefined,
  };
}

function renderSidebar(qc: QueryClient): ReturnType<typeof render> {
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7/compare"]}>
        <Routes>
          <Route
            path="/experiments/:eid/compare"
            element={
              <ComparisonSidebar
                experimentId={7}
                scope="experiment"
                activeComparisonId={undefined}
              />
            }
          />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("Sidebar draft indicator — Compare UX F-2", () => {
  beforeEach(() => {
    vi.spyOn(Date, "now").mockReturnValue(Date.parse("2026-05-14T12:00:00Z"));
    useAppState.setState({ activeDraft: null, username: undefined });
  });
  afterEach(() => {
    vi.restoreAllMocks();
    useAppState.setState({ activeDraft: null });
  });

  it("renders • + (draft) suffix on the row matching the active draft id", () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparisons(7), [
      makeRow({ id: 5, title: "Editing this" }),
      makeRow({ id: 6, title: "Untouched" }),
    ]);
    useAppState.setState({ activeDraft: makeDraft(5) });
    renderSidebar(qc);

    const rows = screen.getAllByTestId("comparison-list-item");
    const draftRow = rows.find((r) => r.dataset.comparisonId === "5")!;
    const otherRow = rows.find((r) => r.dataset.comparisonId === "6")!;

    expect(draftRow).toHaveTextContent("Editing this (draft)");
    expect(draftRow.querySelector('[data-testid="sidebar-draft-dot"]'))
      .not.toBeNull();
    expect(otherRow).not.toHaveTextContent("(draft)");
    expect(otherRow.querySelector('[data-testid="sidebar-draft-dot"]'))
      .toBeNull();
  });

  it("shows a 'by you · just now' byline on the draft row", () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparisons(7), [
      makeRow({ id: 5, author_username: "someone-else" }),
    ]);
    useAppState.setState({ activeDraft: makeDraft(5) });
    renderSidebar(qc);
    expect(screen.getByText(/by you · just now/)).toBeInTheDocument();
  });

  it("renders no draft dot when there is no active draft", () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparisons(7), [makeRow({ id: 5 })]);
    renderSidebar(qc);
    expect(screen.queryByTestId("sidebar-draft-dot")).not.toBeInTheDocument();
  });

  it("renders no draft dot when the active draft is unsaved (id undefined)", () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparisons(7), [makeRow({ id: 5 })]);
    useAppState.setState({ activeDraft: makeDraft(undefined) });
    renderSidebar(qc);
    expect(screen.queryByTestId("sidebar-draft-dot")).not.toBeInTheDocument();
  });
});
