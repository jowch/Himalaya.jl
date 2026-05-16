/**
 * ComparisonSidebar — listing-projection row layout (Compare UX Phase F, F-1).
 *
 * Asserts the redesigned sidebar row consumes the Phase A projection:
 * phase summary (`member_phases` · `member_count` traces), author byline
 * ("by you" / "by <name>"), relative age from `last_event_at`, and a stale
 * ⚠ marker. Also pins the client-side sort to `last_event_at` desc.
 */
import { describe, it, expect, beforeEach, afterEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
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
    id: 1,
    title: "Cubic vs Hex",
    description: null,
    content_hash: "h",
    created_by: 5,
    created_at: null,
    updated_at: "2026-05-14T11:00:00Z",
    forked_from_id: null,
    forked_at_hash: null,
    view_grouping_mode: null,
    view_show_peak_ticks: null,
    view_show_peak_labels: null,
    last_event_at: "2026-05-14T11:00:00Z",
    author_username: "alice",
    member_count: 4,
    member_phases: ["Pn3m", "Hex", "Lam"],
    member_phase_count: 3,
    has_stale_members: false,
    ...over,
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

describe("Sidebar row projection — Compare UX F-1", () => {
  beforeEach(() => {
    // Freeze "now" 1h after the fixture's last_event_at so relative age is
    // deterministic ("1h ago").
    vi.spyOn(Date, "now").mockReturnValue(
      Date.parse("2026-05-14T12:00:00Z"),
    );
    useAppState.setState({ username: undefined });
  });
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("renders phases · traces, author byline, relative age", () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparisons(7), [makeRow()]);
    renderSidebar(qc);
    expect(
      screen.getByText(/Pn3m · Hex · Lam · 4 traces/),
    ).toBeInTheDocument();
    expect(screen.getByText(/by alice/)).toBeInTheDocument();
    expect(screen.getByText(/edited 1h ago/)).toBeInTheDocument();
  });

  it("renders 'by you' when current user is the author", () => {
    const qc = makeQc();
    // useCurrentUserId resolves the Zustand username against the cached
    // ["users"] listing → id 5 === row.created_by → "by you".
    useAppState.setState({ username: "alice" });
    qc.setQueryData(["users"], [{ id: 5, username: "alice" }]);
    qc.setQueryData(queryKeys.comparisons(7), [makeRow({ created_by: 5 })]);
    renderSidebar(qc);
    expect(screen.getByText(/by you/)).toBeInTheDocument();
  });

  it("renders ⚠ when has_stale_members is true", () => {
    const qc = makeQc();
    qc.setQueryData(
      queryKeys.comparisons(7),
      [makeRow({ has_stale_members: true })],
    );
    renderSidebar(qc);
    expect(screen.getByTestId("sidebar-stale-warn")).toBeInTheDocument();
  });

  it("does not render ⚠ when has_stale_members is false", () => {
    const qc = makeQc();
    qc.setQueryData(
      queryKeys.comparisons(7),
      [makeRow({ has_stale_members: false })],
    );
    renderSidebar(qc);
    expect(screen.queryByTestId("sidebar-stale-warn")).not.toBeInTheDocument();
  });

  it("renders '+N more' when member_phase_count exceeds the shown phases", () => {
    const qc = makeQc();
    // The backend caps `member_phases` at 3; `member_phase_count` carries
    // the true distinct-phase total so the client can render the overflow.
    qc.setQueryData(
      queryKeys.comparisons(7),
      [makeRow({
        member_phases: ["Pn3m", "Im3m", "Ia3d"],
        member_phase_count: 5,
        member_count: 5,
      })],
    );
    renderSidebar(qc);
    expect(
      screen.getByText(/Pn3m · Im3m · Ia3d · \+2 more · 5 traces/),
    ).toBeInTheDocument();
  });

  it("falls back to 'N traces' when there are no phases", () => {
    const qc = makeQc();
    qc.setQueryData(
      queryKeys.comparisons(7),
      [makeRow({ member_phases: [], member_phase_count: 0, member_count: 2 })],
    );
    renderSidebar(qc);
    expect(screen.getByText(/^2 traces$/)).toBeInTheDocument();
  });
});

describe("Sidebar sort key — Compare UX F-1 / spec §8.4", () => {
  beforeEach(() => {
    vi.spyOn(Date, "now").mockReturnValue(Date.parse("2026-05-15T00:00:00Z"));
    useAppState.setState({ username: undefined });
  });
  afterEach(() => {
    vi.restoreAllMocks();
  });

  it("sorts pinned-first, then by last_event_at desc (not updated_at)", () => {
    const qc = makeQc();
    // A pinned: last_event_at OLD but updated_at OLDEST.
    // B unpinned: last_event_at NEWEST, updated_at MID.
    // C unpinned: last_event_at MID,    updated_at NEWEST.
    // last_event_at order → B before C. updated_at order would be C before B.
    const rowA = makeRow({
      id: 1, title: "A",
      last_event_at: "2026-05-14T10:00:00Z", updated_at: "2026-05-01T00:00:00Z",
    });
    const rowB = makeRow({
      id: 2, title: "B",
      last_event_at: "2026-05-14T20:00:00Z", updated_at: "2026-05-10T00:00:00Z",
    });
    const rowC = makeRow({
      id: 3, title: "C",
      last_event_at: "2026-05-14T15:00:00Z", updated_at: "2026-05-12T00:00:00Z",
    });
    qc.setQueryData(queryKeys.comparisons(7), [rowA, rowB, rowC]);
    qc.setQueryData(queryKeys.comparisonPins, [1]); // pin A only
    renderSidebar(qc);
    const items = screen.getAllByTestId("comparison-list-item");
    expect(items[0]).toHaveTextContent("A"); // pinned first
    expect(items[1]).toHaveTextContent("B"); // last_event_at newest
    expect(items[2]).toHaveTextContent("C");
  });

  it("orders mixed-format last_event_at correctly (space-sep vs T-sep)", () => {
    const qc = makeQc();
    // `last_event_at` mixes string formats: `MAX(user_actions.timestamp)`
    // is space-separated with no `Z`; the `c.updated_at` COALESCE fallback
    // is `T`-separated with `Z` (see comparisons.jl). A naive `localeCompare`
    // misorders them — the space-sep row here is the more recent one.
    const rowOld = makeRow({
      id: 1, title: "OldTsep", last_event_at: "2026-05-14T00:01:00Z",
    });
    const rowNew = makeRow({
      id: 2, title: "NewSpace", last_event_at: "2026-05-14 23:59:00",
    });
    qc.setQueryData(queryKeys.comparisons(7), [rowOld, rowNew]);
    renderSidebar(qc);
    const items = screen.getAllByTestId("comparison-list-item");
    expect(items[0]).toHaveTextContent("NewSpace");
    expect(items[1]).toHaveTextContent("OldTsep");
  });

  it("sorts rows with null last_event_at last", () => {
    const qc = makeQc();
    const rowTimed = makeRow({
      id: 1, title: "Timed", last_event_at: "2026-05-14T10:00:00Z",
    });
    const rowNull = makeRow({ id: 2, title: "Null", last_event_at: null });
    qc.setQueryData(queryKeys.comparisons(7), [rowNull, rowTimed]);
    renderSidebar(qc);
    const items = screen.getAllByTestId("comparison-list-item");
    expect(items[0]).toHaveTextContent("Timed");
    expect(items[1]).toHaveTextContent("Null");
  });
});
