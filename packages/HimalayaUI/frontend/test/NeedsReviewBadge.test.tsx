/**
 * NeedsReviewBadge tests (Plan §Phase 9, Task 9.6).
 *
 * - Author can click → navigates to the bare comparison URL.
 * - Non-author sees badge as informational; click is a no-op.
 * - Page-level mounting: visible when any member stale, hidden otherwise.
 * - Per-member `data-stale` attribute already covered by MemberMetaRow.test;
 *   here we verify the comparison-level disjunction at the ComparePage level.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { NeedsReviewBadge } from "../src/components/NeedsReviewBadge";
import { ComparePage } from "../src/pages/ComparePage";
import { useAppState } from "../src/state";
import type { Comparison } from "../src/api";

vi.mock("@observablehq/plot", () => ({
  plot: vi.fn(() => {
    const el = document.createElement("div");
    (el as unknown as { scale: (n: string) => unknown }).scale = (n) =>
      n === "x"
        ? { invert: (px: number) => px / 100, apply: (q: number) => q * 100 }
        : undefined;
    return el;
  }),
  line: vi.fn(() => ({})),
  dot:  vi.fn(() => ({})),
  text: vi.fn(() => ({})),
  link: vi.fn(() => ({})),
}));

function makeQc() {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function makeComparison(over: Partial<Comparison> = {}): Comparison {
  return {
    id: 42,
    title: "T",
    description: null,
    content_hash: "sha256:x",
    created_by: 7,
    created_at: "2026-01-01",
    updated_at: "2026-01-01",
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    members: [],
    ...over,
  };
}

describe("NeedsReviewBadge — direct render", () => {
  beforeEach(() => {
    useAppState.setState({ username: undefined });
  });

  it("non-author renders informational badge with data-clickable=false", () => {
    render(
      <MemoryRouter>
        <QueryClientProvider client={makeQc()}>
          <NeedsReviewBadge comparisonId={42} experimentId={7} authorUserId={7} />
        </QueryClientProvider>
      </MemoryRouter>,
    );
    const badge = screen.getByTestId("comparison-needs-review");
    expect(badge).toBeInTheDocument();
    expect(badge).toHaveAttribute("data-clickable", "false");
    expect(badge.textContent).toContain("Needs Review");
  });
});

describe("NeedsReviewBadge — author clickability via ComparePage", () => {
  beforeEach(() => {
    useAppState.setState({ username: undefined });
  });

  it("does not render when no member is stale", async () => {
    const qc = makeQc();
    qc.setQueryData(
      ["comparison", 42] as const,
      makeComparison({
        members: [
          {
            id: 1, comparison_id: 42, exposure_id: null,
            display_order: 0, band_height: 1, y_offset: 0,
            normalization: "qwindow",
            color_override: null, label_override: null,
            q_window_min: null, q_window_max: null, peak_display: null,
            snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
            is_stale: false, created_by: null, created_at: null,
          },
        ],
      }),
    );
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments/7/compare/42"]}>
          <Routes>
            <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    await waitFor(() => screen.getByTestId("compare-review-plot"));
    expect(screen.queryByTestId("comparison-needs-review")).toBeNull();
  });

  it("renders when any member is stale", async () => {
    const qc = makeQc();
    qc.setQueryData(
      ["comparison", 42] as const,
      makeComparison({
        members: [
          {
            id: 1, comparison_id: 42, exposure_id: null,
            display_order: 0, band_height: 1, y_offset: 0,
            normalization: "qwindow",
            color_override: null, label_override: null,
            q_window_min: null, q_window_max: null, peak_display: null,
            snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
            is_stale: false, created_by: null, created_at: null,
          },
          {
            id: 2, comparison_id: 42, exposure_id: null,
            display_order: 1, band_height: 1, y_offset: 0,
            normalization: "qwindow",
            color_override: null, label_override: null,
            q_window_min: null, q_window_max: null, peak_display: null,
            snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
            // STALE
            is_stale: true, created_by: null, created_at: null,
          },
        ],
      }),
    );
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments/7/compare/42"]}>
          <Routes>
            <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    await waitFor(() => screen.getByTestId("comparison-needs-review"));
    const badge = screen.getByTestId("comparison-needs-review");
    expect(badge).toBeInTheDocument();
  });

  it("author (created_by matches current user) → badge is clickable; navigates to the comparison", async () => {
    const qc = makeQc();
    // Pre-cache the users list so `useCurrentUserId` resolves synchronously.
    qc.setQueryData(["users"] as const, [
      { id: 7, username: "alice", first_name: null, last_name: null },
    ]);
    qc.setQueryData(
      ["comparison", 42] as const,
      makeComparison({
        created_by: 7,
        members: [
          {
            id: 1, comparison_id: 42, exposure_id: null,
            display_order: 0, band_height: 1, y_offset: 0,
            normalization: "qwindow",
            color_override: null, label_override: null,
            q_window_min: null, q_window_max: null, peak_display: null,
            snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
            is_stale: true, created_by: null, created_at: null,
          },
        ],
      }),
    );
    useAppState.setState({ username: "alice" });

    function LocationSpy() {
      const loc = useLocation();
      return <div data-testid="current-location">{loc.pathname}</div>;
    }

    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments/7/compare/42"]}>
          {/* Always-mounted so navigation is observable even when it lands
              on the route ComparePage owns — Compare UX Phase B dropped the
              `/edit` segment, so the author badge navigates to the bare URL. */}
          <LocationSpy />
          <Routes>
            <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    const badge = await waitFor(() => screen.getByTestId("comparison-needs-review"));
    await waitFor(() => expect(badge).toHaveAttribute("data-clickable", "true"));
    fireEvent.click(badge);
    await waitFor(() =>
      expect(screen.getByTestId("current-location").textContent)
        .toBe("/experiments/7/compare/42"),
    );
  });

  it("non-author (created_by mismatches) → badge not clickable; click is a no-op", async () => {
    const qc = makeQc();
    qc.setQueryData(["users"] as const, [
      { id: 99, username: "bob", first_name: null, last_name: null },
    ]);
    qc.setQueryData(
      ["comparison", 42] as const,
      makeComparison({
        // Comparison created by user 7 (alice)
        created_by: 7,
        members: [
          {
            id: 1, comparison_id: 42, exposure_id: null,
            display_order: 0, band_height: 1, y_offset: 0,
            normalization: "qwindow",
            color_override: null, label_override: null,
            q_window_min: null, q_window_max: null, peak_display: null,
            snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
            is_stale: true, created_by: null, created_at: null,
          },
        ],
      }),
    );
    // ...but the current user is bob (id 99).
    useAppState.setState({ username: "bob" });

    function LocationSpy() {
      const loc = useLocation();
      return <div data-testid="current-location">{loc.pathname}</div>;
    }

    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments/7/compare/42"]}>
          <LocationSpy />
          <Routes>
            <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    const badge = await waitFor(() => screen.getByTestId("comparison-needs-review"));
    expect(badge).toHaveAttribute("data-clickable", "false");
    fireEvent.click(badge);
    // Non-author click is a no-op — still on the bare comparison page.
    expect(screen.getByTestId("current-location").textContent)
      .toBe("/experiments/7/compare/42");
  });
});

describe("NeedsReviewBadge — per-member data-stale", () => {
  it("MemberMetaRow surfaces data-stale on stale members", async () => {
    const { MemberMetaRow } = await import("../src/components/MemberMetaRow");
    render(
      <MemberMetaRow
        member={{
          id: 1, comparison_id: 42, exposure_id: null,
          display_order: 0, band_height: 1, y_offset: 0,
          normalization: "qwindow",
          color_override: null, label_override: null,
          q_window_min: null, q_window_max: null, peak_display: null,
          snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
          is_stale: true, created_by: null, created_at: null,
        }}
        top={0} height={50} mode="review" displayLabel="row-label"
      />,
    );
    const row = screen.getByTestId("member-meta-row");
    expect(row).toHaveAttribute("data-stale");
    expect(screen.getByTestId("member-meta-stale-icon")).toBeInTheDocument();
  });

  it("MemberMetaRow data-stale absent when not stale", async () => {
    const { MemberMetaRow } = await import("../src/components/MemberMetaRow");
    render(
      <MemberMetaRow
        member={{
          id: 1, comparison_id: 42, exposure_id: null,
          display_order: 0, band_height: 1, y_offset: 0,
          normalization: "qwindow",
          color_override: null, label_override: null,
          q_window_min: null, q_window_max: null, peak_display: null,
          snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
          is_stale: false, created_by: null, created_at: null,
        }}
        top={0} height={50} mode="review" displayLabel="row-label"
      />,
    );
    expect(screen.getByTestId("member-meta-row")).not.toHaveAttribute("data-stale");
  });
});
