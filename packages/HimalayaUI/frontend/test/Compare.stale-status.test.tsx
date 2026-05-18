/**
 * Compare stale-status tests (originally Plan §Phase 9, Task 9.6).
 *
 * Compare UX C-16: the standalone `NeedsReviewBadge` component was deleted
 * once C-12 folded the stale signal into `CompareStatusSurface`. The
 * direct-render describe that exercised the deleted component is gone; what
 * remains here is page-level coverage of the comparison-level stale
 * disjunction (via `Compare` → `CompareStatusSurface`) and the per-member
 * `data-stale` attribute (via `MemberMetaRow`).
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { Compare } from "../src/pages/Compare";
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

// Compare UX C-12: the ComparePage review header no longer mounts
// `NeedsReviewBadge`. The stale signal is now surfaced by
// `CompareStatusSurface` (`compare-status-surface` + a `compare-status-
// resnapshot` button). The new banner is NOT author-gated — anyone viewing a
// stale comparison sees the re-snapshot affordance, so the prior author/
// non-author clickability split collapses. These tests migrate the
// page-level assertions to the new selectors.
describe("Compare review header — stale status surface via ComparePage", () => {
  beforeEach(() => {
    useAppState.setState({ username: undefined });
  });

  it("does not render the status surface when no member is stale", async () => {
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
            <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    await waitFor(() => screen.getByTestId("compare-review-plot"));
    expect(screen.queryByTestId("compare-status-surface")).toBeNull();
    expect(screen.queryByTestId("comparison-needs-review")).toBeNull();
  });

  it("renders the status surface when any member is stale", async () => {
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
            <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    const surface = await waitFor(() =>
      screen.getByTestId("compare-status-surface"),
    );
    expect(surface).toBeInTheDocument();
    expect(screen.getByTestId("compare-status-resnapshot")).toBeInTheDocument();
  });

  it("re-snapshot click navigates to the bare comparison URL", async () => {
    const qc = makeQc();
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
              `/edit` segment, so re-snapshot navigates to the bare URL. */}
          <LocationSpy />
          <Routes>
            <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    const btn = await waitFor(() =>
      screen.getByTestId("compare-status-resnapshot"),
    );
    fireEvent.click(btn);
    await waitFor(() =>
      expect(screen.getByTestId("current-location").textContent)
        .toBe("/experiments/7/compare/42"),
    );
  });
});

describe("MemberMetaRow — per-member data-stale", () => {
  it("MemberMetaRow surfaces data-stale on stale members", async () => {
    const { MemberMetaRow } = await import("../src/components/MemberMetaRow");
    render(
      <MemberMetaRow
        member={{
          id: 1, series_id: 42, exposure_id: null,
          display_order: 0, band_height: 1, y_offset: 0,
          normalization: "qwindow",
          color_override: null, label_override: null,
          q_window_min: null, q_window_max: null, peak_display: null,
          snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
          is_stale: true, created_by: null, created_at: null,
        }}
        top={0} height={50} mode="review" displayLabel="row-label"
        expanded={false} onToggleExpand={() => {}}
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
          id: 1, series_id: 42, exposure_id: null,
          display_order: 0, band_height: 1, y_offset: 0,
          normalization: "qwindow",
          color_override: null, label_override: null,
          q_window_min: null, q_window_max: null, peak_display: null,
          snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
          is_stale: false, created_by: null, created_at: null,
        }}
        top={0} height={50} mode="review" displayLabel="row-label"
        expanded={false} onToggleExpand={() => {}}
      />,
    );
    expect(screen.getByTestId("member-meta-row")).not.toHaveAttribute("data-stale");
  });
});
