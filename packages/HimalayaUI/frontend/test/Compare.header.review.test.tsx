/**
 * Compare header (review mode) — Compare UX C-12.
 *
 * Verifies that ComparePage's review-mode header is now built from the new
 * leaf components (`CompareTitleStrip`, `CompareStatusSurface`,
 * `CompareToolbar`) introduced in issue #139, and that the legacy callsites
 * (`LineageBadge`, `NeedsReviewBadge`, `ForksPopover`, `EditOrForkButton`)
 * have been removed from that header.
 *
 * The component files themselves are NOT deleted here (C-16) — only their
 * ComparePage callsites; so the "absence" assertions target the testids
 * those legacy components emit.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { ComparePage } from "../src/pages/ComparePage";
import { useAppState } from "../src/state";
import type { Comparison } from "../src/api";

// MultiTracePlot pulls Observable Plot — stub it to avoid jsdom weirdness.
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

function makeQc(): QueryClient {
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
    title: "Cubic vs Hex",
    description: "Some description",
    content_hash: "sha256:parent",
    created_by: 7,
    created_at: "2026-01-01",
    updated_at: "2026-01-01",
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    view_grouping_mode: null,
    view_show_peak_ticks: null,
    view_show_peak_labels: null,
    last_event_at: "2026-05-07T00:00:00",
    members: [
      {
        id: 100, comparison_id: 42, exposure_id: null,
        display_order: 0, band_height: 1, y_offset: 0,
        normalization: "qwindow",
        color_override: null, label_override: null,
        q_window_min: null, q_window_max: null, peak_display: null,
        snapshot: {
          effective_peaks: [],
          confirmed_index: null,
          analysis_inputs_hash: "h",
        },
        is_stale: false, created_by: null, created_at: null,
      },
    ],
    ...over,
  };
}

function renderReview(qc: QueryClient) {
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/experiments/7/compare/42"]}>
        <Routes>
          <Route path="/experiments/:eid/compare/:id" element={<ComparePage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("Compare header (review mode) — Compare UX C-12", () => {
  beforeEach(() => {
    useAppState.setState({ username: undefined });
    useAppState.getState().discardDraft();
  });

  it("renders the new title strip + toolbar (not the legacy badges)", async () => {
    const qc = makeQc();
    qc.setQueryData(["comparison", 42] as const, makeComparison());

    renderReview(qc);

    expect(await screen.findByTestId("compare-title-strip")).toBeInTheDocument();
    expect(await screen.findByTestId("compare-toolbar")).toBeInTheDocument();

    // Legacy callsites are gone from the review header.
    expect(screen.queryByTestId("comparison-lineage")).toBeNull();
    expect(screen.queryByTestId("comparison-forks-trigger")).toBeNull();
    expect(screen.queryByTestId("comparison-needs-review")).toBeNull();
    expect(screen.queryByTestId("comparison-edit")).toBeNull();
    expect(screen.queryByTestId("comparison-fork")).toBeNull();
  });

  it("surfaces the needs-review status when a member is stale", async () => {
    const qc = makeQc();
    qc.setQueryData(
      ["comparison", 42] as const,
      makeComparison({
        members: [
          {
            id: 100, comparison_id: 42, exposure_id: null,
            display_order: 0, band_height: 1, y_offset: 0,
            normalization: "qwindow",
            color_override: null, label_override: null,
            q_window_min: null, q_window_max: null, peak_display: null,
            snapshot: {
              effective_peaks: [],
              confirmed_index: null,
              analysis_inputs_hash: "h",
            },
            is_stale: true, created_by: null, created_at: null,
          },
        ],
      }),
    );

    renderReview(qc);

    expect(await screen.findByTestId("compare-status-surface")).toBeInTheDocument();
    expect(screen.getByTestId("compare-status-resnapshot")).toBeInTheDocument();
  });
});
