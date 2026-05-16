/**
 * Compare viewing-stale detection — Compare UX C-15 Step 4b.
 *
 * `useStaleAgainstHash` tracks the comparison `content_hash` the user last
 * "saw". While a foreign update drifts it (and no draft is active), the review
 * header surfaces the server-update banner (`compare-status-server-update`).
 *
 * Three behavioral pins:
 *   (1) flips to viewing-stale on content_hash drift while no draft is active;
 *   (2) banner clears on Acknowledge (rebases the acked hash);
 *   (3) does NOT flash during own-op save completion (draft non-null → null
 *       rebases the acked hash so the save's own hash bump is not a "drift").
 */
import { describe, it, expect, beforeEach, vi, afterEach } from "vitest";
import { render, screen, waitFor, act } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { Compare } from "../src/pages/Compare";
import { useAppState } from "../src/state";
import { queryKeys } from "../src/queries";
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
    description: null,
    content_hash: "h1",
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
        snapshot: { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h" },
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
          <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("Compare viewing-stale — Compare UX C-15", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    sessionStorage.clear();
    localStorage.clear();
    useAppState.setState({ activeDraft: null, username: undefined });
    vi.stubGlobal("ResizeObserver", vi.fn(() => ({
      observe: vi.fn(), disconnect: vi.fn(),
    })));
  });
  afterEach(() => {
    vi.unstubAllGlobals();
  });

  it("flips to viewing-stale when content_hash drifts while no draft is active", async () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparison(42), makeComparison({ content_hash: "h1" }));
    renderReview(qc);

    // First sighting anchors the acked hash — no banner on cold load.
    await screen.findByTestId("compare-review-plot");
    expect(screen.queryByTestId("compare-status-server-update")).toBeNull();

    // A foreign SSE update drifts the cached comparison's content_hash.
    act(() => {
      qc.setQueryData(queryKeys.comparison(42), makeComparison({ content_hash: "h2" }));
    });

    // The server-update banner surfaces.
    await waitFor(() => {
      expect(screen.getByTestId("compare-status-server-update")).toBeInTheDocument();
    });
  });

  it("banner clears when Acknowledge is clicked", async () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparison(42), makeComparison({ content_hash: "h1" }));
    renderReview(qc);
    await screen.findByTestId("compare-review-plot");

    act(() => {
      qc.setQueryData(queryKeys.comparison(42), makeComparison({ content_hash: "h2" }));
    });
    const ack = await screen.findByTestId("compare-status-acknowledge");

    act(() => { ack.click(); });

    // Acknowledge rebases the acked hash to the current value → banner clears.
    await waitFor(() => {
      expect(screen.queryByTestId("compare-status-server-update")).toBeNull();
    });
  });

  it("does NOT flash the banner during own-op save completion", async () => {
    const qc = makeQc();
    qc.setQueryData(queryKeys.comparison(42), makeComparison({ content_hash: "h1" }));
    renderReview(qc);
    await screen.findByTestId("compare-review-plot");

    // Simulate an own-op edit session: a draft becomes active, then the save
    // lands — the comparison's content_hash bumps AND the draft clears in the
    // same beat. The draft non-null → null transition must rebase the acked
    // hash so the save's own bump is not surfaced as a foreign drift.
    act(() => {
      useAppState.setState({
        activeDraft: {
          id: 42, baseHash: "h1", title: "T", description: "",
          members: [],
          forkedFromId: undefined, forkedAtHash: undefined,
          viewGroupingMode: undefined, viewShowPeakTicks: undefined,
          viewShowPeakLabels: undefined,
        },
      });
    });
    // Save lands: content_hash bumps and the draft is discarded.
    act(() => {
      qc.setQueryData(queryKeys.comparison(42), makeComparison({ content_hash: "h2" }));
      useAppState.setState({ activeDraft: null });
    });

    // No banner — the own-op completion rebased acked to h2.
    await waitFor(() => {
      expect(screen.getByTestId("compare-review-plot")).toBeInTheDocument();
    });
    expect(screen.queryByTestId("compare-status-server-update")).toBeNull();
  });
});
