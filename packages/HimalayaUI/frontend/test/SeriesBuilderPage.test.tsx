import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import type { Series, SeriesMember } from "../src/api";
import { SeriesBuilderPage } from "../src/pages/SeriesBuilderPage";

const h = vi.hoisted(() => ({
  seriesQ: {} as { data?: Series; isLoading: boolean; isError: boolean },
}));
vi.mock("../src/queries", () => ({
  useSeries: () => h.seriesQ,
  useMemberTraces: () => new Map(),
  useMemberTracesLoading: () => false,
  useMemberExposures: () => new Map(),
  useMemberSamples: () => new Map(),
}));
// MultiTracePlot touches Observable Plot / ResizeObserver; stub it — the
// page's read/state behavior does not depend on its internal render.
vi.mock("../src/components/MultiTracePlot", () => ({
  MultiTracePlot: () => <div data-testid="mock-multi-trace-plot" />,
  COMPARE_PLOT_ASPECT: 0.3,
}));

function member(over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id: 1, series_id: 5, exposure_id: 101, display_order: 0,
    band_height: 1, y_offset: 0, normalization: "max",
    color_override: null, label_override: null,
    q_window_min: null, q_window_max: null, peak_display: null,
    snapshot: null, is_stale: false, created_by: 1, created_at: null, ...over,
  };
}

function series(over: Partial<Series> = {}): Series {
  return {
    id: 5, title: "LL37 titration", description: null, content_hash: "h",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, forked_from_title: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    ordering_variable: "LL37 : lipid ratio", order_rule: "ascending",
    state: "committed", members: [], samples: [], ...over,
  };
}

function renderAt(id = "5") {
  return render(
    <QueryClientProvider client={makeClient()}>
      <MemoryRouter initialEntries={[`/series/${id}`]}>
        <Routes>
          <Route path="/series/:id" element={<SeriesBuilderPage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("SeriesBuilderPage — read + states", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.seriesQ = { data: undefined, isLoading: false, isError: false };
  });

  it("renders the page shell + title from the series", () => {
    h.seriesQ = { data: series(), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("series-builder-page")).toBeInTheDocument();
    expect(screen.getByText("LL37 titration")).toBeInTheDocument();
  });

  it("shows a skeleton while loading", () => {
    h.seriesQ = { data: undefined, isLoading: true, isError: false };
    const { container } = renderAt();
    expect(container.querySelector('[data-boneyard="series-builder"]')).not.toBeNull();
    expect(screen.queryByTestId("mock-multi-trace-plot")).not.toBeInTheDocument();
  });

  it("shows a not-found state when the series query errors", () => {
    h.seriesQ = { data: undefined, isLoading: false, isError: true };
    renderAt();
    expect(screen.getByTestId("series-builder-error")).toBeInTheDocument();
  });

  it("shows an untitled fallback when title is empty", () => {
    h.seriesQ = { data: series({ title: "" }), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByText(/untitled series/i)).toBeInTheDocument();
  });

  it("composes MultiTracePlot with the series members once loaded", () => {
    h.seriesQ = {
      data: series({ members: [member()] }),
      isLoading: false, isError: false,
    };
    renderAt();
    expect(screen.getByTestId("mock-multi-trace-plot")).toBeInTheDocument();
    // grouping-mode + annotation toggles compose alongside the plot
    expect(screen.getByTestId("grouping-mode")).toBeInTheDocument();
    expect(screen.getByTestId("annotation-toggles")).toBeInTheDocument();
  });

  it("shows the empty-plate state when the series has no members", () => {
    h.seriesQ = { data: series({ members: [] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("series-builder-empty")).toBeInTheDocument();
    expect(screen.queryByTestId("mock-multi-trace-plot")).not.toBeInTheDocument();
  });
});
