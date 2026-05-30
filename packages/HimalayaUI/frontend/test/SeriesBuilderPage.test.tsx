import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import type { Series, SeriesMember } from "../src/api";
import { SeriesBuilderPage } from "../src/pages/SeriesBuilderPage";
import { useAppState } from "../src/state";

const h = vi.hoisted(() => ({
  seriesQ: {} as { data: Series | undefined; isLoading: boolean; isError: boolean },
  // Captures the most recent props passed to the on-screen MultiTracePlot so a
  // test can assert the annotation toggles are forwarded to the plot (not just
  // the export spec) — the regression guard for the round-2 blocking bug.
  lastPlotProps: undefined as undefined | {
    showPeakTicks?: boolean; showPeakLabels?: boolean;
    xType?: "log" | "linear"; workingBandFraction?: number;
    representation?: "waterfall" | "heatmap";
    showCrossTraceTracking?: boolean;
  },
}));
vi.mock("../src/queries", () => ({
  useSeries: () => h.seriesQ,
  useMemberTraces: () => new Map(),
  useMemberTracesLoading: () => false,
  useMemberExposures: () => new Map(),
  useMemberSamples: () => new Map(),
}));
// MultiTracePlot touches Observable Plot / ResizeObserver; stub it — but the
// stub records the props it receives so we can assert prop forwarding.
vi.mock("../src/components/MultiTracePlot", () => ({
  MultiTracePlot: (props: {
    showPeakTicks?: boolean; showPeakLabels?: boolean;
    xType?: "log" | "linear"; workingBandFraction?: number;
    representation?: "waterfall" | "heatmap";
    showCrossTraceTracking?: boolean;
  }) => {
    h.lastPlotProps = props;
    return (
      <div
        data-testid="mock-multi-trace-plot"
        data-show-peak-ticks={String(props.showPeakTicks)}
        data-show-peak-labels={String(props.showPeakLabels)}
        data-x-type={String(props.xType)}
        data-representation={String(props.representation)}
        data-cross-trace-tracking={String(props.showCrossTraceTracking)}
      />
    );
  },
  COMPARE_PLOT_ASPECT: 0.3,
  // The page imports offsetToBandFraction from the same module; keep the real
  // implementation so the workingBandFraction assertions exercise the mapping.
  offsetToBandFraction: (offset: number) =>
    0.45 + Math.min(1, Math.max(0, (offset - 0.4) / 1)) * 0.5,
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

  it("mounts the phases-present reading + member rows in the rail (E-9)", () => {
    const indexed = member({
      id: 7, exposure_id: 101, display_order: 0,
      snapshot: {
        effective_peaks: [{ id: 1, q: 0.043, intensity: 1, sharpness: 1, source: "auto" }],
        confirmed_index: { id: 70, phase: "Pn3m", lattice_d: 205, r_squared: 0.99, ngc: -1.5, peak_ids: [1] },
        confirmed_phases: [{ phase: "Pn3m", lattice_d: 205 }],
        assignment_state: "indexed",
        analysis_inputs_hash: "h",
      },
    });
    h.seriesQ = { data: series({ members: [indexed] }), isLoading: false, isError: false };
    renderAt();
    // The derived phases-present reading mounts…
    expect(screen.getByTestId("series-reading")).toBeInTheDocument();
    expect(screen.getByTestId("reading-phase-row")).toBeInTheDocument();
    // …and a member row.
    expect(screen.getByTestId("series-member-row")).toBeInTheDocument();
  });

  it("mounts the rail and toggles full-bleed collapse", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("series-builder-rail")).toBeInTheDocument();
    expect(screen.getByTestId("representation-toggle")).toBeInTheDocument();
    fireEvent.click(screen.getByTestId("rail-collapse-toggle"));
    expect(screen.getByTestId("rail-restore")).toBeInTheDocument();
    expect(screen.queryByTestId("series-builder-rail")).not.toBeInTheDocument();
  });

  it("forwards the annotation toggles to the on-screen MultiTracePlot (round-2 regression)", () => {
    // Known starting point: both annotation flags on.
    useAppState.setState({ showPeakTicks: true, showPeakLabels: true });
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    // Initially forwarded as true.
    expect(screen.getByTestId("mock-multi-trace-plot"))
      .toHaveAttribute("data-show-peak-ticks", "true");
    // Toggling "Peak ticks" off must flow into the PLOT, not just the export
    // spec — the bug was the on-screen plot omitting the prop (defaulting true).
    fireEvent.click(screen.getByTestId("annotation-toggle-peaks"));
    expect(h.lastPlotProps?.showPeakTicks).toBe(false);
    expect(screen.getByTestId("mock-multi-trace-plot"))
      .toHaveAttribute("data-show-peak-ticks", "false");
  });

  it("renders the figure-as-plate container with kicker tags and a caption", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("series-builder-plate")).toBeInTheDocument();
    expect(screen.getByTestId("fig-tags")).toBeInTheDocument();
    expect(screen.getByTestId("fig-caption")).toBeInTheDocument();
  });

  it("renders only the plate header — the outer page header is retired (R8-N1)", () => {
    // Round-2 finding R8-N1: the outer `<header>` at `SeriesBuilderPage.tsx:99-124`
    // duplicated the figure-plate kicker+title at `:284-301`, diluting the
    // figure-as-plate metaphor (~80px of redundant stack). Mockup
    // `series-builder.html:386-396` has only the plate header.
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.queryByTestId("series-builder-header")).not.toBeInTheDocument();
    // Title still renders, but on the plate (not in a separate outer header).
    const title = screen.getByText("LL37 titration");
    expect(title.closest('[data-testid="series-builder-plate"]')).not.toBeNull();
    // Edit button still discoverable in read mode (re-homed into the kicker row).
    expect(screen.getByTestId("series-builder-edit")).toBeInTheDocument();
    // And it lives inside the plate kicker row, not in a separate header strip.
    expect(screen.getByTestId("series-builder-edit").closest('[data-testid="series-builder-plate"]'))
      .not.toBeNull();
  });

  it("forwards the default scale (log) and offset to the plot", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("mock-multi-trace-plot")).toHaveAttribute("data-x-type", "log");
    expect(h.lastPlotProps?.workingBandFraction).toBeCloseTo(0.85, 2);
  });

  it("flips the plot to linear when the scale toggle is switched", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    fireEvent.click(screen.getByTestId("scale-linear"));
    expect(h.lastPlotProps?.xType).toBe("linear");
    expect(screen.getByTestId("fig-tags")).toHaveTextContent("linear q");
  });

  it("changes the plot offset (workingBandFraction) when the slider moves", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    fireEvent.change(screen.getByTestId("offset-slider"), { target: { value: "0.4" } });
    expect(h.lastPlotProps?.workingBandFraction).toBeCloseTo(0.45, 5);
  });

  it("shows the floating offset dock only when the rail is collapsed", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.queryByTestId("offset-dock")).not.toBeInTheDocument();
    fireEvent.click(screen.getByTestId("rail-collapse-toggle"));
    expect(screen.getByTestId("offset-dock")).toBeInTheDocument();
  });

  it("forwards the default representation 'waterfall' to MultiTracePlot", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("mock-multi-trace-plot"))
      .toHaveAttribute("data-representation", "waterfall");
  });

  it("flips MultiTracePlot to heatmap when the heatmap toggle is clicked (#208 wiring)", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.getByTestId("repr-heatmap")).not.toBeDisabled();
    fireEvent.click(screen.getByTestId("repr-heatmap"));
    expect(h.lastPlotProps?.representation).toBe("heatmap");
    expect(screen.getByTestId("mock-multi-trace-plot"))
      .toHaveAttribute("data-representation", "heatmap");
    // The figure-tag also follows the representation.
    expect(screen.getByTestId("fig-tags")).toHaveTextContent("intensity map");
  });

  it("shows a rotated ordering-variable axis label in heatmap mode (R3-Y07)", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    fireEvent.click(screen.getByTestId("repr-heatmap"));
    const axis = screen.getByTestId("heatmap-axis-title");
    expect(axis).toHaveTextContent("LL37 : lipid ratio");
  });

  it("omits the heatmap axis label in waterfall mode (R3-Y07)", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(screen.queryByTestId("heatmap-axis-title")).toBeNull();
  });

  it("toggles cross-trace tracking when the Track-reflections checkbox is clicked (#208 wiring)", () => {
    h.seriesQ = { data: series({ members: [member()] }), isLoading: false, isError: false };
    renderAt();
    expect(h.lastPlotProps?.showCrossTraceTracking).toBe(false);
    fireEvent.click(screen.getByTestId("track-toggle-input"));
    expect(h.lastPlotProps?.showCrossTraceTracking).toBe(true);
  });

  it("changes the coloring mode via GroupingModeToggle (setGroupingMode wired)", () => {
    h.seriesQ = {
      data: series({ view_grouping_mode: null, members: [member()] }),
      isLoading: false, isError: false,
    };
    renderAt();
    // Default seeds to "bySample" (matches effectiveGroupingMode's hard default).
    expect(screen.getByTestId("grouping-mode")).toHaveAttribute("data-mode", "bySample");
    fireEvent.click(screen.getByRole("radio", { name: /by phase/i }));
    // setGroupingMode is wired: the container reflects the new local mode and
    // the plot stays mounted (no crash on re-render with a new groupingMode).
    expect(screen.getByTestId("grouping-mode")).toHaveAttribute("data-mode", "byPhase");
    expect(screen.getByTestId("mock-multi-trace-plot")).toBeInTheDocument();
  });
});
