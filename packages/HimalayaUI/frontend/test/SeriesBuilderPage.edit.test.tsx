import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import type { Series, SeriesMember, CorpusSample } from "../src/api";
import { SeriesBuilderPage } from "../src/pages/SeriesBuilderPage";
import { useAppState } from "../src/state";
import { __resetOptimisticIdForTest } from "../src/lib/queue/optimisticId";

const h = vi.hoisted(() => ({
  seriesQ: {} as { data: Series | undefined; isLoading: boolean; isError: boolean },
  corpusData: [] as CorpusSample[],
  save: { mutate: vi.fn(), isPending: false, isSuccess: false, error: null as unknown },
  commit: { mutate: vi.fn(), isPending: false, isSuccess: false, error: null as unknown, data: undefined as unknown },
}));

vi.mock("../src/queries", () => ({
  useSeries: () => h.seriesQ,
  useMemberTraces: () => new Map(),
  useMemberTracesLoading: () => false,
  useMemberExposures: () => new Map(),
  useMemberSamples: () => new Map(),
  useCorpusSamples: () => ({ data: h.corpusData }),
  useSaveSeries: () => h.save,
  useCommitSeriesPlate: () => h.commit,
}));
vi.mock("../src/components/MultiTracePlot", () => ({
  MultiTracePlot: () => <div data-testid="mock-multi-trace-plot" />,
  COMPARE_PLOT_ASPECT: 0.3,
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
    id: 5, title: "LL37 titration", description: null, content_hash: "sha256:base",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, forked_from_title: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    ordering_variable: "LL37 : lipid ratio", order_rule: "ascending",
    state: "committed", members: [member()],
    samples: [{ id: 11, series_id: 5, sample_id: 10, position: 0, pinned: false, excluded: false }],
    ...over,
  };
}

function corpusSample(id: number, name: string): CorpusSample {
  return {
    id, experiment_id: 1, name, display_name: name, notes: null, tags: [], q_units: "A-1",
  } as CorpusSample;
}

function renderPage() {
  const client = makeClient();
  return render(
    <QueryClientProvider client={client}>
      <MemoryRouter initialEntries={["/series/5"]}>
        <Routes>
          <Route path="/series/:id" element={<SeriesBuilderPage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("SeriesBuilderPage — edit mode (I3.5b)", () => {
  beforeEach(() => {
    __resetOptimisticIdForTest();
    sessionStorage.clear();
    h.seriesQ = { data: series(), isLoading: false, isError: false };
    h.corpusData = [corpusSample(10, "JC042"), corpusSample(20, "JC050")];
    h.save = { mutate: vi.fn(), isPending: false, isSuccess: false, error: null };
    h.commit = { mutate: vi.fn(), isPending: false, isSuccess: false, error: null, data: undefined };
    act(() => useAppState.getState().discardSeriesDraft());
  });

  it("shows an Edit button in read mode; no recipe editor", () => {
    renderPage();
    expect(screen.getByTestId("series-builder-edit")).toBeInTheDocument();
    expect(screen.queryByTestId("series-recipe-editor")).toBeNull();
  });

  it("clicking Edit seeds the draft and reveals the recipe editor", () => {
    renderPage();
    fireEvent.click(screen.getByTestId("series-builder-edit"));
    expect(useAppState.getState().seriesDraft?.id).toBe(5);
    expect(screen.getByTestId("series-recipe-editor")).toBeInTheDocument();
    expect(screen.getByTestId("series-builder-editing-badge")).toBeInTheDocument();
    // The recipe seeds from the series' samples.
    expect(screen.getAllByTestId("recipe-row")).toHaveLength(1);
  });

  it("adding a sample appends a placeholder row optimistically", () => {
    renderPage();
    fireEvent.click(screen.getByTestId("series-builder-edit"));
    fireEvent.change(screen.getByTestId("recipe-add-sample"), { target: { value: "20" } });
    const rows = screen.getAllByTestId("recipe-row");
    expect(rows).toHaveLength(2);
    // The newly-added recipe row carries a negative placeholder id.
    const recipe = useAppState.getState().seriesDraft!.recipe;
    expect(recipe[1].sample_id).toBe(20);
    expect(recipe[1].id).toBeLessThan(0);
  });

  it("Save recipe flushes useSaveSeries with a PATCH-shaped payload (no expected_content_hash)", () => {
    renderPage();
    fireEvent.click(screen.getByTestId("series-builder-edit"));
    fireEvent.click(screen.getByTestId("recipe-save"));
    expect(h.save.mutate).toHaveBeenCalledTimes(1);
    const arg = h.save.mutate.mock.calls[0][0];
    expect(arg.id).toBe(5);
    expect(arg.samples).toEqual([
      { sample_id: 10, position: 0, pinned: false, excluded: false },
    ]);
    expect("expected_content_hash" in arg).toBe(false);
  });

  it("Commit plate flushes useCommitSeriesPlate WITHOUT expected_content_hash (LWW relax, Plan 6a)", () => {
    renderPage();
    fireEvent.click(screen.getByTestId("series-builder-edit"));
    fireEvent.click(screen.getByTestId("recipe-commit"));
    expect(h.commit.mutate).toHaveBeenCalledTimes(1);
    const arg = h.commit.mutate.mock.calls[0][0];
    expect(arg.id).toBe(5);
    expect("expected_content_hash" in arg).toBe(false);
    expect(arg.members[0].exposure_id).toBe(101);
    expect("id" in arg.members[0]).toBe(false);
  });

  it("Cancel discards the draft and returns to read mode", () => {
    renderPage();
    fireEvent.click(screen.getByTestId("series-builder-edit"));
    fireEvent.click(screen.getByTestId("recipe-cancel"));
    expect(useAppState.getState().seriesDraft).toBeNull();
    expect(screen.queryByTestId("series-recipe-editor")).toBeNull();
    expect(screen.getByTestId("series-builder-edit")).toBeInTheDocument();
  });

  it("export controls remain present in edit mode", () => {
    renderPage();
    fireEvent.click(screen.getByTestId("series-builder-edit"));
    expect(screen.getByTestId("rail-export")).toBeInTheDocument();
  });

  it("a zero-member series is editable: the recipe editor mounts and the first sample can be added", () => {
    // Round-1 blocking fix: the empty-plate placeholder must not lock out edit
    // mode. Edit is reachable, the editor mounts, and the first sample adds.
    h.seriesQ = { data: series({ members: [], samples: [] }), isLoading: false, isError: false };
    renderPage();
    expect(screen.getByTestId("series-builder-empty")).toBeInTheDocument();
    fireEvent.click(screen.getByTestId("series-builder-edit"));
    expect(screen.getByTestId("series-recipe-editor")).toBeInTheDocument();
    // The placeholder persists in the plot area while editing.
    expect(screen.getByTestId("series-builder-empty")).toBeInTheDocument();
    fireEvent.change(screen.getByTestId("recipe-add-sample"), { target: { value: "20" } });
    expect(screen.getAllByTestId("recipe-row")).toHaveLength(1);
  });

  it("the empty-plate CTA is a second door into the recipe editor", () => {
    // Besides the header Edit button, the empty plate offers an in-context
    // 'Add the first sample' CTA — it must open the same recipe-edit flow.
    h.seriesQ = { data: series({ members: [], samples: [] }), isLoading: false, isError: false };
    renderPage();
    fireEvent.click(screen.getByTestId("series-builder-empty-cta"));
    expect(useAppState.getState().seriesDraft?.id).toBe(5);
    expect(screen.getByTestId("series-recipe-editor")).toBeInTheDocument();
    expect(screen.getByTestId("series-builder-editing-badge")).toBeInTheDocument();
  });

  it("Commit is disabled while a recipe save is in flight (stale-plate guard)", () => {
    h.save = { mutate: vi.fn(), isPending: true, isSuccess: false, error: null };
    renderPage();
    fireEvent.click(screen.getByTestId("series-builder-edit"));
    expect(screen.getByTestId("recipe-commit")).toBeDisabled();
  });

  it("title + description inputs drive the draft (surfacing the previously-dead setters)", () => {
    renderPage();
    fireEvent.click(screen.getByTestId("series-builder-edit"));
    fireEvent.change(screen.getByTestId("recipe-title"), { target: { value: "Renamed" } });
    fireEvent.change(screen.getByTestId("recipe-description"), { target: { value: "Notes" } });
    const d = useAppState.getState().seriesDraft!;
    expect(d.title).toBe("Renamed");
    expect(d.description).toBe("Notes");
  });
});
