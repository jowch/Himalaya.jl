import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { Series, SeriesMember, SeriesSample, CorpusSample, Trace } from "../../src/api";
import { useAppState } from "../../src/state";

// ── mock data plane ──────────────────────────────────────────────────────────
type MutResult = {
  mutate: ReturnType<typeof vi.fn>;
  isSuccess: boolean;
  isPending: boolean;
  data: unknown;
  error: unknown;
  reset: ReturnType<typeof vi.fn>;
};

function emptyMut(): MutResult {
  return { mutate: vi.fn(), isSuccess: false, isPending: false, data: undefined, error: null, reset: vi.fn() };
}

const state = {
  seriesById: new Map<number, Series>(),
  traces: {} as Record<number, Trace>,
  corpus: [] as CorpusSample[],
  loading: false,
  error: false,
  save: emptyMut(),
  commit: emptyMut(),
};

vi.mock("../../src/queries", () => ({
  useSeries: (id: number | undefined) => ({
    data: id !== undefined ? state.seriesById.get(id) : undefined,
    isLoading: state.loading,
    isError: state.error,
  }),
  useSeriesTraces: (_id: number | undefined) => ({ data: state.traces, isLoading: false }),
  useCorpusSamples: () => ({ data: state.corpus, isLoading: false, isError: false }),
  useSaveSeries: () => state.save,
  useCommitSeriesPlate: () => state.commit,
}));

// boneyard Skeleton: render children when not loading.
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children, loading, fallback }: {
    children: React.ReactNode; loading: boolean; fallback?: React.ReactNode;
  }) => (loading ? <>{fallback}</> : <>{children}</>),
}));

import { SeriesBuilderPage } from "../../src/print/pages/SeriesBuilderPage";

// ── fixtures ─────────────────────────────────────────────────────────────────
function member(id: number, over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id,
    series_id: 10,
    exposure_id: id,
    display_order: id,
    band_height: 1,
    y_offset: 0,
    normalization: "none",
    color_override: null,
    label_override: `ratio ${id}`,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot: {
      effective_peaks: [{ id: 100 + id, q: 0.05 + id * 0.001, intensity: 1000, sharpness: 0.5, source: "auto" }],
      confirmed_index: { id: id, phase: "Ia3d", lattice_d: 195, r_squared: 0.99, ngc: 3, peak_ids: [100 + id] },
      confirmed_phases: [{ phase: "Ia3d", lattice_d: 195 }],
      analysis_inputs_hash: "h",
    },
    is_stale: false,
    created_by: null,
    created_at: null,
    ...over,
  };
}

function seriesSample(id: number, sampleId: number, position: number): SeriesSample {
  return { id, series_id: 10, sample_id: sampleId, position, pinned: false, excluded: false };
}

function corpusSample(id: number, name: string): CorpusSample {
  return {
    id, experiment_id: 1, name, display_name: name, notes: "", tags: [], q_units: "Å⁻¹",
  };
}

function baseSeries(over: Partial<Series> = {}): Series {
  return {
    id: 10,
    title: "LL37 ratio series",
    description: "ratio sweep",
    content_hash: "abc",
    created_by: null,
    created_at: null,
    updated_at: null,
    forked_from_id: null,
    forked_at_hash: null,
    forked_from_title: null,
    view_grouping_mode: null,
    view_show_peak_ticks: null,
    view_show_peak_labels: null,
    ordering_variable: "ratio",
    order_rule: "ascending",
    state: "committed",
    members: [member(1), member(2)],
    samples: [seriesSample(101, 1, 0), seriesSample(102, 2, 1)],
    ...over,
  };
}

function renderPage(): { rerender: () => void } {
  const qc = new QueryClient();
  const tree = (): JSX.Element => (
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/series/10"]}>
        <Routes>
          <Route path="/series/:id" element={<SeriesBuilderPage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>
  );
  const result = render(tree());
  return { rerender: () => result.rerender(tree()) };
}

function resetDraft(): void {
  act(() => useAppState.getState().discardSeriesDraft());
}

beforeEach(() => {
  vi.clearAllMocks();
  state.seriesById = new Map([[10, baseSeries()]]);
  state.traces = {};
  state.corpus = [corpusSample(1, "A"), corpusSample(2, "B"), corpusSample(3, "C")];
  state.loading = false;
  state.error = false;
  state.save = emptyMut();
  state.commit = emptyMut();
  resetDraft();
});

describe("SeriesBuilderPage", () => {
  it("renders the committed read-state: plate title + member rows, no draft, Confirm inert (no Cancel)", () => {
    renderPage();
    expect(screen.getByTestId("series-plate")).toBeInTheDocument();
    // Title surfaces in the editable plate-title field.
    expect((screen.getByLabelText(/series title/i) as HTMLInputElement).value).toBe("LL37 ratio series");
    // MemberList rows for the committed members.
    expect(screen.getAllByTestId("series-member-row")).toHaveLength(2);
    // Read state: "Adjust" is the live entry point; no Cancel (no draft).
    expect(screen.getByRole("button", { name: /adjust/i })).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: /^cancel$/i })).toBeNull();
    expect(useAppState.getState().seriesDraft).toBeNull();
    // The rail's baked "Confirm series" is present but INERT in read state —
    // clicking it does not fire the save chain.
    fireEvent.click(screen.getByRole("button", { name: /confirm series/i }));
    expect(state.save.mutate).not.toHaveBeenCalled();
  });

  it("surfaces a load error distinctly (no plate)", () => {
    state.error = true;
    state.seriesById = new Map();
    renderPage();
    expect(screen.getByText(/couldn't load|could not load/i)).toBeInTheDocument();
    expect(screen.queryByTestId("series-plate")).toBeNull();
  });

  it("editing the title STARTS a draft and reflects the edit", () => {
    renderPage();
    const titleInput = screen.getByLabelText(/series title/i) as HTMLInputElement;
    fireEvent.change(titleInput, { target: { value: "New title" } });
    const draft = useAppState.getState().seriesDraft;
    expect(draft).not.toBeNull();
    expect(draft!.id).toBe(10);
    expect(draft!.title).toBe("New title");
  });

  it("adding a sample STARTS a draft and appends the recipe row", () => {
    renderPage();
    // Add sample C (id 3), the only corpus sample not already in the recipe.
    const select = screen.getByTestId("builder-add-sample-select") as HTMLSelectElement;
    fireEvent.change(select, { target: { value: "3" } });
    const draft = useAppState.getState().seriesDraft;
    expect(draft).not.toBeNull();
    expect(draft!.recipe.map((r) => r.sample_id)).toContain(3);
  });

  it("view controls (offset / scale) change local state and do NOT start a draft", () => {
    renderPage();
    // q-scale toggle on the rail (two segmented controls exist — plate + rail;
    // either is a local-only view control). Use the offset slider, which is
    // unique to the rail.
    const slider = screen.getByLabelText(/trace offset/i);
    fireEvent.change(slider, { target: { value: "0.8" } });
    // A plate scale toggle is also local-only.
    fireEvent.click(screen.getAllByRole("button", { name: /linear q/i })[0]!);
    expect(useAppState.getState().seriesDraft).toBeNull();
  });

  it("Cancel discards the draft with NO request", () => {
    renderPage();
    // Start a draft via the title edit.
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    expect(useAppState.getState().seriesDraft).not.toBeNull();
    fireEvent.click(screen.getByRole("button", { name: /^cancel$/i }));
    expect(useAppState.getState().seriesDraft).toBeNull();
    expect(state.save.mutate).not.toHaveBeenCalled();
    expect(state.commit.mutate).not.toHaveBeenCalled();
  });

  it("Confirm chain: save THEN commit (members from the save response) THEN discard draft", () => {
    const { rerender } = renderPage();
    // Start a draft.
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    const draft = useAppState.getState().seriesDraft!;
    // Click Confirm series.
    fireEvent.click(screen.getByRole("button", { name: /confirm series/i }));
    // save.mutate called with buildSeriesSaveBody(draft) + id.
    expect(state.save.mutate).toHaveBeenCalledTimes(1);
    const saveArg = state.save.mutate.mock.calls[0]![0] as { id: number; title: string };
    expect(saveArg.id).toBe(10);
    expect(saveArg.title).toBe(draft.title);
    // Commit not yet fired.
    expect(state.commit.mutate).not.toHaveBeenCalled();

    // Flip save to success with a FRESH server response whose members differ
    // from the stale draft plate (different exposure ids prove provenance).
    const savedSeries = baseSeries({
      title: "Edited",
      members: [member(7), member(8)],
    });
    state.save = { ...state.save, isSuccess: true, data: savedSeries };
    act(() => rerender());

    // Commit fired with the SAVE RESPONSE members (display_order 7,8), not 1,2.
    expect(state.commit.mutate).toHaveBeenCalledTimes(1);
    const commitArg = state.commit.mutate.mock.calls[0]![0] as { id: number; members: Array<{ exposure_id: number }> };
    expect(commitArg.id).toBe(10);
    expect(commitArg.members.map((m) => m.exposure_id).sort()).toEqual([7, 8]);

    // Flip commit to success → draft discarded, stay on the page.
    state.commit = { ...state.commit, isSuccess: true };
    act(() => rerender());
    expect(useAppState.getState().seriesDraft).toBeNull();
  });

  it("removing a recipe row mutates the draft", () => {
    renderPage();
    // Start a draft.
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    const before = useAppState.getState().seriesDraft!.recipe.length;
    const dismissers = screen.getAllByTestId("builder-recipe-remove");
    fireEvent.click(dismissers[0]!);
    expect(useAppState.getState().seriesDraft!.recipe.length).toBe(before - 1);
  });

  it("reordering a recipe row mutates the draft order", () => {
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    const firstSampleId = useAppState.getState().seriesDraft!.recipe[0]!.sample_id;
    // Move first row down.
    fireEvent.click(screen.getAllByTestId("builder-recipe-down")[0]!);
    expect(useAppState.getState().seriesDraft!.recipe[1]!.sample_id).toBe(firstSampleId);
  });

  it("shows a commit-error notice (role=alert) on commit failure", () => {
    const { rerender } = renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    state.commit = { ...state.commit, error: new Error("boom") };
    act(() => rerender());
    const alert = screen.getByRole("alert");
    expect(alert.textContent).toMatch(/could.?n.?t|failed|error/i);
  });
});
