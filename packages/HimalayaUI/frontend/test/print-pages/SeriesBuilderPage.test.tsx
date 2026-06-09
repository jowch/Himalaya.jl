import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, act, within } from "@testing-library/react";
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
    // controls-don't-lie: the rail's baked "Confirm series" is visibly DISABLED
    // in read state (no live draft → no onConfirm), not a live accent button.
    const confirm = screen.getByRole("button", { name: /confirm series/i });
    expect(confirm).toBeDisabled();
    // And clicking it does not fire the save chain.
    fireEvent.click(confirm);
    expect(state.save.mutate).not.toHaveBeenCalled();
  });

  it("lays out as a full-height work·rail grid (rail is a direct grid child, mirrors Focus)", () => {
    renderPage();
    // Item 1: the success body is the full-bleed [work 1fr · rail] grid (no
    // capped PageFrame wrapper), so the rail — a direct grid child — stretches
    // to the row height flush under the header, matching FocusPage.
    const workspace = screen.getByTestId("builder-workspace");
    expect(workspace).toBeInTheDocument();
    // The plate and the rail are both DIRECT children of the grid container so
    // the default align-items:stretch lets the rail fill the grid-row height.
    const rail = screen.getByTestId("builder-rail");
    expect(rail.parentElement).toBe(workspace);
  });

  it("read-state ordering-variable Field is static read-only (no interactive trigger)", () => {
    renderPage();
    const field = screen.getByTestId("field");
    // controls-don't-lie: nothing wires an order-variable list, so the Field is
    // a plain read-only value — not a clickable dropdown trigger.
    expect(field.tagName).not.toBe("BUTTON");
    expect(field).not.toHaveTextContent("▾");
    expect(field).toHaveTextContent("ratio");
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
    // The offset slider is unique to the rail's Display section.
    const slider = screen.getByLabelText(/trace offset/i);
    fireEvent.change(slider, { target: { value: "0.8" } });
    // The single q-scale toggle lives on the plate (the rail's redundant copy
    // was removed) — it's local-only too.
    fireEvent.click(screen.getByRole("button", { name: /linear q/i }));
    expect(useAppState.getState().seriesDraft).toBeNull();
  });

  it("renders exactly ONE q-scale toggle (on the plate, not duplicated on the rail)", () => {
    renderPage();
    // Item 3: the redundant rail Display scale toggle was removed; only the
    // plate's contextual toggle remains.
    expect(screen.getAllByRole("button", { name: /log q/i })).toHaveLength(1);
    expect(screen.getAllByRole("button", { name: /linear q/i })).toHaveLength(1);
  });

  it("with a draft live, 'Adjust' is hidden (redundant no-op), Confirm+Cancel remain", () => {
    renderPage();
    // Read state: Adjust present.
    expect(screen.getByRole("button", { name: /adjust/i })).toBeInTheDocument();
    // Start a draft.
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    // Item 2: Adjust withheld; Confirm enabled + Cancel present.
    expect(screen.queryByRole("button", { name: /adjust/i })).toBeNull();
    expect(screen.getByRole("button", { name: /confirm series/i })).not.toBeDisabled();
    expect(screen.getByRole("button", { name: /^cancel$/i })).toBeInTheDocument();
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
    // With a live draft the rail's "Confirm series" is now ENABLED.
    const confirm = screen.getByRole("button", { name: /confirm series/i });
    expect(confirm).not.toBeDisabled();
    // Click Confirm series.
    fireEvent.click(confirm);
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

  it("each recipe row is draggable and carries a grip handle (label tells the truth)", () => {
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    const rows = screen.getAllByTestId("builder-recipe-row");
    expect(rows.length).toBeGreaterThan(1);
    for (const row of rows) {
      expect(row).toHaveAttribute("draggable", "true");
      expect(within(row).getByTestId("grip-handle")).toBeInTheDocument();
    }
  });

  it("dragging row 1 onto row 0 reorders the draft recipe (additive to the ▲▼ buttons)", () => {
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    const secondSampleId = useAppState.getState().seriesDraft!.recipe[1]!.sample_id;
    const rows = screen.getAllByTestId("builder-recipe-row");
    // Native HTML5 drag: grab row index 1, drop it on row index 0.
    const dataTransfer = {
      effectAllowed: "",
      dropEffect: "",
      setData: vi.fn(),
      getData: vi.fn(),
    };
    fireEvent.dragStart(rows[1]!, { dataTransfer });
    fireEvent.dragOver(rows[0]!, { dataTransfer });
    fireEvent.drop(rows[0]!, { dataTransfer });
    // Row that was at index 1 is now at index 0.
    expect(useAppState.getState().seriesDraft!.recipe[0]!.sample_id).toBe(secondSampleId);
  });

  it("annotation toggles are ARMED when the series has indexed peaks", () => {
    renderPage();
    const group = screen.getByTestId("annotation-toggles");
    const ticks = within(group).getByRole("button", { name: /peak ticks/i });
    const labels = within(group).getByRole("button", { name: /peak labels/i });
    // baseSeries members are indexed with peak anchors → armed + enabled.
    expect(ticks).toHaveAttribute("aria-pressed", "true");
    expect(labels).toHaveAttribute("aria-pressed", "true");
    expect(ticks).not.toBeDisabled();
    expect(labels).not.toBeDisabled();
  });

  it("annotation toggles are INERT (disabled, not armed) on a peakless series", () => {
    // Item 4: a form-factor series has no indexed-peak anchors → nothing to
    // annotate, so the toggles default OFF and disabled regardless of the
    // global showPeakTicks/showPeakLabels defaults (which are true).
    const peakless = baseSeries({
      members: [
        member(1, { snapshot: { effective_peaks: [], confirmed_index: null, confirmed_phases: [], assignment_state: "form_factor", analysis_inputs_hash: "h" } }),
        member(2, { snapshot: { effective_peaks: [], confirmed_index: null, confirmed_phases: [], assignment_state: "form_factor", analysis_inputs_hash: "h" } }),
      ],
    });
    state.seriesById = new Map([[10, peakless]]);
    renderPage();
    const group = screen.getByTestId("annotation-toggles");
    const ticks = within(group).getByRole("button", { name: /peak ticks/i });
    const labels = within(group).getByRole("button", { name: /peak labels/i });
    expect(ticks).toBeDisabled();
    expect(labels).toBeDisabled();
    expect(ticks).not.toHaveAttribute("aria-pressed", "true");
    expect(labels).not.toHaveAttribute("aria-pressed", "true");
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
