import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, act, within, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { Series, SeriesMember, SeriesSample, CorpusSample, Trace, Experiment } from "../../src/api";
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

// Minimal corpus-picker row: the page reads ONLY sample.id +
// indexing_exposure_id (the recipe→plate resolution source, BU-RECIPENOOP).
type PickerRow = { sample: { id: number }; indexing_exposure_id: number | null };
function pickerRow(sampleId: number, exposureId: number | null): PickerRow {
  return { sample: { id: sampleId }, indexing_exposure_id: exposureId };
}

const state = {
  seriesById: new Map<number, Series>(),
  // Mirrors the TanStack query's dataUpdatedAt for the series key. The page's
  // Save→Commit chain watermarks this at Confirm and only commits once the
  // cache delivers a series NEWER than the watermark — so chain tests bump it
  // when they simulate the queue's onSuccess (setQueryData / invalidate→refetch).
  seriesUpdatedAt: 1000,
  traces: {} as Record<number, Trace>,
  tracesLoading: false,
  tracesError: false,
  corpus: [] as CorpusSample[],
  experiments: [] as Experiment[],
  // Corpus picker projection (sample → indexing exposure). `undefined`
  // simulates a not-yet-loaded picker (Confirm must stay gated).
  picker: undefined as PickerRow[] | undefined,
  pickerError: false,
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
    dataUpdatedAt: state.seriesUpdatedAt,
  }),
  useSeriesTraces: (_id: number | undefined) => ({
    data: state.traces,
    isLoading: state.tracesLoading,
    isError: state.tracesError,
  }),
  useCorpusSamples: () => ({ data: state.corpus, isLoading: false, isError: false }),
  useCorpusPickerSamples: () => ({
    data: state.picker,
    isLoading: state.picker === undefined && !state.pickerError,
    isError: state.pickerError,
  }),
  useExperiments: () => ({ data: state.experiments }),
  useSaveSeries: () => state.save,
  useCommitSeriesPlate: () => state.commit,
}));

// Export-spec spy (BU-TOGGLELIE single-source pin): the page builds its
// ExportSpec through this adapter; mocking it lets the test assert the SAME
// Zustand flags that drive the plate annotations reach the export call.
// Only evaluated at export-click time, so the mock is inert everywhere else.
const { specSpy } = vi.hoisted(() => ({
  // A minimally render-viable spec (title included — renderer reads
  // spec.title.primary) so the export path under test completes on the
  // success branch, not via the renderer's error toast.
  specSpy: vi.fn((_args: Record<string, unknown>) => ({
    width: 10,
    height: 10,
    marks: [],
    title: { primary: "spec under test" },
  })),
}));
vi.mock("../../src/lib/figure-export/adapters/multiTraceAdapter", () => ({
  buildMultiTraceExportSpec: specSpy,
}));

// boneyard Skeleton: render children when not loading.
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children, loading, fallback }: {
    children: React.ReactNode; loading: boolean; fallback?: React.ReactNode;
  }) => (loading ? <>{fallback}</> : <>{children}</>),
}));

import { SeriesBuilderPage } from "../../src/print/pages/SeriesBuilderPage";
import { setToastImpl } from "../../src/lib/toast";

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
  state.seriesUpdatedAt = 1000;
  state.traces = {};
  state.tracesLoading = false;
  state.tracesError = false;
  state.corpus = [corpusSample(1, "A"), corpusSample(2, "B"), corpusSample(3, "C")];
  // Default resolution map mirrors the member() fixture (exposure_id == the
  // sample's id) so recipe [sample 1, sample 2] resolves to exposures [1, 2].
  state.picker = [pickerRow(1, 1), pickerRow(2, 2), pickerRow(3, 3)];
  state.pickerError = false;
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
    // Read state: "Edit" is the live entry point; no Cancel (no draft).
    expect(screen.getByRole("button", { name: /^edit$/i })).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: /^cancel$/i })).toBeNull();
    expect(useAppState.getState().seriesDraft).toBeNull();
    // controls-don't-lie: the rail's baked "Save changes" is visibly DISABLED
    // in read state (no live draft → no onConfirm), not a live accent button.
    const confirm = screen.getByRole("button", { name: /save changes/i });
    expect(confirm).toBeDisabled();
    // And clicking it does not fire the save chain.
    fireEvent.click(confirm);
    expect(state.save.mutate).not.toHaveBeenCalled();
  });

  it("normal page exposes the series title as a real h1, not an interactive control nested in a heading (BU-NOHEAD)", () => {
    renderPage();
    // The built-series page's sole top-level heading is the series title, so it
    // is reachable by heading navigation (WCAG 1.3.1 / 2.4.6). Previously the
    // only <h1> WRAPPED the editable Input, so its accessible name was empty and
    // an interactive control was nested inside a heading.
    expect(
      screen.getByRole("heading", { level: 1, name: "LL37 ratio series" }),
    ).toBeInTheDocument();
    // The editable title field stays a distinct control, NOT a heading
    // descendant — correcting the interactive-inside-heading anti-pattern.
    expect(screen.getByLabelText(/series title/i).closest("h1")).toBeNull();
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

  it("adding a sample via the picker STARTS a draft and appends the recipe row (BU-PICKER)", () => {
    renderPage();
    // Open the search-first picker, then add sample C (id 3) — the only corpus
    // sample not already in the recipe.
    fireEvent.click(screen.getByTestId("builder-add-sample"));
    fireEvent.click(screen.getByTestId("add-opt-3"));
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

  it("with a draft live, 'Edit' is hidden (redundant no-op), Save changes+Cancel remain", () => {
    renderPage();
    // Read state: Edit present.
    expect(screen.getByRole("button", { name: /^edit$/i })).toBeInTheDocument();
    // Start a draft.
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    // Item 2: Edit withheld; Save changes enabled + Cancel present.
    expect(screen.queryByRole("button", { name: /^edit$/i })).toBeNull();
    expect(screen.getByRole("button", { name: /save changes/i })).not.toBeDisabled();
    expect(screen.getByRole("button", { name: /^cancel$/i })).toBeInTheDocument();
  });

  it("read state: no draft plate notice; the rail asserts the WYSIWYG default", () => {
    renderPage();
    // No draft → the plate IS the committed figure; no honest-state notice.
    expect(screen.queryByText(/editing the recipe/i)).toBeNull();
    // The rail's default caption may assert WYSIWYG because it's true here.
    expect(
      screen.getByText(
        "The plate above is the figure as it will export. What you compose is what you publish.",
      ),
    ).toBeInTheDocument();
    expect(screen.queryByText(/last confirmed figure/i)).toBeNull();
  });

  it("draft state: the plate carries the lazy-draft notice and the rail caption stops claiming WYSIWYG", () => {
    renderPage();
    // Start a draft (title edit is the lazy-draft entry).
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    // Honest state: the plate says membership/order edits are not on it yet.
    expect(
      screen.getByText(
        "Editing the recipe. Membership and order changes appear here after you confirm.",
      ),
    ).toBeInTheDocument();
    // The rail caption is the honest draft variant; the WYSIWYG default is gone.
    expect(
      screen.getByText(
        "Membership and order on the plate are from the last confirmed figure. Confirm the series to publish your edits.",
      ),
    ).toBeInTheDocument();
    expect(screen.queryByText(/what you compose is what you publish/i)).toBeNull();
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

  it("arms a beforeunload guard while the draft has unsaved changes, disarms on cancel (BU-NAVAWAY-DRAFT)", () => {
    renderPage();
    // Pristine read-state (no draft): closing the tab is not guarded.
    const clean = new Event("beforeunload", { cancelable: true });
    window.dispatchEvent(clean);
    expect(clean.defaultPrevented).toBe(false);

    // A real unsaved change (title edit) arms the guard — the handler
    // preventDefault()s, which is what triggers the browser's leave prompt.
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    const dirty = new Event("beforeunload", { cancelable: true });
    window.dispatchEvent(dirty);
    expect(dirty.defaultPrevented).toBe(true);

    // Cancel discards the draft → the guard disarms (no false prompt on a clean exit).
    fireEvent.click(screen.getByRole("button", { name: /^cancel$/i }));
    const afterCancel = new Event("beforeunload", { cancelable: true });
    window.dispatchEvent(afterCancel);
    expect(afterCancel.defaultPrevented).toBe(false);
  });

  it("Confirm chain: save THEN commit (plate resolved from the saved RECIPE) THEN discard draft", () => {
    const { rerender } = renderPage();
    // Start a draft.
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    const draft = useAppState.getState().seriesDraft!;
    // With a live draft the rail's "Save changes" is now ENABLED.
    const confirm = screen.getByRole("button", { name: /save changes/i });
    expect(confirm).not.toBeDisabled();
    // Click Save changes.
    fireEvent.click(confirm);
    // save.mutate called with buildSeriesSaveBody(draft) + id.
    expect(state.save.mutate).toHaveBeenCalledTimes(1);
    const saveArg = state.save.mutate.mock.calls[0]![0] as { id: number; title: string };
    expect(saveArg.id).toBe(10);
    expect(saveArg.title).toBe(draft.title);
    // Commit not yet fired.
    expect(state.commit.mutate).not.toHaveBeenCalled();

    // Flip save to success with a FRESH server response. Its CACHED members
    // (exposure 7, 8) deliberately differ from the picker resolution of its
    // recipe (samples 1, 2 → exposures 1, 2): the commit body following the
    // RECIPE proves provenance (BU-RECIPENOOP — the PATCH does not rebuild
    // members, so the cached plate is exactly what must NOT be re-posted).
    // HTTP-wins path: saveSeriesMutator.onSuccess setQueryData's the full
    // response into the series cache, bumping dataUpdatedAt — simulate that
    // cache write alongside the mutation flip (the page commits from the
    // CACHE, the only source that is correct on both race outcomes).
    const savedSeries = baseSeries({
      title: "Edited",
      members: [member(7), member(8)],
    });
    state.save = { ...state.save, isSuccess: true, data: savedSeries };
    state.seriesById.set(10, savedSeries);
    state.seriesUpdatedAt = 2000;
    act(() => rerender());

    // Commit fired with the plate resolved from the saved RECIPE (samples 1, 2
    // → exposures 1, 2), not the stale cached members (7, 8).
    expect(state.commit.mutate).toHaveBeenCalledTimes(1);
    const commitArg = state.commit.mutate.mock.calls[0]![0] as { id: number; members: Array<{ exposure_id: number }> };
    expect(commitArg.id).toBe(10);
    expect(commitArg.members.map((m) => m.exposure_id)).toEqual([1, 2]);

    // Flip commit to success → draft discarded, stay on the page.
    state.commit = { ...state.commit, isSuccess: true };
    act(() => rerender());
    expect(useAppState.getState().seriesDraft).toBeNull();
  });

  it("the confirm action names the CURRENT chain step (visible progress), not a static 'Saving…' (BU-PROGRESS)", () => {
    const { rerender } = renderPage();
    // Enter a draft so the chain is live.
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    // Phase 1 — the recipe save is in flight.
    state.save = { ...emptyMut(), isPending: true };
    act(() => rerender());
    expect(
      screen.getByRole("button", { name: /saving order/i }),
    ).toBeInTheDocument();
    // Phase 2 — the save landed; the figure commit is publishing.
    state.save = { ...emptyMut(), isSuccess: true };
    state.commit = { ...emptyMut(), isPending: true };
    act(() => rerender());
    expect(
      screen.getByRole("button", { name: /publishing the figure/i }),
    ).toBeInTheDocument();
  });

  it("a save-step failure and a commit-step failure show DIFFERENT, cause-named errors (BU-PROGRESS / F-ERRSILENT)", () => {
    const { rerender } = renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    // The save (recipe persist) failed → nothing was published.
    state.save = { ...emptyMut(), error: new Error("save boom") };
    act(() => rerender());
    expect(screen.getByRole("alert").textContent).toMatch(/couldn.t save/i);
    // The commit (figure publish) failed AFTER the save landed → the order is
    // saved, only the figure didn't rebuild. Distinct recovery, distinct words.
    state.save = { ...emptyMut(), isSuccess: true };
    state.commit = { ...emptyMut(), error: new Error("commit boom") };
    act(() => rerender());
    expect(screen.getByRole("alert").textContent).toMatch(
      /publishing the figure failed/i,
    );
  });

  it("HTTP-wins ordering: a fresh cache BEFORE save.isSuccess does not commit; the later isSuccess flip does", () => {
    // Pins the merged-single-effect rationale: on HTTP-wins the mutator's
    // setQueryData lands (cache fresh, dataUpdatedAt advanced) BEFORE the
    // mutation flips to success. A premature commit here would publish before
    // the save settled; a split effect keyed only on dataUpdatedAt would
    // conversely stall. Two renders, one assertion each.
    const { rerender } = renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    fireEvent.click(screen.getByRole("button", { name: /save changes/i }));

    // Render N: cache already fresh, mutation still pending.
    const savedSeries = baseSeries({ title: "Edited", members: [member(7), member(8)] });
    state.seriesById.set(10, savedSeries);
    state.seriesUpdatedAt = 2000;
    act(() => rerender());
    expect(state.commit.mutate).not.toHaveBeenCalled();

    // Render N+1: the success flip arrives; the same effect re-checks the
    // (already fresh) cache and commits.
    state.save = { ...state.save, isSuccess: true, data: savedSeries };
    act(() => rerender());
    expect(state.commit.mutate).toHaveBeenCalledTimes(1);
  });

  it("P0 BU-CONFIRMSTALL: SSE-wins partial save response still commits once the cache refetch delivers the fresh series", () => {
    // Live repro: the queue's own-op SSE confirmation resolves save.data as the
    // SYNTHESIZED PARTIAL ({id}, no members/state). The old code gated commit
    // on isFullSeries(save.data) and silently stalled. The fix commits from the
    // SERIES QUERY CACHE: saveSeriesMutator.onSuccess invalidates the key, the
    // refetch lands the canonical full series, and the page watches
    // dataUpdatedAt past the Confirm-time watermark.
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      const { rerender } = renderPage();
      fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
      fireEvent.click(screen.getByRole("button", { name: /save changes/i }));
      expect(state.save.mutate).toHaveBeenCalledTimes(1);

      // SSE wins the race: save.data is the memberless partial shape.
      state.save = { ...state.save, isSuccess: true, data: { id: 10 } };
      act(() => rerender());
      // The partial alone must NOT be committed (no members to commit from).
      expect(state.commit.mutate).not.toHaveBeenCalled();

      // The mutator's invalidate → TanStack refetch → canonical full series
      // lands in the cache with a newer dataUpdatedAt.
      const fresh = baseSeries({ title: "Edited", members: [member(7), member(8)] });
      state.seriesById.set(10, fresh);
      state.seriesUpdatedAt = 2000;
      act(() => rerender());

      // Commit fired with the plate resolved from the FRESH cache's RECIPE
      // (samples 1, 2 → picker exposures 1, 2) — not from the cached members
      // (7, 8), which the PATCH never rebuilds (BU-RECIPENOOP).
      expect(state.commit.mutate).toHaveBeenCalledTimes(1);
      const commitArg = state.commit.mutate.mock.calls[0]![0] as { id: number; members: Array<{ exposure_id: number }> };
      expect(commitArg.id).toBe(10);
      expect(commitArg.members.map((m) => m.exposure_id)).toEqual([1, 2]);

      // Commit success → draft discarded + terminal toast (chain completes).
      state.commit = { ...state.commit, isSuccess: true };
      act(() => rerender());
      expect(useAppState.getState().seriesDraft).toBeNull();
      expect(toast).toHaveBeenCalledWith("Series confirmed", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("awaiting-fresh + series query error: stage resets with a truthful toast, draft preserved, Confirm re-armed", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      const { rerender } = renderPage();
      fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
      fireEvent.click(screen.getByRole("button", { name: /save changes/i }));

      // SSE-wins partial → the chain is awaiting the cache refetch…
      state.save = { ...state.save, isSuccess: true, data: { id: 10 } };
      act(() => rerender());
      expect(state.commit.mutate).not.toHaveBeenCalled();

      // …but the refetch ERRORS. The chain must not stall silently: it resets
      // and tells the truth (the PATCH landed; only the confirm step failed).
      state.error = true;
      act(() => rerender());
      expect(toast).toHaveBeenCalledWith(
        "Order saved, but confirming failed. Try Confirm again.",
        "error",
      );
      expect(state.commit.mutate).not.toHaveBeenCalled();
      // Draft preserved so the user keeps their edits.
      expect(useAppState.getState().seriesDraft).not.toBeNull();

      // Query recovers → Confirm is re-enabled and re-fires the save chain.
      state.error = false;
      act(() => rerender());
      const confirm = screen.getByRole("button", { name: /save changes/i });
      expect(confirm).not.toBeDisabled();
      fireEvent.click(confirm);
      expect(state.save.mutate).toHaveBeenCalledTimes(2);
    } finally {
      setToastImpl(null);
    }
  });

  it("BU-PROGRESS: while the Confirm chain runs the rail control names the live step with aria-busy, still disabled", () => {
    const { rerender } = renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    fireEvent.click(screen.getByRole("button", { name: /save changes/i }));

    // The save is in flight → the control tells the truth: the step-named
    // progressive register + aria-busy, and no live control (no double-submit).
    state.save = { ...state.save, isPending: true };
    act(() => rerender());
    const busyBtn = screen.getByRole("button", { name: /saving order/i });
    expect(busyBtn).toHaveAttribute("aria-busy", "true");
    expect(busyBtn).toBeDisabled();
    expect(screen.queryByRole("button", { name: /save changes/i })).not.toBeInTheDocument();

    // The awaiting-fresh GAP (save settled, commit not yet fired because the
    // cache is still stale) is part of the chain → still busy, now "Confirming…".
    state.save = { ...state.save, isPending: false, isSuccess: true, data: { id: 10 } };
    act(() => rerender());
    expect(screen.getByRole("button", { name: /confirming/i })).toHaveAttribute(
      "aria-busy",
      "true",
    );
  });

  it("BU-PROGRESS: the busy register reverts to the resting label after commit success", () => {
    const { rerender } = renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    fireEvent.click(screen.getByRole("button", { name: /save changes/i }));

    // Save lands + fresh cache → commit fires (the chain is mid-flight).
    const savedSeries = baseSeries({ title: "Edited" });
    state.save = { ...state.save, isSuccess: true, data: savedSeries };
    state.seriesById.set(10, savedSeries);
    state.seriesUpdatedAt = 2000;
    act(() => rerender());
    expect(state.commit.mutate).toHaveBeenCalledTimes(1);
    expect(screen.getByRole("button", { name: /confirming/i })).toBeInTheDocument();

    // Commit success → draft discarded, the control reverts to the resting
    // register with no aria-busy.
    state.commit = { ...state.commit, isSuccess: true };
    act(() => rerender());
    const resting = screen.getByRole("button", { name: /save changes/i });
    expect(resting).not.toHaveAttribute("aria-busy");
    expect(
      screen.queryByRole("button", { name: /saving order|confirming|publishing/i }),
    ).not.toBeInTheDocument();
  });

  it("BU-PROGRESS: the busy register reverts on the stall-exit error (chain reset, Confirm re-armed at rest)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      const { rerender } = renderPage();
      fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
      fireEvent.click(screen.getByRole("button", { name: /save changes/i }));

      // SSE-wins partial → awaiting the cache refetch: busy ("Confirming…").
      state.save = { ...state.save, isSuccess: true, data: { id: 10 } };
      act(() => rerender());
      expect(screen.getByRole("button", { name: /confirming/i })).toBeInTheDocument();

      // The refetch ERRORS (stall exit) and then recovers → the control is
      // back at rest: resting label, no aria-busy, enabled for the retry.
      state.error = true;
      act(() => rerender());
      state.error = false;
      act(() => rerender());
      const resting = screen.getByRole("button", { name: /save changes/i });
      expect(resting).not.toHaveAttribute("aria-busy");
      expect(resting).not.toBeDisabled();
      expect(
        screen.queryByRole("button", { name: /saving order|confirming|publishing/i }),
      ).not.toBeInTheDocument();
    } finally {
      setToastImpl(null);
    }
  });

  it("BU-PROGRESS: a series-query error PRE-DATING the save settle does not stick saving… (saving-clause exit)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      const { rerender } = renderPage();
      fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
      fireEvent.click(screen.getByRole("button", { name: /save changes/i }));

      // The series query errors FIRST (stage ref still "saving")...
      state.error = true;
      act(() => rerender());
      // ...then the save settles as the SSE-wins partial in the same state.
      // The stage ref flips saving -> awaiting-fresh -> idle inside effects
      // that trigger no re-render of their own, so the saving clause itself
      // must exit on (save.isSuccess && seriesQ.isError).
      state.save = { ...state.save, isSuccess: true, data: { id: 10 } };
      act(() => rerender());
      // The page-level error branch (BU-BADID) replaces the body entirely
      // whenever the series query errors, so no rail — and therefore no
      // stuck saving… — can render in any seriesQ.isError state. The
      // saving-clause exit term keeps the DERIVATION truthful regardless.
      expect(
        screen.queryByRole("button", { name: /saving order|confirming|publishing/i }),
      ).not.toBeInTheDocument();
      expect(screen.getByText(/couldn't load this series/i)).toBeInTheDocument();
    } finally {
      setToastImpl(null);
    }
  });

  it("BU-PROGRESS: the busy register reverts IN THE SAME RENDER a save error surfaces (no stuck saving…)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      const { rerender } = renderPage();
      fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
      fireEvent.click(screen.getByRole("button", { name: /save changes/i }));
      state.save = { ...state.save, isPending: true };
      act(() => rerender());
      expect(screen.getByRole("button", { name: /saving order/i })).toBeInTheDocument();

      // The save FAILS. The render that surfaces the role=alert must already
      // show the resting register — the stage ref resets in an effect that
      // triggers no re-render of its own, so a derivation lagging one render
      // would leave a lying busy label next to the failure notice.
      state.save = { ...state.save, isPending: false, error: new Error("boom") };
      act(() => rerender());
      expect(screen.getByRole("alert").textContent).toMatch(/couldn.t save/i);
      const resting = screen.getByRole("button", { name: /save changes/i });
      expect(resting).not.toHaveAttribute("aria-busy");
      expect(
        screen.queryByRole("button", { name: /saving order|confirming|publishing/i }),
      ).not.toBeInTheDocument();
    } finally {
      setToastImpl(null);
    }
  });

  it("stale-cache guard: a series query older than the Confirm watermark does NOT trigger the commit", () => {
    const { rerender } = renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    fireEvent.click(screen.getByRole("button", { name: /save changes/i }));

    // Save lands (SSE-wins partial), but the cache still holds the PRE-save
    // series: dataUpdatedAt has not advanced past the Confirm-time watermark.
    state.save = { ...state.save, isSuccess: true, data: { id: 10 } };
    state.seriesById.set(10, baseSeries({ members: [member(7), member(8)] }));
    // state.seriesUpdatedAt stays at its beforeEach value (== watermark).
    act(() => rerender());
    act(() => rerender());

    // No premature commit from stale data.
    expect(state.commit.mutate).not.toHaveBeenCalled();
  });

  // ── P0 BU-RECIPENOOP: the commit body is the plate RESOLVED FROM THE SAVED
  // RECIPE (fresh.samples + picker), never the cached members the PATCH does
  // not rebuild. ──────────────────────────────────────────────────────────────

  /** Drive a full Confirm → fresh-cache cycle and return the commit body's
   *  members plus the rerender handle (for toast/terminal assertions). */
  function runConfirmChain(fresh: Series): {
    members: Array<{ exposure_id: number; display_order: number } & Record<string, unknown>>;
    rerender: () => void;
  } {
    const { rerender } = renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    fireEvent.click(screen.getByRole("button", { name: /save changes/i }));
    expect(state.save.mutate).toHaveBeenCalledTimes(1);
    // Save lands; the fresh full series (recipe just saved, members stale)
    // reaches the cache past the watermark.
    state.save = { ...state.save, isSuccess: true, data: fresh };
    state.seriesById.set(10, fresh);
    state.seriesUpdatedAt = 2000;
    act(() => rerender());
    expect(state.commit.mutate).toHaveBeenCalledTimes(1);
    const arg = state.commit.mutate.mock.calls[0]![0] as {
      members: Array<{ exposure_id: number; display_order: number } & Record<string, unknown>>;
    };
    return { members: arg.members, rerender };
  }

  it("BU-RECIPENOOP repro: a saved REORDER reaches the commit body in recipe order, not the stale plate order", () => {
    // The PATCH persisted recipe [sample 2, sample 1]; the cached members
    // still carry the OLD order (exposure 1 then 2). The old code posted the
    // members verbatim → byte-identical old plate; the fix follows the recipe.
    const fresh = baseSeries({
      title: "Edited",
      samples: [seriesSample(102, 2, 0), seriesSample(101, 1, 1)],
      members: [member(1, { band_height: 2.5 }), member(2)],
    });
    const { members } = runConfirmChain(fresh);
    expect(members.map((m) => m.exposure_id)).toEqual([2, 1]);
    expect(members.map((m) => m.display_order)).toEqual([0, 1]);
  });

  it("BU-RECIPENOOP repro: an ADDED sample (3 recipe rows, 2 old members) becomes a third member with defaults + resolved exposure", () => {
    const fresh = baseSeries({
      title: "Edited",
      samples: [seriesSample(101, 1, 0), seriesSample(102, 2, 1), seriesSample(103, 3, 2)],
      members: [member(1), member(2)],
    });
    const { members } = runConfirmChain(fresh);
    expect(members.map((m) => m.exposure_id)).toEqual([1, 2, 3]);
    // The new member sends identity only; the commit route fills the same
    // defaults the scoping create path gets (band_height 1.0, y_offset 0,
    // normalization "none") and computes the snapshot server-side.
    expect(Object.keys(members[2]!).sort()).toEqual(["display_order", "exposure_id"]);
  });

  it("BU-RECIPENOOP: a REMOVED sample's member drops from the commit body", () => {
    const fresh = baseSeries({
      title: "Edited",
      samples: [seriesSample(102, 2, 0)],
      members: [member(1), member(2)],
    });
    const { members } = runConfirmChain(fresh);
    expect(members.map((m) => m.exposure_id)).toEqual([2]);
  });

  it("BU-RECIPENOOP carry-over: an old member's label_override/band_height survive a reorder", () => {
    const fresh = baseSeries({
      title: "Edited",
      samples: [seriesSample(102, 2, 0), seriesSample(101, 1, 1)],
      members: [member(1, { label_override: "ratio 1:50", band_height: 2.5 }), member(2)],
    });
    const { members } = runConfirmChain(fresh);
    // Exposure 1 moved to display slot 1 but kept its display props.
    const moved = members.find((m) => m.exposure_id === 1)!;
    expect(moved.display_order).toBe(1);
    expect(moved.label_override).toBe("ratio 1:50");
    expect(moved.band_height).toBe(2.5);
  });

  it("BU-RECIPENOOP policy: an unresolvable sample is left out and the terminal toast says so honestly", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      // Sample 3 has NO usable exposure (picker resolves it to null).
      state.picker = [pickerRow(1, 1), pickerRow(2, 2), pickerRow(3, null)];
      const fresh = baseSeries({
        title: "Edited",
        samples: [seriesSample(101, 1, 0), seriesSample(102, 2, 1), seriesSample(103, 3, 2)],
        members: [member(1), member(2)],
      });
      const { members, rerender } = runConfirmChain(fresh);
      // The resolvable members still commit (blocking would dead-end the
      // series; the backend create path skips such samples the same way).
      expect(members.map((m) => m.exposure_id)).toEqual([1, 2]);
      state.commit = { ...state.commit, isSuccess: true };
      act(() => rerender());
      // The toast must not announce a clean confirm over a partial plate.
      expect(toast).toHaveBeenCalledWith(
        "Confirmed. 1 sample has no usable exposure and was left out.",
        "success",
      );
    } finally {
      setToastImpl(null);
    }
  });

  it("BU-RECIPENOOP gate: Confirm stays DISABLED until the picker (resolution source) has loaded", () => {
    state.picker = undefined; // picker still in flight
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    // Draft is live, but the recipe cannot be resolved yet → Confirm gated.
    const confirm = screen.getByRole("button", { name: /save changes/i });
    expect(confirm).toBeDisabled();
    fireEvent.click(confirm);
    expect(state.save.mutate).not.toHaveBeenCalled();
  });

  it("BU-RECIPENOOP gate: a picker LOAD ERROR surfaces a truthful alert beside the disabled Confirm", () => {
    state.picker = undefined;
    state.pickerError = true;
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
    expect(screen.getByRole("button", { name: /save changes/i })).toBeDisabled();
    expect(screen.getByText(/couldn't load exposure data/i)).toBeInTheDocument();
  });

  it("confirm success announces a 'Series confirmed' toast", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      const { rerender } = renderPage();
      fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "Edited" } });
      fireEvent.click(screen.getByRole("button", { name: /save changes/i }));
      const savedSeries = baseSeries({ title: "Edited" });
      state.save = { ...state.save, isSuccess: true, data: savedSeries };
      // Simulate the mutator's HTTP-wins setQueryData (cache delivers fresh).
      state.seriesById.set(10, savedSeries);
      state.seriesUpdatedAt = 2000;
      act(() => rerender());
      state.commit = { ...state.commit, isSuccess: true };
      act(() => rerender());
      expect(toast).toHaveBeenCalledWith("Series confirmed", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("confirm failure announces an error toast (in addition to the role=alert notice)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      const { rerender } = renderPage();
      fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
      // Click Confirm so the chain is in-flight (stage != "idle"); only then is
      // the stage-guarded error effect armed.
      fireEvent.click(screen.getByRole("button", { name: /save changes/i }));
      state.commit = { ...state.commit, error: new Error("boom") };
      act(() => rerender());
      expect(toast).toHaveBeenCalledWith(expect.stringMatching(/confirm/i), "error");
    } finally {
      setToastImpl(null);
    }
  });

  it("error toast fires once even if BOTH save and commit error resolve together", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      const { rerender } = renderPage();
      fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
      fireEvent.click(screen.getByRole("button", { name: /save changes/i }));
      // Both errors land in the same render cycle.
      state.save = { ...state.save, error: new Error("s") };
      state.commit = { ...state.commit, error: new Error("c") };
      act(() => rerender());
      const errorCalls = toast.mock.calls.filter((c) => c[1] === "error");
      expect(errorCalls).toHaveLength(1);
    } finally {
      setToastImpl(null);
    }
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

  it("BU-EMPTYREMOVE: the last member's Remove is disabled (a series keeps ≥1)", () => {
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    // Remove rows until one remains; the final Remove must be inert so the
    // draft can never be emptied to zero and "saved" as an empty series.
    let removers = screen.getAllByTestId("builder-recipe-remove");
    while (removers.length > 1) {
      const enabled = removers.find((b) => !(b as HTMLButtonElement).disabled);
      expect(enabled).toBeTruthy();
      fireEvent.click(enabled!);
      removers = screen.getAllByTestId("builder-recipe-remove");
    }
    expect(useAppState.getState().seriesDraft!.recipe.length).toBe(1);
    expect(removers[0]).toBeDisabled();
  });

  it("reordering a recipe row mutates the draft order", () => {
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    const firstSampleId = useAppState.getState().seriesDraft!.recipe[0]!.sample_id;
    // Move first row down.
    fireEvent.click(screen.getAllByTestId("builder-recipe-down")[0]!);
    expect(useAppState.getState().seriesDraft!.recipe[1]!.sample_id).toBe(firstSampleId);
  });

  // BU-REORDER-ALT: the unified Alt+↑/↓ reorder power-gesture (shared registry
  // reorderUp/Down) mirrors clicking the row's ▲▼ buttons, so the same gesture
  // reorders on the Builder as on the Scoping worksheet.
  it("Alt+↓ on a recipe row reorders it, mirroring the ▼ Move-down button", () => {
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    const firstSampleId = useAppState.getState().seriesDraft!.recipe[0]!.sample_id;
    // Alt+ArrowDown fired from a control in the top visual row mirrors clicking
    // that row's ▼ Move-down.
    fireEvent.keyDown(screen.getAllByTestId("builder-recipe-down")[0]!, {
      key: "ArrowDown",
      altKey: true,
    });
    expect(useAppState.getState().seriesDraft!.recipe[1]!.sample_id).toBe(firstSampleId);
  });

  it("each recipe row shows the figure trace's label so it cross-references the plate (BU-NAMES)", () => {
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    const rows = screen.getAllByTestId("builder-recipe-row");
    expect(rows).toHaveLength(2);
    // The figure (committed members) labels traces by memberRowLabel — here the
    // members carry label_override "ratio N". The recipe rows must echo that same
    // token (not just the sample name) so a row maps to its plate trace. Visual
    // order is reversed (plate-top first): sample 2 row, then sample 1 row.
    const figureLabels = rows.map(
      (r) => within(r).queryByTestId("builder-recipe-figure-label")?.textContent,
    );
    expect(figureLabels).toEqual(["ratio 2", "ratio 1"]);
    // The sample name is still the primary label (not replaced by the token).
    expect(within(rows[0]!).getByText("B")).toBeInTheDocument();
    expect(within(rows[1]!).getByText("A")).toBeInTheDocument();
  });

  it("a draft-only sample not yet in the committed plate carries NO figure label (BU-NAMES)", () => {
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    // Add sample 3 (picker resolves it to exposure 3, but no committed member has
    // exposure 3, so there is no plate trace to reference yet).
    fireEvent.click(screen.getByTestId("builder-add-sample"));
    fireEvent.click(screen.getByTestId("add-opt-3"));
    const rows = screen.getAllByTestId("builder-recipe-row");
    // The newly-added sample sits at recipe position 2 → visual row 0 (top).
    const newRow = rows[0]!;
    expect(within(newRow).getByText("C")).toBeInTheDocument();
    expect(
      within(newRow).queryByTestId("builder-recipe-figure-label"),
    ).toBeNull();
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

  it("the plate head exposes the figure-export split button (Copy + download menu)", () => {
    // Loaded traces for both members → the figure has data → the gate is open.
    state.traces = {
      1: { q: [0.05, 0.06], I: [10, 9], sigma: [1, 1] },
      2: { q: [0.05, 0.06], I: [8, 7], sigma: [1, 1] },
    };
    renderPage();
    const plate = screen.getByTestId("series-plate");
    expect(within(plate).getByTestId("export-button")).toBeInTheDocument();
    expect(within(plate).getByTestId("export-copy")).toBeInTheDocument();
    // With member traces loaded there's a figure to export → the page gate is
    // open. (export-copy itself stays disabled in JSDOM — no clipboard /
    // OffscreenCanvas capability — so assert the gate on the menu trigger.)
    expect(within(plate).getByTestId("export-menu-trigger")).not.toBeDisabled();
  });

  it("'Copy as PNG' no longer renders anywhere (single contextual export on the plate)", () => {
    renderPage();
    expect(screen.queryByText(/copy as png/i)).toBeNull();
  });

  it("the export control is disabled when the series has no members (nothing to export)", () => {
    state.seriesById = new Map([[10, baseSeries({ members: [], samples: [] })]]);
    renderPage();
    // export-copy would be vacuous here: it's ALWAYS disabled in JSDOM via the
    // capability probes — the menu trigger is the assertion that carries the
    // page gate.
    expect(screen.getByTestId("export-menu-trigger")).toBeDisabled();
  });

  it("the export control is disabled while members have NO trace data (loading / failed fetch)", () => {
    // The beforeEach trace mock is empty: members exist but every waterfall row
    // is padded with EMPTY_TRACE — the exported figure would contain no data,
    // so the WYSIWYG-honest gate stays closed (rows.length alone would lie).
    renderPage();
    expect(screen.getAllByTestId("series-member-row")).toHaveLength(2);
    expect(screen.getByTestId("export-menu-trigger")).toBeDisabled();
  });

  // BU-EXPORT-EMPTYALLOWED: a disabled Export names WHY, distinguishing the
  // three causes the WYSIWYG gate folds together.
  it("a disabled Export reads 'no trace data' when members have empty traces (settled)", () => {
    renderPage(); // beforeEach: empty traces, not loading, not error
    expect(screen.getByTestId("export-disabled-reason")).toHaveTextContent(
      "This series has no trace data to export.",
    );
  });

  it("a disabled Export reads 'still loading' while traces are loading", () => {
    state.tracesLoading = true;
    renderPage();
    expect(screen.getByTestId("export-disabled-reason")).toHaveTextContent(
      "Traces are still loading.",
    );
  });

  it("a disabled Export reads 'couldn't load' when the trace fetch errored", () => {
    state.tracesError = true;
    renderPage();
    expect(screen.getByTestId("export-disabled-reason")).toHaveTextContent(
      "Traces couldn't load, so there's nothing to export.",
    );
  });

  it("an ENABLED Export (traces present) shows no disabled-reason", () => {
    state.traces = {
      1: { q: [0.05, 0.06], I: [10, 9], sigma: [1, 1] },
      2: { q: [0.05, 0.06], I: [8, 7], sigma: [1, 1] },
    };
    renderPage();
    expect(screen.queryByTestId("export-disabled-reason")).toBeNull();
  });

  // ── BU-INVERT: the rail mirrors the plate's vertical order ────────────────

  it("BU-INVERT read mode: the rail lists members REVERSED so its top row is the plate's top trace", () => {
    renderPage();
    // Members display order is [ratio 1, ratio 2]; the waterfall paints display
    // order bottom-up, so the plate TOP is "ratio 2". The MemberList contract
    // ("page reverses so top = high variable") means the rail leads with it.
    const rows = screen.getAllByTestId("series-member-row");
    expect(rows[0]).toHaveTextContent("ratio 2");
    expect(rows[1]).toHaveTextContent("ratio 1");
    // The plate's geometry agrees: member 2's wf-row sits ABOVE member 1's.
    const plate = screen.getByTestId("series-plate");
    const topOf = (key: string): number =>
      Number((plate.querySelector(`[data-role="wf-row"][data-key="${key}"]`) as HTMLElement)
        .style.top.replace("px", ""));
    expect(topOf("2")).toBeLessThan(topOf("1"));
  });

  it("BU-INVERT hover-sync: hovering the rail's TOP row lights the plate's TOP trace", () => {
    renderPage();
    const first = screen.getAllByTestId("series-member-row")[0]!;
    fireEvent.mouseEnter(first);
    const hot = document.querySelector('[data-role="wf-row"][data-hot="true"]');
    // Member 2 is the plate top (largest display order, painted highest).
    expect(hot?.getAttribute("data-key")).toBe("2");
  });

  it("BU-INVERT draft mode: recipe rows render plate-top-first; 'Move down' on the FIRST VISUAL row moves the trace toward the plate BOTTOM", () => {
    // Three rows so the permutation discriminates: the pre-fix top-down editor
    // moving its first row (A) down produced [B,A,C]; the plate-mirroring
    // editor's first visual row is C (recipe last) and moving it down gives
    // [A,C,B].
    state.seriesById = new Map([[10, baseSeries({
      samples: [seriesSample(101, 1, 0), seriesSample(102, 2, 1), seriesSample(103, 3, 2)],
    })]]);
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    const rows = screen.getAllByTestId("builder-recipe-row");
    // Visual order mirrors the plate: recipe position 2 (sample C) on top.
    expect(rows[0]).toHaveTextContent("C");
    expect(rows[1]).toHaveTextContent("B");
    expect(rows[2]).toHaveTextContent("A");
    // Endpoint guards in VISUAL terms: top row can't move up, bottom can't move down.
    expect(within(rows[0]!).getByTestId("builder-recipe-up")).toBeDisabled();
    expect(within(rows[2]!).getByTestId("builder-recipe-down")).toBeDisabled();
    // Move down on the top visual row (C): recipe [A,B,C] → [A,C,B].
    fireEvent.click(within(rows[0]!).getByTestId("builder-recipe-down"));
    expect(useAppState.getState().seriesDraft!.recipe.map((r) => r.sample_id)).toEqual([1, 3, 2]);
  });

  it("BU-INVERT draft drag shares the same visual mapping (drop B onto C → B becomes the plate top)", () => {
    state.seriesById = new Map([[10, baseSeries({
      samples: [seriesSample(101, 1, 0), seriesSample(102, 2, 1), seriesSample(103, 3, 2)],
    })]]);
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    const rows = screen.getAllByTestId("builder-recipe-row");
    const dataTransfer = { effectAllowed: "", dropEffect: "", setData: vi.fn(), getData: vi.fn() };
    // Visual rows are [C, B, A]; drag B (visual 1) onto C (visual 0) = move B
    // to the visual top = recipe end: [A,B,C] → [A,C,B].
    fireEvent.dragStart(rows[1]!, { dataTransfer });
    fireEvent.dragOver(rows[0]!, { dataTransfer });
    fireEvent.drop(rows[0]!, { dataTransfer });
    expect(useAppState.getState().seriesDraft!.recipe.map((r) => r.sample_id)).toEqual([1, 3, 2]);
  });

  // ── BU-MOVEFOCUS: a move that self-disables a button keeps focus (WCAG 2.4.3) ──

  it("BU-MOVEFOCUS: moving a row to the TOP redirects focus from the now-disabled 'Move up' to its 'Move down' (not <body>)", () => {
    // Visual order mirrors the plate: recipe [A,B,C] renders [C,B,A]. Press
    // 'Move up' on visual row 1 (B) → B reaches the top, its 'Move up'
    // disables. A keyboard user must not be dumped onto document.body; focus
    // moves to B's still-actionable sibling, 'Move down'.
    state.seriesById = new Map([[10, baseSeries({
      samples: [seriesSample(101, 1, 0), seriesSample(102, 2, 1), seriesSample(103, 3, 2)],
    })]]);
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    let rows = screen.getAllByTestId("builder-recipe-row");
    // Visual [C, B, A]; move B (visual index 1) up.
    const bUp = within(rows[1]!).getByTestId("builder-recipe-up");
    act(() => bUp.focus());
    fireEvent.click(bUp);
    // Recipe [A,B,C] → [A,C,B]; B is now the top visual row.
    rows = screen.getAllByTestId("builder-recipe-row");
    expect(rows[0]).toHaveTextContent("B");
    const bUpNow = within(rows[0]!).getByTestId("builder-recipe-up");
    const bDownNow = within(rows[0]!).getByTestId("builder-recipe-down");
    expect(bUpNow).toBeDisabled();
    expect(document.activeElement).toBe(bDownNow);
    expect(document.activeElement).not.toBe(document.body);
  });

  it("BU-MOVEFOCUS: moving a row to the BOTTOM redirects focus from the now-disabled 'Move down' to its 'Move up' (not <body>)", () => {
    // Symmetric bottom case: press 'Move down' on visual row 1 (B) → B reaches
    // the bottom, its 'Move down' disables, focus moves to 'Move up'.
    state.seriesById = new Map([[10, baseSeries({
      samples: [seriesSample(101, 1, 0), seriesSample(102, 2, 1), seriesSample(103, 3, 2)],
    })]]);
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    let rows = screen.getAllByTestId("builder-recipe-row");
    // Visual [C, B, A]; move B (visual index 1) down.
    const bDown = within(rows[1]!).getByTestId("builder-recipe-down");
    act(() => bDown.focus());
    fireEvent.click(bDown);
    // Recipe [A,B,C] → [B,A,C]; visual order becomes [C,A,B]; B is now bottom.
    rows = screen.getAllByTestId("builder-recipe-row");
    expect(rows[2]).toHaveTextContent("B");
    const bUpNow = within(rows[2]!).getByTestId("builder-recipe-up");
    const bDownNow = within(rows[2]!).getByTestId("builder-recipe-down");
    expect(bDownNow).toBeDisabled();
    expect(document.activeElement).toBe(bUpNow);
    expect(document.activeElement).not.toBe(document.body);
  });

  it("BU-MOVEFOCUS: a non-extreme move does NOT steal focus to a sibling (focus stays on the activated 'Move up')", () => {
    // Four rows so a move can land mid-stack; the just-activated button stays
    // actionable, so the focus-redirect must NOT fire — the row keeps focus on
    // the button the user pressed.
    state.seriesById = new Map([[10, baseSeries({
      samples: [
        seriesSample(101, 1, 0), seriesSample(102, 2, 1),
        seriesSample(103, 3, 2), seriesSample(104, 4, 3),
      ],
    })]]);
    state.corpus = [corpusSample(1, "A"), corpusSample(2, "B"), corpusSample(3, "C"), corpusSample(4, "D")];
    state.picker = [pickerRow(1, 1), pickerRow(2, 2), pickerRow(3, 3), pickerRow(4, 4)];
    renderPage();
    fireEvent.change(screen.getByLabelText(/series title/i), { target: { value: "x" } });
    let rows = screen.getAllByTestId("builder-recipe-row");
    // Visual [D, C, B, A]; move C (visual index 1) up to visual index 0... that
    // WOULD be an extreme. Instead move B (visual index 2) up to index 1 — a
    // non-extreme landing, so 'Move up' stays enabled and keeps focus.
    const bUp = within(rows[2]!).getByTestId("builder-recipe-up");
    act(() => bUp.focus());
    fireEvent.click(bUp);
    rows = screen.getAllByTestId("builder-recipe-row");
    // B moved up one visual slot → now visual index 1.
    expect(rows[1]).toHaveTextContent("B");
    const bUpNow = within(rows[1]!).getByTestId("builder-recipe-up");
    expect(bUpNow).not.toBeDisabled();
    expect(document.activeElement).toBe(bUpNow);
  });

  // ── BU-TOGGLELIE: the annotation toggles drive the PLATE ──────────────────

  // Traces covering the member anchors (q = 0.051 / 0.052) so glyphs render.
  function loadTraces(): void {
    state.traces = {
      1: { q: [0.04, 0.051, 0.07], I: [100, 60, 12], sigma: [1, 1, 1] },
      2: { q: [0.04, 0.052, 0.07], I: [90, 55, 10], sigma: [1, 1, 1] },
    };
  }

  it("BU-TOGGLELIE: switching 'Peak ticks' OFF removes the peak glyphs from the plate", () => {
    loadTraces();
    renderPage();
    try {
      const plate = screen.getByTestId("series-plate");
      expect(plate.querySelectorAll('[data-role="peak-glyph"]').length).toBeGreaterThan(0);
      fireEvent.click(
        within(screen.getByTestId("annotation-toggles")).getByRole("button", { name: /peak ticks/i }),
      );
      expect(plate.querySelectorAll('[data-role="peak-glyph"]')).toHaveLength(0);
    } finally {
      // The flags live in the real (module-singleton) Zustand store; restore.
      act(() => useAppState.setState({ showPeakTicks: true }));
    }
  });

  it("BU-TOGGLELIE: 'Peak labels' ON renders labels at the anchors; OFF removes them", () => {
    loadTraces();
    renderPage();
    try {
      const plate = screen.getByTestId("series-plate");
      // Zustand default is ON → labels render (one anchor per member → "1").
      const labels = plate.querySelectorAll('[data-role="peak-label"]');
      expect(labels.length).toBeGreaterThan(0);
      expect(labels[0]!.textContent).toBe("1");
      fireEvent.click(
        within(screen.getByTestId("annotation-toggles")).getByRole("button", { name: /peak labels/i }),
      );
      expect(plate.querySelectorAll('[data-role="peak-label"]')).toHaveLength(0);
    } finally {
      act(() => useAppState.setState({ showPeakLabels: true }));
    }
  });

  it("BU-TOGGLELIE single source: the export spec is built from the SAME flags that drive the plate", async () => {
    loadTraces();
    renderPage();
    // JSDOM has no URL.createObjectURL; stub it so the download path completes
    // on the success branch instead of logging the renderer's error toast.
    const createObjectURL = vi.fn(() => "blob:test");
    const revokeObjectURL = vi.fn();
    Object.assign(URL, { createObjectURL, revokeObjectURL });
    try {
      // Toggle labels OFF on the plate, then export: the spec must carry the
      // post-toggle flags (ticks true, labels false) — plate and export cannot
      // diverge on these axes because both read one Zustand pair.
      fireEvent.click(
        within(screen.getByTestId("annotation-toggles")).getByRole("button", { name: /peak labels/i }),
      );
      expect(screen.getByTestId("series-plate").querySelectorAll('[data-role="peak-label"]')).toHaveLength(0);
      fireEvent.click(screen.getByTestId("export-menu-trigger"));
      fireEvent.click(screen.getByRole("menuitem", { name: /download as svg/i }));
      await waitFor(() => expect(specSpy).toHaveBeenCalled());
      expect(specSpy.mock.calls.at(-1)![0]).toMatchObject({
        showPeakTicks: true,
        showPeakLabels: false,
      });
    } finally {
      act(() => useAppState.setState({ showPeakLabels: true }));
    }
  });

  it("BU-EXPORTDIVERGE: the export spec carries the plate's grouping/scale/offset/domain/label register", async () => {
    loadTraces();
    renderPage();
    const createObjectURL = vi.fn(() => "blob:test");
    const revokeObjectURL = vi.fn();
    Object.assign(URL, { createObjectURL, revokeObjectURL });
    // Default plate state: phase-colored traces, log q, offset ×1.20.
    fireEvent.click(screen.getByTestId("export-menu-trigger"));
    fireEvent.click(screen.getByRole("menuitem", { name: /download as svg/i }));
    await waitFor(() => expect(specSpy).toHaveBeenCalled());
    const args = specSpy.mock.calls.at(-1)![0] as Record<string, unknown>;
    // The plate footnote promises WYSIWYG — every axis the plate renders
    // flows into the spec from the same sources the plate reads.
    expect(args).toMatchObject({
      groupingMode: "byPhase",
      xType: "log",
      offsetScale: 1.2,
      preset: "clean",
    });
    // Label register: the plate's (label_override ?? "exp N"), per member.
    const labels = args.displayLabelByMemberId as Map<number, string>;
    expect(labels.get(1)).toBe("ratio 1");
    expect(labels.get(2)).toBe("ratio 2");
    // The plate's padded q-domain, not the old hardcoded null.
    expect(Array.isArray(args.xDomain)).toBe(true);

    // Flip the plate to linear q and export again — the spec follows.
    fireEvent.click(screen.getByRole("button", { name: /linear q/i }));
    fireEvent.click(screen.getByTestId("export-menu-trigger"));
    fireEvent.click(screen.getByRole("menuitem", { name: /download as svg/i }));
    await waitFor(() =>
      expect((specSpy.mock.calls.at(-1)![0] as { xType?: string }).xType).toBe("linear"),
    );
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

// ── BU-BADID: a non-numeric :id must exit honestly, not hang on the skeleton ──
// A NaN id is knowable synchronously from the route param, so the page must
// render the SAME not-found surface a missing numeric id gets (the FocusPage
// routeStatus "unknown" precedent), never the perpetual loading skeleton.
describe("SeriesBuilderPage with a non-numeric route id", () => {
  function renderAt(path: string): void {
    const qc = new QueryClient();
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={[path]}>
          <Routes>
            <Route path="/series" element={<div data-testid="folio-stub" />} />
            <Route path="/series/:id" element={<SeriesBuilderPage />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
  }

  it("/series/abc renders the not-found EmptyState, not the loading skeleton", () => {
    renderAt("/series/abc");
    // The honest exit: the shared not-found surface (same branch as a missing
    // numeric id), with a way back to the folio.
    expect(screen.getByTestId("builder-not-found")).toBeInTheDocument();
    expect(screen.getByText(/couldn't load|could not load/i)).toBeInTheDocument();
    // Red-check pins: the old code sat on the skeleton fallback forever
    // (disabled query → isLoading/data never settle) and never showed a plate.
    expect(screen.queryByText(/loading series/i)).toBeNull();
    expect(screen.queryByTestId("series-plate")).toBeNull();
  });

  it("not-found exposes its sole heading as an h1, not a level-skipped h2 (FO-NFHEAD)", () => {
    // The not-found branch returns before the SeriesPlate (whose editable title
    // is an Input, not a heading), so EmptyState's title is the page's only
    // heading and must be the top level (WCAG 1.3.1).
    renderAt("/series/abc");
    expect(
      screen.getByRole("heading", { level: 1, name: /couldn't load this series/i }),
    ).toBeInTheDocument();
    expect(screen.queryByRole("heading", { level: 2, name: /couldn't load this series/i })).toBeNull();
  });

  it("the not-found CTA navigates back to the series folio (/series)", () => {
    renderAt("/series/abc");
    fireEvent.click(screen.getByRole("button", { name: /back to the series folio/i }));
    expect(screen.getByTestId("folio-stub")).toBeInTheDocument();
  });

  it("a numeric-missing id (query error) renders the same not-found surface", () => {
    state.error = true;
    state.seriesById = new Map();
    renderAt("/series/999");
    expect(screen.getByTestId("builder-not-found")).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /back to the series folio/i })).toBeInTheDocument();
    expect(screen.queryByTestId("series-plate")).toBeNull();
  });
});
