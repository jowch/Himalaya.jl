// test/SeriesBuilderPage.actions.test.tsx
//
// Phase 5.2 — SeriesBuilderPage interaction-architecture migration tests.
// Mounts the full TestShell (keyboard layer + InteractionDock) so dock buttons
// AND keyboard-triggered actions flow through the REAL registry path:
//   usePageActions → registry → InteractionDock / useKeyboardLayer → action.run
// Arrow keys (↑/↓) fire on the scope container (not window).
// Alt+↑/↓ and ⌘Z/⌘⇧Z fire on window (keyboard layer).
import { describe, it, expect, vi, beforeEach } from "vitest";
import type React from "react";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { Series, SeriesMember, SeriesSample, CorpusSample, Experiment } from "../src/api";
import { SeriesBuilderPage } from "../src/print/pages/SeriesBuilderPage";
import { InteractionDock } from "../src/print/interaction/InteractionDock";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../src/print/interaction/registry";
import { useListCursor } from "../src/print/interaction/useListCursor";
import { DockStepper } from "../src/print/ui";
import { useAppState } from "../src/state";
import { runCursorContract } from "./interaction/cursorContract";

// ── mock state ───────────────────────────────────────────────────────────────
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

type PickerRow = { sample: { id: number; name?: string; experiment_id?: number }; indexing_exposure_id: number | null };
function pickerRow(sampleId: number, exposureId: number | null, name = `S${sampleId}`, experiment_id = 1): PickerRow {
  return { sample: { id: sampleId, name, experiment_id }, indexing_exposure_id: exposureId };
}

const st = {
  seriesById: new Map<number, Series>(),
  seriesUpdatedAt: 1000,
  corpus: [] as CorpusSample[],
  experiments: [] as Experiment[],
  picker: undefined as PickerRow[] | undefined,
  pickerError: false,
  loading: false,
  save: emptyMut(),
  commit: emptyMut(),
};

vi.mock("../src/queries", () => ({
  useSeries: (id: number | undefined) => ({
    data: id !== undefined ? st.seriesById.get(id) : undefined,
    isLoading: st.loading,
    isError: false,
    dataUpdatedAt: st.seriesUpdatedAt,
  }),
  useSeriesTraces: () => ({ data: {}, isLoading: false, isError: false }),
  useCorpusSamples: () => ({ data: st.corpus, isLoading: false, isError: false }),
  useCorpusPickerSamples: () => ({
    data: st.picker,
    isLoading: st.picker === undefined && !st.pickerError,
    isError: st.pickerError,
  }),
  useExperiments: () => ({ data: st.experiments }),
  useSaveSeries: () => st.save,
  useCommitSeriesPlate: () => st.commit,
}));

vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children, loading, fallback }: { children: React.ReactNode; loading: boolean; fallback?: React.ReactNode }) =>
    loading ? <>{fallback}</> : <>{children}</>,
}));

const navigateSpy = vi.fn();
vi.mock("react-router-dom", async (importOriginal) => {
  const actual = await importOriginal<typeof import("react-router-dom")>();
  return { ...actual, useNavigate: () => navigateSpy };
});

// ── fixtures ─────────────────────────────────────────────────────────────────
function member(id: number): SeriesMember {
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
      effective_peaks: [],
      confirmed_index: null,
      confirmed_phases: [],
      analysis_inputs_hash: "h",
    },
    is_stale: false,
    created_by: null,
    created_at: null,
  };
}

function sample(id: number, sampleId: number, position: number): SeriesSample {
  return { id, series_id: 10, sample_id: sampleId, position, pinned: false, excluded: false };
}

function corpus(id: number, name: string): CorpusSample {
  return { id, experiment_id: 1, name, notes: "", tags: [], q_units: "Å⁻¹" };
}

/** 3-sample series (more room to navigate than 2 for reorder tests). */
function baseSeries(over: Partial<Series> = {}): Series {
  return {
    id: 10,
    title: "LL37 series",
    description: "",
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
    // 3 members, 3 samples in recipe order (position 0 = BOTTOM visually per BU-INVERT)
    members: [member(1), member(2), member(3)],
    samples: [sample(101, 1, 0), sample(102, 2, 1), sample(103, 3, 2)],
    ...over,
  };
}

// ── shell ─────────────────────────────────────────────────────────────────────
function TestShell({ children }: { children: React.ReactNode }): JSX.Element {
  useKeyboardLayer();
  return (
    <>
      {children}
      <InteractionDock />
    </>
  );
}

function renderBuilder() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/series/10"]}>
        <Routes>
          <Route path="/series/:id" element={<TestShell><SeriesBuilderPage /></TestShell>} />
          <Route path="/sample/:id" element={<div data-testid="focus-stub" />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

function getScope(): HTMLElement {
  const scope = document.querySelector("[data-interaction-scope]") as HTMLElement | null;
  if (!scope) throw new Error("No [data-interaction-scope] found");
  return scope;
}

function startEdit(): void {
  act(() => { fireEvent.click(screen.getByRole("button", { name: /^edit$/i })); });
}

function resetDraft(): void {
  act(() => useAppState.getState().discardSeriesDraft());
}

const IDS = [1, 2, 3];

beforeEach(() => {
  navigateSpy.mockClear();
  useInteraction.getState().clearPage();
  st.seriesById = new Map([[10, baseSeries()]]);
  st.seriesUpdatedAt = 1000;
  // Include sample 4 in corpus so AddSamplePicker renders (addable.length > 0).
  st.corpus = [corpus(1, "A"), corpus(2, "B"), corpus(3, "C"), corpus(4, "D")];
  st.picker = [pickerRow(1, 1, "A"), pickerRow(2, 2, "B"), pickerRow(3, 3, "C"), pickerRow(4, 4, "D")];
  st.pickerError = false;
  st.loading = false;
  st.save = emptyMut();
  st.commit = emptyMut();
  resetDraft();
});

// ── tests ─────────────────────────────────────────────────────────────────────
describe("SeriesBuilderPage interaction (task 5.2)", () => {
  // ── cursor identity ─────────────────────────────────────────────────────────
  it("cursor is by sample_id, surviving a recipe edit (was selectedIndex bug)", () => {
    // Verify the ID-based cursor distinguishes itself from old index-based behavior.
    // Set up: cursor on sample_id=1 (position 1/3), then move to sample_id=2 (position 2/3).
    // Test: remove sample_id=1 (first row, recipe index 0).
    // Expected: cursor stays on sample_id=2 (now position 1/2).
    // Old bug: index 1 would remain index 1 (out of bounds after removal of index 0),
    // so it might shift or clamp to sample_id=2 anyway, but at index 0 (WRONG, loses position).
    renderBuilder();
    startEdit();

    // Cursor starts on sample_id=1 (position 1 / 3).
    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("1 / 3");
    expect(screen.getByTestId("dock-identity-name")).toHaveTextContent("A");

    // Move cursor to sample_id=2 (position 2 / 3) by pressing ArrowDown once.
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowDown" }); });

    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("2 / 3");
    expect(screen.getByTestId("dock-identity-name")).toHaveTextContent("B");

    // Remove sample_id=1 (bottom visual row = recipe[0], last Remove button due to BU-INVERT).
    // With old index-based cursor: index 1 would stay at index 1, but now there's no
    // index 1 (only 0 left), so it would clamp/wrap, staying on sample_id=2 but
    // at a different index/position (WRONG — wrong index even if same sample).
    // With ID-based cursor: we stay on sample_id=2 at position 1/2 (CORRECT).
    const removeBtns = screen.getAllByRole("button", { name: /remove/i });
    // Visual order is inverted (BU-INVERT), so removeBtns[2] is the bottom row (recipe[0]).
    act(() => { fireEvent.click(removeBtns[2]!); });

    // After removing sample_id=1, recipe has [sample_id=2, sample_id=3] — 2 samples.
    // Cursor should still be on sample_id=2 (now position 1 / 2, since sample_id=2 moves to recipe[0]).
    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("1 / 2");
    expect(screen.getByTestId("dock-identity-name")).toHaveTextContent("B");
  });

  // ── Enter → Focus ────────────────────────────────────────────────────────────
  it("dock-primary (Focus) opens the cursor member's sample with ?from=series", () => {
    renderBuilder();
    // Cursor starts on sample_id=1.
    act(() => { fireEvent.click(screen.getByTestId("dock-primary")); });
    expect(navigateSpy).toHaveBeenCalledWith("/sample/1?from=series");
  });

  it("ArrowDown on scope moves cursor to next sample, dock-primary opens that sample", () => {
    renderBuilder();
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowDown" }); });

    // Cursor is now on sample_id=2.
    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("2 / 3");

    act(() => { fireEvent.click(screen.getByTestId("dock-primary")); });
    expect(navigateSpy).toHaveBeenCalledWith("/sample/2?from=series");
  });

  // ── identity readout (dockExtra) ─────────────────────────────────────────────
  it("dock renders identity readout (dock-identity-name) via dockExtra when picker provides names", () => {
    renderBuilder();
    // dock-identity-name shows the current cursor sample's name.
    expect(screen.getByTestId("dock-identity-name")).toHaveTextContent("A");
  });

  it("identity updates when cursor advances to the next sample", () => {
    renderBuilder();
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowDown" }); });
    expect(screen.getByTestId("dock-identity-name")).toHaveTextContent("B");
  });

  it("identity shows 'from <experiment>' when experiments are provided", () => {
    st.picker = [
      pickerRow(1, 1, "LL2", 1),
      pickerRow(2, 2, "LL4", 1),
      pickerRow(3, 3, "LL6", 1),
    ] as unknown as typeof st.picker;
    st.experiments = [{ id: 1, name: "Titration A" } as unknown as Experiment];
    renderBuilder();
    expect(screen.getByTestId("dock-identity")).toHaveTextContent(/from Titration A/);
  });

  // ── a key → add-sample picker ─────────────────────────────────────────────
  it("'a' key triggers the add-sample picker (builder-add-sample button clicked)", () => {
    renderBuilder();
    startEdit();
    // The add-sample button must be present (builder is in draft mode).
    const addBtn = screen.getByTestId("builder-add-sample");
    const clickSpy = vi.fn();
    addBtn.addEventListener("click", clickSpy);

    act(() => { fireEvent.keyDown(window, { key: "a" }); });
    expect(clickSpy).toHaveBeenCalledTimes(1);
  });

  // ── Mod+Enter → confirm ───────────────────────────────────────────────────
  it("Mod+Enter is declined when no draft is active (stage is idle but no liveDraft)", () => {
    renderBuilder();
    // No draft yet — confirm should be a no-op.
    act(() => { fireEvent.keyDown(window, { key: "Enter", metaKey: true }); });
    expect(st.save.mutate).not.toHaveBeenCalled();
    expect(st.commit.mutate).not.toHaveBeenCalled();
  });

  it("Mod+Enter triggers confirm when draft is ready and picker is loaded", () => {
    renderBuilder();
    startEdit();

    // Draft is now active and picker is pre-populated (resolverReady = true).
    // Click the on-screen Confirm button to verify the path exists, then test
    // keyboard Mod+Enter triggers the same action.
    const confirmBtn = screen.getByRole("button", { name: /^confirm/i });
    expect(confirmBtn).not.toBeDisabled();

    act(() => { fireEvent.keyDown(window, { key: "Enter", metaKey: true }); });
    // The confirm action fires the save mutation.
    expect(st.save.mutate).toHaveBeenCalledTimes(1);
  });

  // ── Alt+↑/↓ reorder ───────────────────────────────────────────────────────
  it("Alt+ArrowUp moves the cursored recipe row UP visually (recipe index +1)", () => {
    renderBuilder();
    startEdit();

    // Cursor starts at sample_id=1, recipe index 0 (BOTTOM visually).
    // Alt+ArrowUp → recipe 0 → recipe 1 (moves sample_id=1 UP visually).
    act(() => { fireEvent.keyDown(window, { key: "ArrowUp", altKey: true }); });

    const draft = useAppState.getState().seriesDraft!;
    // sample_id=1 should now be at recipe index 1 (moved UP in the figure stack).
    expect(draft.recipe[1]!.sample_id).toBe(1);
    // sample_id=2 should have moved DOWN to recipe index 0.
    expect(draft.recipe[0]!.sample_id).toBe(2);
  });

  it("Alt+ArrowDown moves the cursored recipe row DOWN visually (recipe index -1)", () => {
    renderBuilder();
    startEdit();

    // Move cursor to sample_id=2 (recipe index 1, MIDDLE) first.
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowDown" }); });
    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("2 / 3");

    // Alt+ArrowDown → recipe 1 → recipe 0 (moves sample_id=2 DOWN visually).
    act(() => { fireEvent.keyDown(window, { key: "ArrowDown", altKey: true }); });

    const draft = useAppState.getState().seriesDraft!;
    // sample_id=2 should now be at recipe index 0 (moved DOWN in the figure stack).
    expect(draft.recipe[0]!.sample_id).toBe(2);
    // sample_id=1 should have moved UP to recipe index 1.
    expect(draft.recipe[1]!.sample_id).toBe(1);
  });

  it("cursor follows sample_id after Alt+ArrowUp reorder (stays on same sample)", () => {
    renderBuilder();
    startEdit();

    // Cursor at sample_id=1 (recipe 0). After Alt+↑, sample_id=1 moves to recipe 1.
    act(() => { fireEvent.keyDown(window, { key: "ArrowUp", altKey: true }); });

    // Cursor should still be on sample_id=1 (dock-sample-count reflects new position).
    // sample_id=1 is now at ids index 1 → readout = "2 / 3"
    expect(screen.getByTestId("dock-sample-count")).toHaveTextContent("2 / 3");
  });

  // ── empty list shows only up-link ─────────────────────────────────────────
  it("an empty member list shows ONLY the up-link (no stepper / identity / Focus)", () => {
    st.seriesById = new Map([[10, baseSeries({ members: [], samples: [] })]]);
    renderBuilder();
    expect(screen.getByTestId("dock-up-link")).toBeInTheDocument();
    expect(screen.queryByTestId("dock-sample-count")).toBeNull();
    expect(screen.queryByTestId("dock-identity")).toBeNull();
    expect(screen.queryByTestId("dock-primary")).toBeNull();
  });
});

// ── cursor contract (sample cursor) ──────────────────────────────────────────
// Runs the standard cursor contract on a standalone useListCursor over IDS.
// Builder's cursor is HEADLESS on the page (no rowProps on recipe rows), so
// the contract probe renders its own rows with rowProps for the harness.
runCursorContract("Builder sample cursor", () => {
  const onActivate = vi.fn();
  return {
    ui: (capture) => {
      function Probe(): JSX.Element {
        const cursor = useListCursor({
          ids: IDS,
          onActivate,
          stepperLabel: "Sample",
          stepperTestIdBase: "sample",
          axis: "vertical",
        });
        capture({ cursor, ids: IDS, onActivate });
        return (
          <div data-interaction-scope="">
            {IDS.map((id) => (
              <div key={id} {...cursor.rowProps(id)}>
                sample {id}
              </div>
            ))}
            <DockStepper {...cursor.stepperProps()} />
          </div>
        );
      }
      return <Probe />;
    },
  };
});
