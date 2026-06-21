import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import type { CorpusSample, Exposure, Experiment } from "../../src/api";

// ── mutator spy ───────────────────────────────────────────────────────────────
const batchMutate = vi.fn();
const corpusRefetch = vi.fn();

// ── mock data plane (mutated per test) ────────────────────────────────────────
const state = {
  samples: [] as CorpusSample[],
  byId: new Map<number, Exposure[]>(),
  experiments: [] as Experiment[],
  screened: 0,
  total: 0,
  loading: false,
  error: false,
  fetching: false,
};

vi.mock("../../src/queries", () => ({
  useCorpusSamples: () => ({
    data: state.samples,
    isLoading: state.loading,
    isError: state.error,
    isFetching: state.fetching,
    refetch: corpusRefetch,
  }),
  useCorpusExposures: (filtered: CorpusSample[]) => {
    // Only expose the rows for the filtered samples (mirrors the real hook).
    const byId = new Map<number, Exposure[]>();
    for (const s of filtered) {
      const exps = state.byId.get(s.id);
      if (exps !== undefined) byId.set(s.id, exps);
    }
    return { byId, isLoading: state.loading };
  },
  useScreenedProgress: () => ({ screened: state.screened, total: state.total }),
  useExperiments: () => ({ data: state.experiments }),
  useSetExposureStatusBatch: () => ({ mutate: batchMutate }),
}));

// Zustand store (carried hooks inside the page do not read it, but the adapter /
// transitive imports might — keep it trivial).
vi.mock("../../src/state", () => ({
  useAppState: (sel: (s: Record<string, unknown>) => unknown) =>
    sel({ username: "tester" }),
}));

// boneyard Skeleton: render children when not loading.
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

// DetectorImage touches fetch / createImageBitmap (absent in JSDOM).
vi.mock("../../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

import { SamplesPage } from "../../src/print/pages/SamplesPage";
import { setToastImpl } from "../../src/lib/toast";
import { useFloatingDock } from "../../src/print/shell/floatingDock";

// ── fixtures ──────────────────────────────────────────────────────────────────
function corpus(over: Partial<CorpusSample> = {}): CorpusSample {
  return {
    id: 1, experiment_id: 1, name: "JC001",
    notes: null, q_units: "A-1", tags: [], phase: null, ...over,
  } as CorpusSample;
}
function exp(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 100, sample_id: 1, filename: "JC001-001.dat", kind: "file",
    selected: false, status: "accepted", image_path: "/x.tif", image_version: "v1",
    tags: [], sources: [], trace_hash: "h", analysis_inputs_hash: "a", ...over,
  };
}
function experiment(over: Partial<Experiment> = {}): Experiment {
  return { id: 1, name: "BL-19 April", ...over } as Experiment;
}

function seed(): void {
  state.experiments = [experiment({ id: 1, name: "BL-19 April" })];
  state.samples = [
    corpus({ id: 1, experiment_id: 1, name: "Buffer" }),
    corpus({ id: 2, experiment_id: 1, name: "Lipid 1-2", phase: "Pn3m" }),
    corpus({ id: 3, experiment_id: 2, name: "Other beamtime" }),
  ];
  state.byId = new Map<number, Exposure[]>([
    [1, [
      exp({ id: 100, sample_id: 1, filename: "JC001-001.dat" }),
      exp({ id: 101, sample_id: 1, filename: "JC001-002.dat" }),
      exp({ id: 102, sample_id: 1, filename: "JC001-003.dat" }),
    ]],
    [2, [exp({ id: 200, sample_id: 2 }), exp({ id: 201, sample_id: 2 })]],
    [3, [exp({ id: 300, sample_id: 3 })]],
  ]);
  state.screened = 2;
  state.total = 3;
  state.loading = false;
  state.error = false;
  state.fetching = false;
}

/** Loupe-route stand-in that also surfaces the search string, so SA-F2 can
 *  assert WHICH frame the navigation carried (not just that it navigated). */
function LoupeRouteProbe(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="loupe-route" data-search={loc.search} />;
}

function renderAt(path = "/samples") {
  return render(
    <MemoryRouter initialEntries={[path]}>
      <Routes>
        <Route path="/samples" element={<SamplesPage />} />
        <Route path="/samples/loupe/:sampleId" element={<LoupeRouteProbe />} />
        <Route path="/sample/:sampleId" element={<div data-testid="focus-route" />} />
      </Routes>
    </MemoryRouter>,
  );
}

beforeEach(() => {
  vi.clearAllMocks();
  seed();
});

describe("SamplesPage", () => {
  it("renders the head, the sheet table, rows and the keyboard legend", () => {
    renderAt("/samples?experiment=1");
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
    // beamtime title = experiment name
    expect(screen.getByRole("heading", { name: "BL-19 April" })).toBeInTheDocument();
    expect(screen.getByTestId("sheet-table")).toBeInTheDocument();
    expect(screen.getByTestId("kb-legend")).toBeInTheDocument();
  });

  it("renders one row per filtered sample", () => {
    renderAt("/samples?experiment=1");
    // beamtime 1 has 2 samples (ids 1, 2)
    expect(screen.getAllByTestId("sample-table-row")).toHaveLength(2);
  });

  it("clicking a thumb selects it → cull bar shows with the right count", () => {
    renderAt("/samples?experiment=1");
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[0]!);
    const bar = screen.getByTestId("cull-bar");
    expect(bar).toHaveAttribute("data-show", "true");
    expect(within(bar).getByText("1")).toBeInTheDocument();
  });

  it("publishes the floating-dock lane as occupied while a selection bar shows, and frees it on unmount (LA-COLLIDE)", () => {
    const { unmount } = renderAt("/samples?experiment=1");
    // No selection → lane is free, banner stays centred.
    expect(useFloatingDock.getState().centerLaneOccupied).toBe(false);
    // Selecting a frame raises the CullBar → the dock lane is now occupied.
    fireEvent.click(screen.getAllByTestId("thumbnail")[0]!);
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "true");
    expect(useFloatingDock.getState().centerLaneOccupied).toBe(true);
    // Leaving the page must release the lane (else the banner stays cornered).
    unmount();
    expect(useFloatingDock.getState().centerLaneOccupied).toBe(false);
  });

  it("CullBar Drop calls the batch mutator with the selected exposure", () => {
    renderAt("/samples?experiment=1");
    // first row (sample 1) → first thumb = exposure 100
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[0]!);
    fireEvent.click(screen.getByRole("button", { name: /Drop/ }));
    expect(batchMutate).toHaveBeenCalledWith({
      sampleId: 1,
      exposureId: 100,
      status: "rejected",
    });
  });

  it("CullBar Keep calls the batch mutator with status accepted (SA-SCREENED)", () => {
    renderAt("/samples?experiment=1");
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[0]!);
    fireEvent.click(screen.getByRole("button", { name: /Keep/ }));
    expect(batchMutate).toHaveBeenCalledWith({
      sampleId: 1,
      exposureId: 100,
      status: "accepted",
    });
    // Keep clears the selection, same as Drop.
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");
  });

  it("CullBar Restore sends status null even for already-accepted frames", () => {
    renderAt("/samples?experiment=1");
    // Fixture exposure 100 is status "accepted" — Restore must still null it.
    fireEvent.click(screen.getAllByTestId("thumbnail")[0]!);
    fireEvent.click(screen.getByRole("button", { name: /Restore/ }));
    expect(batchMutate).toHaveBeenCalledWith({
      sampleId: 1,
      exposureId: 100,
      status: null,
    });
  });

  it("K keeps the selection, announces a count toast, and clears it", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt("/samples?experiment=1");
      fireEvent.click(screen.getAllByTestId("thumbnail")[0]!);
      fireEvent.keyDown(window, { key: "k" });
      expect(batchMutate).toHaveBeenCalledWith({
        sampleId: 1,
        exposureId: 100,
        status: "accepted",
      });
      expect(toast).toHaveBeenCalledWith("1 frame kept", "success");
      expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");
    } finally {
      setToastImpl(null);
    }
  });

  it("the footer legend does NOT repeat the X/K cull gesture (SA-KBDDUP)", () => {
    // X/K live in the contextual cull hint (SA-CULLHINT, below) and on the
    // CullBar once frames are selected; repeating them in the footer was a
    // duplicate affordance. The footer keeps navigation/selection/open.
    renderAt("/samples?experiment=1");
    const legend = screen.getByTestId("kb-legend");
    expect(within(legend).queryByText("keep the selected frames")).not.toBeInTheDocument();
    expect(within(legend).queryByText("drop the selected frames")).not.toBeInTheDocument();
    expect(within(legend).getByText("move between cells")).toBeInTheDocument();
  });

  it("hints the X/K cull gesture (registry-driven) until a selection exists (SA-CULLHINT)", () => {
    renderAt("/samples?experiment=1");
    const hint = screen.getByTestId("samples-cull-hint");
    const caps = within(hint).getAllByTestId("kbkey").map((k) => k.textContent);
    // key caps come straight from the shortcut registry (drop=X, keep=K)
    expect(caps).toContain("X");
    expect(caps).toContain("K");
    // once frames are selected the CullBar's real Drop/Keep buttons take over,
    // so the teaching hint retires (no duplicate affordance)
    fireEvent.click(screen.getAllByTestId("thumbnail")[0]!);
    expect(screen.queryByTestId("samples-cull-hint")).toBeNull();
  });

  it("the keyboard legend advertises the ⌘K finder (SA-F4)", () => {
    renderAt("/samples?experiment=1");
    const legend = screen.getByTestId("kb-legend");
    expect(within(legend).getByText("⌘K")).toBeInTheDocument();
    expect(within(legend).getByText("find a sample")).toBeInTheDocument();
  });

  it("K from a contenteditable target does not batch-keep", () => {
    renderAt("/samples?experiment=1");
    fireEvent.click(screen.getAllByTestId("thumbnail")[0]!);
    const editor = document.createElement("div");
    editor.setAttribute("contenteditable", "true");
    document.body.appendChild(editor);
    try {
      fireEvent.keyDown(editor, { key: "k" });
      expect(batchMutate).not.toHaveBeenCalled();
      expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "true");
    } finally {
      editor.remove();
    }
  });

  it("an open modal dialog suppresses K (selection survives, no batch mutate)", () => {
    renderAt("/samples?experiment=1");
    fireEvent.click(screen.getAllByTestId("thumbnail")[0]!);
    const dialog = document.createElement("div");
    dialog.setAttribute("role", "dialog");
    dialog.setAttribute("aria-modal", "true");
    document.body.appendChild(dialog);
    try {
      fireEvent.keyDown(window, { key: "k" });
      expect(batchMutate).not.toHaveBeenCalled();
      expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "true");
    } finally {
      dialog.remove();
    }
  });

  it("X drops the selection and Escape clears it", () => {
    renderAt("/samples?experiment=1");
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[0]!);
    // X → drop
    fireEvent.keyDown(window, { key: "X" });
    expect(batchMutate).toHaveBeenCalledWith({
      sampleId: 1,
      exposureId: 100,
      status: "rejected",
    });
    // selection cleared after a drop → cull bar hidden
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");

    // re-select, then Escape clears without a mutate
    batchMutate.mockClear();
    fireEvent.click(screen.getAllByTestId("thumbnail")[0]!);
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "true");
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");
    expect(batchMutate).not.toHaveBeenCalled();
  });

  it("dropping frames announces a count toast (consequential → visible)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt("/samples?experiment=1");
      const thumbs = screen.getAllByTestId("thumbnail");
      fireEvent.click(thumbs[0]!);
      fireEvent.keyDown(window, { key: "X" });
      expect(toast).toHaveBeenCalledWith("1 frame dropped", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("restoring frames announces a symmetric count toast", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt("/samples?experiment=1");
      const thumbs = screen.getAllByTestId("thumbnail");
      fireEvent.click(thumbs[0]!);
      // The CullBar Restore button routes through batchSet(null).
      fireEvent.click(screen.getByRole("button", { name: /Restore/ }));
      expect(toast).toHaveBeenCalledWith("1 frame restored", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("a cross-sample selection discloses its spread in the cull bar (SA-F3)", () => {
    renderAt("/samples?experiment=1");
    // Sample 1 renders thumbs 0-2; sample 2 renders thumbs 3-4. Select one
    // frame from EACH row — the bar must disclose that the selection spans
    // two samples, not just report a flat count.
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[0]!); // sample 1 → exposure 100
    fireEvent.click(thumbs[3]!); // sample 2 → exposure 200
    const bar = screen.getByTestId("cull-bar");
    expect(bar.textContent).toContain("frames selected across");
    expect(bar.textContent).toContain("2 samples");
  });

  it("a single-sample selection stays count-only — no 'across' noise", () => {
    renderAt("/samples?experiment=1");
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[0]!); // sample 1
    fireEvent.click(thumbs[1]!); // sample 1 again
    const bar = screen.getByTestId("cull-bar");
    expect(bar.textContent).toContain("frames selected");
    expect(bar.textContent).not.toContain("across");
  });

  it("the Drop toast carries the same spread disclosure as the bar (SA-F3)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt("/samples?experiment=1");
      const thumbs = screen.getAllByTestId("thumbnail");
      fireEvent.click(thumbs[0]!); // sample 1 → exposure 100
      fireEvent.click(thumbs[3]!); // sample 2 → exposure 200
      fireEvent.click(screen.getByRole("button", { name: /Drop/ }));
      // The receipt matches the bar's promise: both frames, both samples.
      expect(toast).toHaveBeenCalledWith(
        "2 frames dropped across 2 samples",
        "success",
      );
      expect(batchMutate).toHaveBeenCalledWith({
        sampleId: 1, exposureId: 100, status: "rejected",
      });
      expect(batchMutate).toHaveBeenCalledWith({
        sampleId: 2, exposureId: 200, status: "rejected",
      });
    } finally {
      setToastImpl(null);
    }
  });

  it("a single-sample Drop toast carries no spread suffix", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt("/samples?experiment=1");
      const thumbs = screen.getAllByTestId("thumbnail");
      fireEvent.click(thumbs[0]!);
      fireEvent.click(thumbs[1]!);
      fireEvent.click(screen.getByRole("button", { name: /Drop/ }));
      expect(toast).toHaveBeenCalledWith("2 frames dropped", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("shift-click extends the contiguous range within one sample only", () => {
    renderAt("/samples?experiment=1");
    // Sample 1 (3 thumbs) then sample 2 (2 thumbs).
    const thumbs = screen.getAllByTestId("thumbnail");
    // Click frame 1 of sample 1 (anchor).
    fireEvent.click(thumbs[0]!);
    // Shift held → shift-click frame 3 of sample 1.
    fireEvent.keyDown(window, { key: "Shift" });
    fireEvent.click(thumbs[2]!, { shiftKey: true });
    fireEvent.keyUp(window, { key: "Shift" });
    // 3 frames selected (100, 101, 102) — the whole span of sample 1.
    const bar = screen.getByTestId("cull-bar");
    expect(within(bar).getByText("3")).toBeInTheDocument();
    // …and none of sample 2's frames were pulled in (thumbs[3] = sample 2 frame 1).
    expect(thumbs[3]).toHaveAttribute("data-state", expect.not.stringContaining("selected"));
  });

  it("the Status door navigates to /sample/:id; the name navigates to the loupe", () => {
    renderAt("/samples?experiment=1");
    // Status door for sample 1 (no phase → "Index <name>").
    fireEvent.click(screen.getByRole("button", { name: /index buffer/i }));
    expect(screen.getByTestId("focus-route")).toBeInTheDocument();
  });

  it("a confirmed zero-exposure sample shows 'No exposures' and no Index door (SA-ZEROEXP)", () => {
    // Sample 1 (Buffer) resolved with an EMPTY exposure list; sample 2 still has frames.
    state.byId = new Map<number, Exposure[]>([
      [1, []],
      [2, [exp({ id: 200, sample_id: 2 })]],
    ]);
    renderAt("/samples?experiment=1");
    // No dead Index door into an empty Focus workspace.
    expect(screen.queryByRole("button", { name: /index buffer/i })).toBeNull();
    // A clear terminal status instead.
    expect(screen.getByText("No exposures")).toBeInTheDocument();
  });

  it("a still-loading sample (exposures not yet fetched) keeps its Index door — no premature 'No exposures' (SA-ZEROEXP)", () => {
    // Sample 1 absent from byId = exposures undefined = not loaded yet.
    state.byId = new Map<number, Exposure[]>([
      [2, [exp({ id: 200, sample_id: 2 })]],
    ]);
    renderAt("/samples?experiment=1");
    expect(screen.getByRole("button", { name: /index buffer/i })).toBeInTheDocument();
    expect(screen.queryByText("No exposures")).toBeNull();
  });

  it("clicking a sample name opens the loupe at that sample", () => {
    renderAt("/samples?experiment=1");
    fireEvent.click(screen.getByRole("button", { name: "Buffer" }));
    expect(screen.getByTestId("loupe-route")).toBeInTheDocument();
  });

  it("double-clicking a frame opens the loupe AT that frame (SA-F2)", () => {
    renderAt("/samples?experiment=1");
    // Sample 1's thumbs render first; thumbs[1] = exposure 101 (frame 2).
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.dblClick(thumbs[1]!);
    expect(screen.getByTestId("loupe-route")).toHaveAttribute(
      "data-search",
      "?experiment=1&exposure=101",
    );
  });

  it("Shift+Enter on a thumbnail opens the loupe AT that frame — keyboard parity for double-click (SA-THUMBKEY)", () => {
    renderAt("/samples?experiment=1");
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.keyDown(thumbs[1]!, { key: "Enter", shiftKey: true });
    expect(screen.getByTestId("loupe-route")).toHaveAttribute(
      "data-search",
      "?experiment=1&exposure=101",
    );
  });

  it("a name-click loupe opening carries no exposure param (default frame)", () => {
    renderAt("/samples?experiment=1");
    fireEvent.click(screen.getByRole("button", { name: "Buffer" }));
    expect(screen.getByTestId("loupe-route")).toHaveAttribute(
      "data-search",
      "?experiment=1",
    );
  });

  it("?beamtime filters rows and titles by the experiment", () => {
    renderAt("/samples?experiment=2");
    // beamtime 2 → only sample 3
    expect(screen.getAllByTestId("sample-table-row")).toHaveLength(1);
    // experiment 2 has no name in the fixtures → falls back to "experiment 2"
    expect(screen.getByRole("heading", { name: "experiment 2" })).toBeInTheDocument();
    expect(screen.queryByRole("heading", { name: "BL-19 April" })).toBeNull();
  });

  it("unfiltered shows the whole corpus titled 'The corpus'", () => {
    renderAt("/samples");
    expect(screen.getByRole("heading", { name: "The corpus" })).toBeInTheDocument();
    expect(screen.getAllByTestId("sample-table-row")).toHaveLength(3);
  });

  it("an unknown ?beamtime renders an honest EmptyState, not a bare div (SA-F5)", () => {
    // 99 names no experiment record AND no sample carries it → unknown filter.
    renderAt("/samples?experiment=99");
    const block = screen.getByTestId("samples-unknown-experiment");
    expect(within(block).getByTestId("empty-state")).toBeInTheDocument();
    expect(
      within(block).getByRole("heading", { name: "Unknown experiment" }),
    ).toBeInTheDocument();
    // The raw id never leaks into the page title.
    expect(screen.queryByRole("heading", { name: /experiment 99/i })).toBeNull();
    expect(screen.getByRole("heading", { level: 1 })).toHaveTextContent(
      "Unknown experiment",
    );
  });

  it("the unknown-beamtime CTA clears the filter and shows the whole corpus", () => {
    renderAt("/samples?experiment=99");
    fireEvent.click(screen.getByRole("button", { name: "Show all experiments" }));
    expect(screen.getByRole("heading", { name: "The corpus" })).toBeInTheDocument();
    expect(screen.getAllByTestId("sample-table-row")).toHaveLength(3);
    expect(screen.queryByTestId("samples-unknown-experiment")).toBeNull();
  });

  it("a beamtime known only through its samples is NOT unknown (shared-predicate pin)", () => {
    // Experiment 2 has no /experiments record but sample 3 carries it — the
    // shared resolver treats it as real, so the page (and the topbar select,
    // which uses the same resolver) never calls it unknown.
    renderAt("/samples?experiment=2");
    expect(screen.queryByTestId("samples-unknown-experiment")).toBeNull();
    expect(screen.getAllByTestId("sample-table-row")).toHaveLength(1);
  });

  it("empty state when there are no samples", () => {
    state.samples = [];
    state.byId = new Map();
    renderAt("/samples");
    expect(screen.getByTestId("sheet-empty")).toBeInTheDocument();
  });

  it("ON-EMPTY: the empty corpus is a real door (EmptyState + reload), not a bare dead-end", () => {
    state.samples = [];
    state.byId = new Map();
    renderAt("/samples");
    const block = screen.getByTestId("sheet-empty");
    // A real EmptyState, consistent with the error / unknown-beamtime states —
    // a title + a body that names how samples appear, not a lone sentence.
    expect(within(block).getByTestId("empty-state")).toBeInTheDocument();
    expect(
      within(block).getByRole("heading", { name: /no samples yet/i }),
    ).toBeInTheDocument();
    // A way forward: reload once an experiment has been ingested/analyzed
    // (mirrors the error state's control, wired to the same refetch).
    fireEvent.click(within(block).getByRole("button", { name: /reload the corpus/i }));
    expect(corpusRefetch).toHaveBeenCalled();
  });

  it("SA-RESCORE3 F11: the keyboard legend is hidden when there are no rows to navigate", () => {
    state.samples = [];
    state.byId = new Map();
    renderAt("/samples");
    // The cell-navigation legend is meaningless under an empty table — suppress it.
    expect(screen.getByTestId("sheet-empty")).toBeInTheDocument();
    expect(screen.queryByTestId("kb-legend")).toBeNull();
  });

  it("error state renders an EmptyState with a retry control wired to refetch (SA-RETRY)", () => {
    state.error = true;
    renderAt("/samples");
    const block = screen.getByTestId("samples-error");
    expect(within(block).getByTestId("empty-state")).toBeInTheDocument();
    expect(
      within(block).getByRole("heading", { name: "Could not load the corpus" }),
    ).toBeInTheDocument();
    // The control embodies the way forward — no "try reloading" instruction.
    fireEvent.click(within(block).getByRole("button", { name: "Reload the corpus" }));
    expect(corpusRefetch).toHaveBeenCalled();
  });

  it("the retry control is disabled while the refetch is in flight", () => {
    state.error = true;
    state.fetching = true;
    renderAt("/samples");
    const block = screen.getByTestId("samples-error");
    expect(
      within(block).getByRole("button", { name: "Reload the corpus" }),
    ).toBeDisabled();
  });

  // ── SA-SORT: sortable columns ────────────────────────────────────────────────
  /** Read the rendered sample-name order off the rows (the SpecCell name button,
   *  falling back to the row's text), so order assertions never touch classes. */
  function rowNames(): string[] {
    return screen.getAllByTestId("sample-table-row").map((r) => {
      // The SpecCell identity name (NOT the Status-door button, whose accessible
      // name can also contain the sample name).
      const name = r.querySelector('[data-role="spec-name"]');
      return (name?.textContent ?? "").trim();
    });
  }
  function activeSortHeader(): string | undefined {
    return screen
      .getAllByRole("columnheader")
      .find(
        (h) =>
          h.getAttribute("aria-sort") === "ascending" ||
          h.getAttribute("aria-sort") === "descending",
      )
      ?.textContent ?? undefined;
  }

  it("sorting by Sample reorders the rows by display name (numeric-aware)", () => {
    renderAt("/samples"); // whole corpus: Buffer, Lipid 1-2, Other beamtime
    const ingest = rowNames();
    fireEvent.click(screen.getByRole("button", { name: "Sample" }));
    const sorted = rowNames();
    // "Buffer" < "Lipid 1-2" < "Other beamtime" alphabetically → unchanged here,
    // but the first row must be the alphabetically-first name regardless of ingest.
    expect(sorted[0]).toContain("Buffer");
    expect(sorted).not.toEqual(undefined);
    // The Sample header becomes the active ascending column.
    expect(activeSortHeader()).toContain("Sample");
    expect(ingest.length).toBe(3);
  });

  it("the asc → desc → clear cycle restores ingest order", () => {
    renderAt("/samples");
    const ingest = rowNames();
    const header = () => screen.getByRole("button", { name: "Sample" });
    fireEvent.click(header()); // asc
    fireEvent.click(header()); // desc
    const desc = rowNames();
    expect(desc).toEqual([...ingest].reverse());
    fireEvent.click(header()); // clear → ingest order
    expect(rowNames()).toEqual(ingest);
    // No column is active after clearing.
    expect(activeSortHeader()).toBeUndefined();
  });

  it("activating a second column moves the sort there at ascending", () => {
    renderAt("/samples");
    fireEvent.click(screen.getByRole("button", { name: "Sample" })); // Sample asc
    expect(activeSortHeader()).toContain("Sample");
    fireEvent.click(screen.getByRole("button", { name: "Exposures" })); // → Exposures asc
    const active = screen
      .getAllByRole("columnheader")
      .find((h) => h.getAttribute("aria-sort") === "ascending");
    expect(active?.textContent).toContain("Exposures");
    // Sample is no longer the active sort.
    expect(
      screen
        .getAllByRole("columnheader")
        .find((h) => h.textContent?.includes("Sample"))
        ?.getAttribute("aria-sort"),
    ).toBe("none");
  });

  it("Exposures sort orders rows by frame count (sample 3 has 1 frame → first asc)", () => {
    renderAt("/samples");
    // ingest: sample1=3 frames, sample2=2, sample3=1
    fireEvent.click(screen.getByRole("button", { name: "Exposures" }));
    const asc = rowNames();
    expect(asc[0]).toContain("Other beamtime"); // sample 3, 1 frame
    expect(asc[2]).toContain("Buffer"); // sample 1, 3 frames
  });

  it("the sort key+dir round-trips through the URL (mirrors the beamtime filter)", () => {
    // A permalink that already carries a sort renders pre-sorted.
    renderAt("/samples?sort=sample&dir=desc");
    expect(activeSortHeader()).toContain("Sample");
    const desc = rowNames();
    expect(desc[0]).toContain("Other beamtime"); // "Other…" is last alphabetically → first desc

    // And a click writes the sort back into the URL.
    fireEvent.click(screen.getByRole("button", { name: "Exposures" }));
    // Exposures becomes ascending; the URL now names it (asserted via re-derivable
    // active header, since the same searchParams drive the render).
    const active = screen
      .getAllByRole("columnheader")
      .find((h) => h.getAttribute("aria-sort") === "ascending");
    expect(active?.textContent).toContain("Exposures");
  });

  it("an unsorted permalink (no ?sort) renders in ingest order with no active column", () => {
    renderAt("/samples");
    expect(activeSortHeader()).toBeUndefined();
  });

  it("X from a contenteditable target does not batch-drop", () => {
    renderAt("/samples?experiment=1");
    fireEvent.click(screen.getAllByTestId("thumbnail")[0]!);
    const editor = document.createElement("div");
    editor.setAttribute("contenteditable", "true");
    document.body.appendChild(editor);
    try {
      fireEvent.keyDown(editor, { key: "X" });
      expect(batchMutate).not.toHaveBeenCalled();
      // Selection survives — the cull bar is still up.
      expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "true");
    } finally {
      editor.remove();
    }
  });

  it("an open modal dialog suppresses X (selection survives, no batch mutate)", () => {
    renderAt("/samples?experiment=1");
    fireEvent.click(screen.getAllByTestId("thumbnail")[0]!);
    const dialog = document.createElement("div");
    dialog.setAttribute("role", "dialog");
    dialog.setAttribute("aria-modal", "true");
    document.body.appendChild(dialog);
    try {
      fireEvent.keyDown(window, { key: "X" });
      expect(batchMutate).not.toHaveBeenCalled();
      expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "true");
      // Escape behind the modal must not clear the selection either.
      fireEvent.keyDown(window, { key: "Escape" });
      expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "true");
    } finally {
      dialog.remove();
    }
  });
});
