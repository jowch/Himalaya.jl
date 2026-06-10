import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
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

// ── fixtures ──────────────────────────────────────────────────────────────────
function corpus(over: Partial<CorpusSample> = {}): CorpusSample {
  return {
    id: 1, experiment_id: 1, name: "JC001", display_name: "JC001 — buffer",
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
    corpus({ id: 1, experiment_id: 1, name: "JC001", display_name: "Buffer" }),
    corpus({ id: 2, experiment_id: 1, name: "JC002", display_name: "Lipid 1-2", phase: "Pn3m" }),
    corpus({ id: 3, experiment_id: 2, name: "JC003", display_name: "Other beamtime" }),
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

function renderAt(path = "/samples") {
  return render(
    <MemoryRouter initialEntries={[path]}>
      <Routes>
        <Route path="/samples" element={<SamplesPage />} />
        <Route path="/samples/loupe/:sampleId" element={<div data-testid="loupe-route" />} />
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
    renderAt("/samples?beamtime=1");
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
    // beamtime title = experiment name
    expect(screen.getByRole("heading", { name: "BL-19 April" })).toBeInTheDocument();
    expect(screen.getByTestId("sheet-table")).toBeInTheDocument();
    expect(screen.getByTestId("kb-legend")).toBeInTheDocument();
  });

  it("renders one row per filtered sample", () => {
    renderAt("/samples?beamtime=1");
    // beamtime 1 has 2 samples (ids 1, 2)
    expect(screen.getAllByTestId("sample-table-row")).toHaveLength(2);
  });

  it("clicking a thumb selects it → cull bar shows with the right count", () => {
    renderAt("/samples?beamtime=1");
    const thumbs = screen.getAllByTestId("thumbnail");
    fireEvent.click(thumbs[0]!);
    const bar = screen.getByTestId("cull-bar");
    expect(bar).toHaveAttribute("data-show", "true");
    expect(within(bar).getByText("1")).toBeInTheDocument();
  });

  it("CullBar Drop calls the batch mutator with the selected exposure", () => {
    renderAt("/samples?beamtime=1");
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

  it("X drops the selection and Escape clears it", () => {
    renderAt("/samples?beamtime=1");
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
      renderAt("/samples?beamtime=1");
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
      renderAt("/samples?beamtime=1");
      const thumbs = screen.getAllByTestId("thumbnail");
      fireEvent.click(thumbs[0]!);
      // The CullBar Restore button routes through batchSet(null).
      fireEvent.click(screen.getByRole("button", { name: /Restore/ }));
      expect(toast).toHaveBeenCalledWith("1 frame restored", "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("shift-click extends the contiguous range within one sample only", () => {
    renderAt("/samples?beamtime=1");
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
    renderAt("/samples?beamtime=1");
    // Status door for sample 1 (no phase → "Index <name>").
    fireEvent.click(screen.getByRole("button", { name: /index buffer/i }));
    expect(screen.getByTestId("focus-route")).toBeInTheDocument();
  });

  it("clicking a sample name opens the loupe at that sample", () => {
    renderAt("/samples?beamtime=1");
    fireEvent.click(screen.getByRole("button", { name: "Buffer" }));
    expect(screen.getByTestId("loupe-route")).toBeInTheDocument();
  });

  it("?beamtime filters rows and titles by the experiment", () => {
    renderAt("/samples?beamtime=2");
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

  it("empty state when there are no samples", () => {
    state.samples = [];
    state.byId = new Map();
    renderAt("/samples");
    expect(screen.getByTestId("sheet-empty")).toBeInTheDocument();
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

  it("X from a contenteditable target does not batch-drop", () => {
    renderAt("/samples?beamtime=1");
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
    renderAt("/samples?beamtime=1");
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
