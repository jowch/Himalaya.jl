import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within, act } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import type {
  CorpusSample, Exposure, IndexEntry, Peak, Trace, Assignment,
} from "../../src/api";

// ── mutator spies ─────────────────────────────────────────────────────────────
const addPeakMutate = vi.fn();
const removePeakMutate = vi.fn();
const setPeakExclMutate = vi.fn();
const addAssignMutate = vi.fn();
const removeAssignMutate = vi.fn();
const setAssignStateMutate = vi.fn();
const commitCustomMutate = vi.fn();

// ── mock data plane (mutated per test) ────────────────────────────────────────
const state = {
  corpus: [] as CorpusSample[],
  samples: [] as CorpusSample[],
  exposures: [] as Exposure[],
  trace: undefined as Trace | undefined,
  peaks: [] as Peak[],
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
  loading: false,
  activeSampleId: undefined as number | undefined,
  activeExposureId: undefined as number | undefined,
};

// Capture props the q-link triple receives so we can assert the shared wire
// without driving a brittle JSDOM pointer round-trip.
const combsProps = { hoveredQ: undefined as number | undefined, onHoverQ: undefined as unknown };

vi.mock("../../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.corpus, isLoading: state.loading }),
  useSamples: () => ({ data: state.samples, isLoading: state.loading }),
  useExperiment: () => ({ data: { id: 1, name: "BL-19", q_units: "A-1" } }),
  useExposures: () => ({ data: state.exposures, isLoading: state.loading }),
  useExposure: () => ({ data: state.exposures.find((e) => e.id === state.activeExposureId) }),
  useTrace: () => ({ data: state.trace, isLoading: state.loading }),
  usePeaks: () => ({ data: state.peaks, isLoading: state.loading }),
  useIndices: () => ({ data: state.indices, isLoading: state.loading }),
  useAssignment: () => ({ data: state.assignment, isLoading: state.loading }),
  useAddPeak: () => ({ mutate: addPeakMutate }),
  useRemovePeak: () => ({ mutate: removePeakMutate }),
  useSetPeakExcluded: () => ({ mutate: setPeakExclMutate }),
  useAddAssignmentPhase: () => ({ mutate: addAssignMutate }),
  useRemoveAssignmentPhase: () => ({ mutate: removeAssignMutate }),
  useSetAssignmentState: () => ({ mutate: setAssignStateMutate }),
  useCommitCustomIndex: () => ({ mutate: commitCustomMutate }),
}));

// useSyncActiveSampleFromRoute runs REAL here (against the mocked queries and
// store): the page consumes its route-resolution status, so shimming it would
// hide the mid-session bogus-URL behaviour (F-STALEURL).
vi.mock("../../src/hooks/useAutoPickExposure", () => ({
  useAutoPickExposure: () => {},
}));

// Zustand store: return the mock active ids.
vi.mock("../../src/state", () => ({
  useAppState: (sel: (s: Record<string, unknown>) => unknown) =>
    sel({
      activeSampleId: state.activeSampleId,
      activeExposureId: state.activeExposureId,
      setActiveSample: (id: number | undefined) => {
        state.activeSampleId = id;
      },
    }),
}));

// boneyard Skeleton: render children when not loading.
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

// DetectorImage touches fetch / createImageBitmap (absent in JSDOM).
vi.mock("../../src/print/detector", async (orig) => {
  const real = (await orig()) as Record<string, unknown>;
  return { ...real, DetectorImage: () => <div data-testid="mock-detector-image" /> };
});
vi.mock("../../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

// Capture the CombsPanel props (q-link wire) while keeping its render trivial.
vi.mock("../../src/print/components/CombsPanel", () => ({
  CombsPanel: (p: { hoveredQ?: number; onHoverQ?: unknown }) => {
    combsProps.hoveredQ = p.hoveredQ;
    combsProps.onHoverQ = p.onHoverQ;
    return <div data-testid="combs-panel" />;
  },
}));

import { FocusPage } from "../../src/print/pages/FocusPage";
import { setAnnounceImpl } from "../../src/lib/announce";
import { setToastImpl } from "../../src/lib/toast";

// ── fixtures ──────────────────────────────────────────────────────────────────
function corpus(over: Partial<CorpusSample> = {}): CorpusSample {
  return {
    id: 42, experiment_id: 1, name: "JC042", display_name: "JC042 — LL37",
    notes: null, q_units: "A-1", tags: [], ...over,
  } as CorpusSample;
}
function exp(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 7, sample_id: 42, filename: "JC042-001.dat", kind: "file",
    selected: true, status: "accepted", image_path: "/x.tif", image_version: "v1",
    tags: [], sources: [], trace_hash: "h", analysis_inputs_hash: "hashA", ...over,
  };
}
function ix(over: Partial<IndexEntry> = {}): IndexEntry {
  return {
    id: 1, exposure_id: 7, phase: "Pn3m", basis: 0.15, score: 0.9,
    r_squared: 0.99, lattice_d: 197, ngc: -1.5, status: "candidate",
    kind: "auto", inputs_hash: "hashA",
    peaks: [{ peak_id: 1, ratio_position: 1, residual: 0.001, q_observed: 0.2 }],
    predicted_q: [0.2, 0.28], ...over,
  };
}

function seedFull(): void {
  state.activeSampleId = 42;
  state.activeExposureId = 7;
  state.corpus = [corpus()];
  state.samples = [corpus({ notes: "watch q ≈ 0.200" })];
  state.exposures = [exp()];
  state.trace = { q: [0.1, 0.2, 0.3], I: [10, 40, 20], sigma: [1, 1, 1] };
  state.peaks = [
    { id: 1, exposure_id: 7, q: 0.2, intensity: 40, prominence: 10, sharpness: 2, source: "auto", excluded: false },
    { id: 2, exposure_id: 7, q: 0.3, intensity: 20, prominence: 8, sharpness: 1.5, source: "auto", excluded: false },
  ];
  state.indices = [
    ix({ id: 1, phase: "Pn3m" }),
    ix({ id: 2, phase: "Lamellar", score: 0.5, predicted_q: [0.3, 0.6],
      peaks: [{ peak_id: 2, ratio_position: 1, residual: 0.001, q_observed: 0.3 }] }),
  ];
  state.assignment = { exposure_id: 7, state: "indexed", members: [1] };
  state.loading = false;
}

function renderAt(sampleId: number | string) {
  return render(
    <MemoryRouter initialEntries={[`/sample/${sampleId}`]}>
      <Routes>
        <Route path="/sample/:sampleId" element={<FocusPage />} />
        <Route path="/samples" element={<div data-testid="sheet">sheet</div>} />
      </Routes>
    </MemoryRouter>,
  );
}

beforeEach(() => {
  vi.clearAllMocks();
  combsProps.hoveredQ = undefined;
  combsProps.onHoverQ = undefined;
  seedFull();
});

describe("FocusPage", () => {
  it("renders the trace plate, assignment rail, detector panel and combs panel", () => {
    renderAt(42);
    expect(screen.getByTestId("trace-plate")).toBeInTheDocument();
    expect(screen.getByTestId("assignment-rail")).toBeInTheDocument();
    expect(screen.getByTestId("detector-panel")).toBeInTheDocument();
    expect(screen.getByTestId("combs-panel")).toBeInTheDocument();
  });

  it("toggling a candidate row fires useAddAssignmentPhase().mutate(indexId)", () => {
    renderAt(42);
    // Lamellar (id 2) is NOT in the active set → its CandidateRow toggles it on.
    const candidate = screen.getByRole("button", { name: /Lamellar/ });
    fireEvent.click(candidate);
    expect(addAssignMutate).toHaveBeenCalledWith(2);
  });

  it("clicking an in-call candidate removes it from the assignment", () => {
    renderAt(42);
    const pn3m = screen.getByRole("button", { name: /Pn3m, in assignment/ });
    fireEvent.click(pn3m);
    expect(removeAssignMutate).toHaveBeenCalledWith(1);
  });

  it("does not render the legacy '+ Add speculative' affordance (custom-index is the hypothesis tool)", () => {
    renderAt(42);
    expect(screen.queryByTestId("add-speculative-button")).toBeNull();
    expect(screen.queryByText("+ Add speculative")).toBeNull();
  });

  it("arming '+ Peak' then clicking the trace fires useAddPeak().mutate(q)", () => {
    const { container } = renderAt(42);
    const addPeakBtn = screen.getByText("+ Peak");
    fireEvent.click(addPeakBtn);
    expect(addPeakBtn).toHaveAttribute("aria-pressed", "true");
    // The plot's add path: click empty plot space → interaction.onAddPeak(q).
    const svg = container.querySelector('svg[data-testid="trace-plate-plot"]')!;
    fireEvent.click(svg, { clientX: 300, clientY: 150 });
    expect(addPeakMutate).toHaveBeenCalledTimes(1);
    expect(typeof addPeakMutate.mock.calls[0]![0]).toBe("number");
  });

  it("wires the same q-link (hoveredQ/onHoverQ) into the combs panel and trace", () => {
    const { container } = renderAt(42);
    // The CombsPanel mock captured the shared onHoverQ; firing it sets hoveredQ
    // on the page, which then propagates back through every panel.
    expect(typeof combsProps.onHoverQ).toBe("function");
    act(() => {
      (combsProps.onHoverQ as (q?: number) => void)(0.2);
    });
    // The trace plate's plot reflects the shared hoveredQ → q-readout chip lights.
    const readout = container.querySelector('[data-role="q-readout"]');
    expect(readout).toBeTruthy();
    // …and the same hoveredQ reaches the combs panel (the third q-link surface).
    expect(combsProps.hoveredQ).toBe(0.2);
  });

  it("captions the detector rings with the assignment's phase for the displayed frame (FO-RING)", () => {
    // assignment members = [1] → Pn3m is the active index; the detector caption
    // must name it in text (the second channel for the hue-only rings) and tie
    // it to the displayed frame.
    renderAt(42);
    const panel = screen.getByTestId("detector-panel");
    const caption = within(panel).getByTestId("detector-ring-caption");
    expect(within(caption).getByTestId("phase-chip")).toHaveTextContent("Pn3m");
    // label leads the reading order
    expect(caption.textContent).toMatch(/^rings:/);
    expect(caption).toHaveTextContent("this frame's indexing");
  });

  it("does not chip a fully-landed custom index that colours zero rings (FO-RING no-lie)", () => {
    // A committed custom index arrives with peaks: [] (insert_custom_index!
    // writes no index_peaks rows). With every predicted_q within tol of an
    // observed peak it emits NO coloured or ghost ring, so the caption must
    // not name it — only the ring-emitting sibling.
    state.indices = [
      ix({ id: 1, phase: "Pn3m" }),
      ix({
        id: 3, phase: "Hexagonal", kind: "speculative",
        peaks: [],
        predicted_q: [0.2, 0.3], // both sit on observed peaks -> zero rings
      }),
    ];
    state.assignment = { exposure_id: 7, state: "indexed", members: [1, 3] };
    renderAt(42);
    const panel = screen.getByTestId("detector-panel");
    const caption = within(panel).getByTestId("detector-ring-caption");
    const chips = within(caption).getAllByTestId("phase-chip");
    expect(chips.map((c) => c.textContent)).toEqual(["Pn3m"]);
  });

  it("omits the detector ring caption when the frame has no assigned phases (FO-RING honest empty)", () => {
    state.assignment = { exposure_id: 7, state: "indexed", members: [] };
    renderAt(42);
    const panel = screen.getByTestId("detector-panel");
    expect(within(panel).queryByTestId("detector-ring-caption")).toBeNull();
  });

  it("'+ custom index…' opens the CustomIndexModal", () => {
    renderAt(42);
    fireEvent.click(screen.getByTestId("custom-index-trigger"));
    expect(screen.getByTestId("custom-index-modal")).toBeInTheDocument();
  });

  it("toggling a candidate phase announces SR-only (frequent → quiet channel)", () => {
    const announce = vi.fn();
    const toast = vi.fn();
    setAnnounceImpl(announce);
    setToastImpl(toast);
    try {
      renderAt(42);
      fireEvent.click(screen.getByRole("button", { name: /Lamellar/ }));
      expect(announce.mock.calls[0]?.[0]).toBe("Lamellar added to the call");
      // a candidate toggle is NOT a visible toast (would be spam)
      expect(toast).not.toHaveBeenCalled();
    } finally {
      setAnnounceImpl(null);
      setToastImpl(null);
    }
  });

  it("adding a peak announces 'Peak added' SR-only", () => {
    const announce = vi.fn();
    setAnnounceImpl(announce);
    try {
      const { container } = renderAt(42);
      fireEvent.click(screen.getByText("+ Peak"));
      const svg = container.querySelector('svg[data-testid="trace-plate-plot"]')!;
      fireEvent.click(svg, { clientX: 300, clientY: 150 });
      expect(announce.mock.calls[0]?.[0]).toBe("Peak added");
    } finally {
      setAnnounceImpl(null);
    }
  });

  it("arming '+ Peak' then adding via the q field fires useAddPeak().mutate(q) and announces", () => {
    const announce = vi.fn();
    setAnnounceImpl(announce);
    try {
      renderAt(42);
      fireEvent.click(screen.getByText("+ Peak"));
      // The keyboard parity for click-empty-space-adds: type a q, press Add.
      const input = screen.getByLabelText("q value for new peak");
      fireEvent.change(input, { target: { value: "0.15" } });
      fireEvent.click(screen.getByRole("button", { name: "Add peak at q" }));
      expect(addPeakMutate).toHaveBeenCalledTimes(1);
      expect(addPeakMutate).toHaveBeenCalledWith(0.15);
      expect(announce.mock.calls[0]?.[0]).toBe("Peak added");
    } finally {
      setAnnounceImpl(null);
    }
  });

  it("arming '+ Peak' then pressing Enter on a focused peak mark removes it and announces", () => {
    const announce = vi.fn();
    setAnnounceImpl(announce);
    try {
      const { container } = renderAt(42);
      fireEvent.click(screen.getByText("+ Peak"));
      const mark = container.querySelector('[data-role="plot-peaks"] [role="button"]')!;
      expect(mark).toBeTruthy();
      fireEvent.keyDown(mark, { key: "Enter" });
      expect(removePeakMutate).toHaveBeenCalledWith(1);
      expect(announce.mock.calls[0]?.[0]).toBe("Peak removed");
      fireEvent.keyDown(mark, { key: "Enter", altKey: true });
      expect(setPeakExclMutate).toHaveBeenCalledWith({ peakId: 1, excluded: true });
    } finally {
      setAnnounceImpl(null);
    }
  });

  it("committing a custom index announces a visible toast (consequential)", () => {
    const toast = vi.fn();
    setToastImpl(toast);
    try {
      renderAt(42);
      fireEvent.click(screen.getByTestId("custom-index-trigger"));
      // The modal's Add commits the custom index.
      const modal = screen.getByTestId("custom-index-modal");
      fireEvent.click(within(modal).getByRole("button", { name: /^Add/ }));
      expect(commitCustomMutate).toHaveBeenCalled();
      expect(toast).toHaveBeenCalledWith(expect.stringContaining("index added"), "success");
    } finally {
      setToastImpl(null);
    }
  });

  it("Escape sequence: an open custom-index modal closes first; only the NEXT Escape disarms '+ Peak' (F7)", () => {
    renderAt(42);
    const addPeakBtn = screen.getByText("+ Peak");
    fireEvent.click(addPeakBtn);
    expect(addPeakBtn).toHaveAttribute("aria-pressed", "true");
    fireEvent.click(screen.getByTestId("custom-index-trigger"));
    expect(screen.getByTestId("custom-index-modal")).toBeInTheDocument();
    // Escape #1: the dialog owns Escape — the modal closes and the armed mode
    // survives (the dialog is still in the DOM during that dispatch, so the
    // TracePlate guard sees it).
    fireEvent.keyDown(document.body, { key: "Escape" });
    expect(screen.queryByTestId("custom-index-modal")).toBeNull();
    expect(addPeakBtn).toHaveAttribute("aria-pressed", "true");
    // Escape #2: no dialog left — this one disarms.
    fireEvent.keyDown(document.body, { key: "Escape" });
    // FO-RESCORE2 F12: the disarmed toggle reads aria-pressed="false" (not
    // dropped) — a toggle button keeps the attribute in both states.
    expect(addPeakBtn).toHaveAttribute("aria-pressed", "false");
  });

  it("not-found: no corpus sample for the id", () => {
    state.corpus = [];
    state.activeSampleId = undefined;
    renderAt(999);
    expect(screen.getByTestId("focus-not-found")).toBeInTheDocument();
  });

  it("not-found renders an EmptyState whose action leads back to the contact sheet (FO-ERR)", () => {
    state.corpus = [];
    state.activeSampleId = undefined;
    renderAt(999);
    const block = screen.getByTestId("focus-not-found");
    expect(within(block).getByTestId("empty-state")).toBeInTheDocument();
    expect(
      within(block).getByRole("heading", { name: "Sample not found" }),
    ).toBeInTheDocument();
    fireEvent.click(
      within(block).getByRole("button", { name: "Back to the contact sheet" }),
    );
    expect(screen.getByTestId("sheet")).toBeInTheDocument();
  });

  it("not-found exposes its sole heading as an h1, not a level-skipped h2 (FO-NFHEAD)", () => {
    // The shell carries no h1 and the not-found branch returns before any
    // TracePlate, so EmptyState's title is the document's only heading — it must
    // be the top level (WCAG 1.3.1 / no missing h1).
    state.corpus = [];
    state.activeSampleId = undefined;
    renderAt(999);
    expect(screen.getByRole("heading", { level: 1, name: "Sample not found" })).toBeInTheDocument();
    expect(screen.queryByRole("heading", { level: 2, name: "Sample not found" })).toBeNull();
  });

  it("mid-session: a bogus numeric /sample/:id renders not-found, not the previous sample (F-STALEURL)", () => {
    // seedFull leaves activeSampleId at 42 with a warm corpus cache; the URL
    // names a sample that does not exist. Rendering sample 42 under
    // /sample/99999 would be a lie: the page must show not-found instead.
    renderAt(99999);
    expect(screen.getByTestId("focus-not-found")).toBeInTheDocument();
    expect(screen.queryByTestId("trace-plate")).toBeNull();
  });

  it("mid-session: a non-numeric /sample/:id renders not-found (F-STALEURL)", () => {
    renderAt("not-a-number");
    expect(screen.getByTestId("focus-not-found")).toBeInTheDocument();
    expect(screen.queryByTestId("trace-plate")).toBeNull();
  });

  it("no-exposure: the sample has no exposures", () => {
    state.exposures = [];
    state.activeExposureId = undefined;
    state.trace = undefined;
    renderAt(42);
    expect(screen.getByText(/no exposures/i)).toBeInTheDocument();
  });

  it("does not render a stale-index banner (peak edits auto-reanalyze server-side)", () => {
    // A speculative index keeps its inputs_hash across reanalysis, so the old
    // hash-mismatch banner was a permanent false alarm; it was removed. Even
    // with a mismatched hash, no alert renders.
    state.indices = [ix({ id: 1, inputs_hash: "OLD" })];
    state.assignment = { exposure_id: 7, state: "indexed", members: [1] };
    renderAt(42);
    expect(screen.queryByRole("alert")).toBeNull();
    expect(screen.queryByText(/stale/i)).toBeNull();
  });

  it("renders the export-button in the trace plate when the sample has trace data", () => {
    renderAt(42);
    const plate = screen.getByTestId("trace-plate");
    expect(within(plate).getByTestId("export-button")).toBeInTheDocument();
    // ariaContext is wired — the Copy button carries the aria-label
    expect(
      within(plate).getByRole("button", { name: /copy trace plot to clipboard/i }),
    ).toBeInTheDocument();
  });

  it("export stays enabled with a rendered trace and zero peaks (F4: WYSIWYG figure)", () => {
    // A peakless trace still renders a full figure; gating export on peaks
    // made the buttons lie. Only trace emptiness disables export. The menu
    // trigger is the probe: its disabled is exactly the page gate (Copy is
    // environment-disabled in JSDOM — no ClipboardItem — so it can't see it).
    state.peaks = [];
    renderAt(42);
    const plate = screen.getByTestId("trace-plate");
    expect(within(plate).getByTestId("export-menu-trigger")).not.toBeDisabled();
  });

  it("export is disabled when the trace is empty or missing (F4 data gate)", () => {
    state.trace = { q: [], I: [], sigma: [] };
    const view = renderAt(42);
    let plate = screen.getByTestId("trace-plate");
    expect(within(plate).getByTestId("export-menu-trigger")).toBeDisabled();

    view.unmount();
    state.trace = undefined;
    renderAt(42);
    plate = screen.getByTestId("trace-plate");
    expect(within(plate).getByTestId("export-menu-trigger")).toBeDisabled();
  });

  it("zero peaks: cart and candidate list agree that peaks come first (F3)", () => {
    state.peaks = [];
    state.indices = [];
    state.assignment = { exposure_id: 7, state: "indexed", members: [] };
    renderAt(42);
    const cartEmpty = screen.getByTestId("assignment-empty");
    expect(cartEmpty).toHaveTextContent(
      "No peaks marked. Find peaks on the trace to start indexing.",
    );
    // The default copy's two lies must be gone with zero peaks.
    expect(screen.queryByText(/Every peak is unindexed/)).toBeNull();
    expect(screen.queryByText(/Check a candidate below/)).toBeNull();
    // The candidate list derives from the SAME peaks-empty predicate.
    expect(
      screen.getByText("Candidates appear once peaks are marked."),
    ).toBeInTheDocument();
    expect(screen.queryByText("No candidate indexings.")).toBeNull();
  });

  it("peaks but zero candidates: cart points at the custom-index tool (F3)", () => {
    state.indices = [];
    state.assignment = { exposure_id: 7, state: "indexed", members: [] };
    renderAt(42);
    expect(screen.getByTestId("assignment-empty")).toHaveTextContent(
      "No phase assigned. No candidate fits these peaks. Try a custom index.",
    );
    // The named next action renders directly below in the cart footer.
    expect(screen.getByTestId("custom-index-trigger")).toBeInTheDocument();
    // Candidate-list line for the peaks-exist branch is unchanged.
    expect(screen.getByText("No candidate indexings.")).toBeInTheDocument();
  });

  it("peaks all excluded: cart names the exclusion, not a failed candidate search (FO-ALLEXCLUDED-CAPTION)", () => {
    state.peaks = [
      { id: 1, exposure_id: 7, q: 0.2, intensity: 40, prominence: 10, sharpness: 2, source: "auto", excluded: true },
      { id: 2, exposure_id: 7, q: 0.3, intensity: 20, prominence: 8, sharpness: 1.5, source: "auto", excluded: true },
    ];
    state.indices = [];
    state.assignment = { exposure_id: 7, state: "indexed", members: [] };
    renderAt(42);
    // Distinct, honest copy — peaks exist but none are indexable.
    expect(screen.getByTestId("assignment-empty")).toHaveTextContent(
      "All peaks are excluded. Restore a peak, or add one, to index.",
    );
    // NOT the no-candidate-fits message (which implies candidates were tried and failed).
    expect(screen.queryByText(/No candidate fits these peaks/)).toBeNull();
    // The candidate list agrees with the same reason (they must never contradict).
    expect(screen.getByText("Candidates appear once a peak is restored.")).toBeInTheDocument();
    expect(screen.queryByText("No candidate indexings.")).toBeNull();
  });

  it("candidates exist: the cart keeps its default empty copy (F3 unchanged branch)", () => {
    state.assignment = { exposure_id: 7, state: "indexed", members: [] };
    renderAt(42);
    expect(screen.getByTestId("assignment-empty")).toHaveTextContent(
      "No phase assigned. Every peak is unindexed. Check a candidate below.",
    );
  });
});
