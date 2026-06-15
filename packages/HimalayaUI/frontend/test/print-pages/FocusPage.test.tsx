import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within, act } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
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
  exposures: [] as Exposure[] | undefined,
  trace: undefined as Trace | undefined,
  peaks: [] as Peak[],
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
  loading: false,
  activeSampleId: undefined as number | undefined,
  activeExposureId: undefined as number | undefined,
  // sibling order for the [ ] sample-step shortcut (useExperimentSiblings mock)
  sibPrev: undefined as { id: number } | undefined,
  sibNext: undefined as { id: number } | undefined,
};

// Capture props the q-link triple receives so we can assert the shared wire
// without driving a brittle JSDOM pointer round-trip.
const combsProps = { hoveredQ: undefined as number | undefined, onHoverQ: undefined as unknown };

vi.mock("../../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.corpus, isLoading: state.loading }),
  useSamples: () => ({ data: state.samples, isLoading: state.loading }),
  useExperiment: () => ({ data: { id: 1, name: "BL-19", q_units: "A-1" } }),
  useExposures: () => ({ data: state.exposures, isLoading: state.loading }),
  useExposure: () => ({ data: state.exposures?.find((e) => e.id === state.activeExposureId) }),
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
// Stub only the effect hook (it touches queries + the store); keep the real
// `acceptableExposures` so the exposure-axis stepper filters rejected frames
// against the same predicate production uses (FO-EXPSKIP).
vi.mock("../../src/hooks/useAutoPickExposure", async (importOriginal) => ({
  ...(await importOriginal<typeof import("../../src/hooks/useAutoPickExposure")>()),
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
      setActiveExposure: (id: number | undefined) => {
        state.activeExposureId = id;
      },
    }),
}));

// useExperimentSiblings drives the [ ] sample-step shortcut; return the seeded
// prev/next so the keyboard test can assert navigation without a real corpus order.
vi.mock("../../src/hooks/useExperimentSiblings", () => ({
  useExperimentSiblings: () => ({
    activeSample: state.activeSampleId !== undefined ? { id: state.activeSampleId } : undefined,
    siblings: [],
    index: 0,
    prev: state.sibPrev,
    next: state.sibNext,
  }),
}));

// boneyard Skeleton: surface the `loading` gate as a data attribute so tests
// can assert WHEN the skeleton holds (the real component swaps children for
// bones; here we keep children mounted and just expose the flag). data-loading
// lets the FO-NAV-SKELETON regression check the exposure-resolution window.
vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children, loading }: { children: React.ReactNode; loading?: boolean }) => (
    <div data-testid="focus-skeleton-gate" data-loading={loading ? "true" : "false"}>
      {children}
    </div>
  ),
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
  state.sibPrev = undefined;
  state.sibNext = undefined;
}

function LocationProbe() {
  const loc = useLocation();
  return <div data-testid="loc" data-path={loc.pathname} />;
}

function focusTreeAt(sampleId: number | string) {
  return (
    <MemoryRouter initialEntries={[`/sample/${sampleId}`]}>
      <LocationProbe />
      <Routes>
        <Route path="/sample/:sampleId" element={<FocusPage />} />
        <Route path="/samples" element={<div data-testid="sheet">sheet</div>} />
      </Routes>
    </MemoryRouter>
  );
}

function renderAt(sampleId: number | string) {
  return render(focusTreeAt(sampleId));
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

  it("the candidate-rail note is the distilled one-sentence guide (DI-FOCUSNOTE)", () => {
    renderAt(42);
    // Distilled to the single load-bearing fact (a sample can be multiphasic),
    // dropping the dense second swap/coexist sentence.
    expect(
      screen.getByText(/A sample can be multiphasic, so check every phase that fits\./),
    ).toBeInTheDocument();
    // the old run-on second sentence is gone
    expect(screen.queryByText(/Candidates that explain the same peaks swap/)).toBeNull();
  });

  // ── keyboard: the two-axis model (shared shortcut library) ───────────────────
  describe("keyboard — two-axis model", () => {
    it("] / [ step the sample (agree with the stepper via useExperimentSiblings)", () => {
      state.sibNext = { id: 43 };
      state.sibPrev = { id: 41 };
      renderAt(42);
      fireEvent.keyDown(document.body, { key: "]" });
      expect(screen.getByTestId("loc")).toHaveAttribute("data-path", "/sample/43");
    });

    it("[ steps to the previous sample", () => {
      state.sibPrev = { id: 41 };
      renderAt(42);
      fireEvent.keyDown(document.body, { key: "[" });
      expect(screen.getByTestId("loc")).toHaveAttribute("data-path", "/sample/41");
    });

    // FO-NAV-STATE: React Router does NOT remount FocusPage on a same-route
    // /sample/:id step, so page-owned interaction state (the "+ Peak" arm, the
    // zoom window) would survive a sample switch — the first click on the next
    // sample's trace would silently mutate ITS peaks. The arm must reset on the
    // sample change.
    it("resets per-sample interaction state (arm + candidate preview) when the active sample changes (FO-NAV-STATE)", () => {
      seedFull();
      state.corpus = [corpus(), corpus({ id: 43, name: "JC043" })];
      const view = renderAt(42);
      // Arm "+ Peak" AND preview a candidate.
      fireEvent.click(screen.getByText("+ Peak"));
      expect(screen.getByText("+ Peak")).toHaveAttribute("aria-pressed", "true");
      fireEvent.keyDown(document.body, { key: "ArrowDown" }); // preview first candidate
      expect(
        screen.getByRole("button", { name: /Pn3m, in assignment/ }),
      ).toHaveAttribute("data-previewed", "true");
      // A [ / ] step changes activeSampleId WITHOUT remounting FocusPage. In
      // production Zustand re-renders on that change; the mocked store is
      // non-reactive, so drive the sample change + re-render the SAME tree
      // (FocusPage is reused, not remounted — its useState survives). The reset
      // effect keyed on activeSampleId must clear the arm, the zoom, AND the
      // candidate preview (a stale preview would otherwise eat the first Escape).
      act(() => {
        state.activeSampleId = 43;
      });
      view.rerender(focusTreeAt(42));
      expect(screen.getByText("+ Peak")).toHaveAttribute("aria-pressed", "false");
      expect(
        screen.getByRole("button", { name: /Pn3m, in assignment/ }),
      ).not.toHaveAttribute("data-previewed");
    });

    it("→ / ← step the active exposure (no wrap)", () => {
      state.exposures = [exp({ id: 7 }), exp({ id: 8, filename: "JC042-002.dat" })];
      state.activeExposureId = 7;
      renderAt(42);
      fireEvent.keyDown(document.body, { key: "ArrowRight" });
      expect(state.activeExposureId).toBe(8);
      fireEvent.keyDown(document.body, { key: "ArrowLeft" });
      expect(state.activeExposureId).toBe(7);
    });

    // FO-EXPSKIP: the exposure axis traverses INDEXABLE exposures only — the
    // same acceptable set useAutoPickExposure pins to. Stepping onto a rejected
    // (dropped) frame would be reverted by the auto-pick (it yanks any
    // non-acceptable active exposure to the representative), so a step onto a
    // dropped frame reads as the axis going dead or jumping to the rep. The
    // stepper must skip rejected frames so ← / → moves between the frames the
    // page can actually hold active.
    it("→ / ← skip rejected (dropped) frames, stepping among acceptable exposures only (FO-EXPSKIP)", () => {
      state.exposures = [
        exp({ id: 7 }),
        exp({ id: 8, filename: "JC042-002.dat", status: "rejected", selected: false }),
        exp({ id: 9, filename: "JC042-003.dat" }),
      ];
      state.activeExposureId = 7;
      renderAt(42);
      // → skips the rejected middle frame (8) and lands on the next acceptable (9)
      fireEvent.keyDown(document.body, { key: "ArrowRight" });
      expect(state.activeExposureId).toBe(9);
      // clamp at the last acceptable frame — no wrap, and never onto a rejected one
      fireEvent.keyDown(document.body, { key: "ArrowRight" });
      expect(state.activeExposureId).toBe(9);
      // ← steps back to the first acceptable frame, also skipping the rejected one
      fireEvent.keyDown(document.body, { key: "ArrowLeft" });
      expect(state.activeExposureId).toBe(7);
    });

    it("↓ / ↑ move the previewed candidate (the arrow cursor), clamped, with a visible marker", () => {
      renderAt(42);
      const pn3m = screen.getByRole("button", { name: /Pn3m, in assignment/ });
      const lam = screen.getByRole("button", { name: /^Lamellar$/ });
      expect(pn3m).not.toHaveAttribute("data-previewed");
      fireEvent.keyDown(document.body, { key: "ArrowDown" }); // none → first
      expect(pn3m).toHaveAttribute("data-previewed", "true");
      fireEvent.keyDown(document.body, { key: "ArrowDown" }); // first → second
      expect(lam).toHaveAttribute("data-previewed", "true");
      expect(pn3m).not.toHaveAttribute("data-previewed");
      fireEvent.keyDown(document.body, { key: "ArrowDown" }); // clamp at last
      expect(lam).toHaveAttribute("data-previewed", "true");
      fireEvent.keyDown(document.body, { key: "ArrowUp" }); // back to first
      expect(pn3m).toHaveAttribute("data-previewed", "true");
    });

    it("Escape is a ladder: clear the candidate preview first, THEN back to the sheet", () => {
      renderAt(42);
      fireEvent.keyDown(document.body, { key: "ArrowDown" });
      const pn3m = () => screen.getByRole("button", { name: /Pn3m, in assignment/ });
      expect(pn3m()).toHaveAttribute("data-previewed", "true");
      // first Escape clears the preview, does NOT navigate
      fireEvent.keyDown(document.body, { key: "Escape" });
      expect(pn3m()).not.toHaveAttribute("data-previewed");
      expect(screen.getByTestId("loc")).toHaveAttribute("data-path", "/sample/42");
      // second Escape (no preview) backs out to the sheet
      fireEvent.keyDown(document.body, { key: "Escape" });
      expect(screen.getByTestId("loc")).toHaveAttribute("data-path", "/samples");
    });
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

  it("hovering a candidate lights its claimed peaks on the trace and dims the rest", () => {
    // Lamellar (index id 2) claims peak_id 2 only; Pn3m (id 1) claims peak 1.
    // Hovering Lamellar must highlight peak 2 and dim peak 1 — the restored
    // candidate-hover preview (now keyed on the candidate's CLAIM, so it works
    // even with an empty durable assignment post auto-group removal).
    const { container } = renderAt(42);
    // No ratio labels at rest (the labels layer is hover-only).
    expect(container.querySelector('[data-role="plot-labels"] [data-role="peak-label"]')).toBeNull();
    const lam = screen.getByRole("button", { name: /Lamellar/ });
    fireEvent.mouseEnter(lam.parentElement!); // the onMouseEnter preview wrapper
    const g1 = container.querySelector('[data-role="plot-peaks"] g[data-peak-id="1"]');
    const g2 = container.querySelector('[data-role="plot-peaks"] g[data-peak-id="2"]');
    expect(g1?.getAttribute("data-dimmed")).toBe("true"); // not claimed → dim
    expect(g2?.getAttribute("data-dimmed")).toBeNull(); // claimed → stays lit
    // The claimed peak now carries its ratio label (Lamellar pos 1 → "1").
    const labels = [...container.querySelectorAll('[data-role="plot-labels"] [data-role="peak-label"]')];
    expect(labels.map((l) => l.textContent)).toContain("1");
  });

  it("does not render the legacy '+ Add speculative' affordance (custom-index is the hypothesis tool)", () => {
    renderAt(42);
    expect(screen.queryByTestId("add-speculative-button")).toBeNull();
    expect(screen.queryByText("+ Add speculative")).toBeNull();
  });

  describe("form-factor declaration", () => {
    it("shows the form-factor row when the call is empty; marking it sets state + announces", () => {
      state.assignment = { exposure_id: 7, state: "indexed", members: [] };
      const announce = vi.fn();
      setAnnounceImpl(announce);
      try {
        renderAt(42);
        const row = screen.getByTestId("form-factor-row");
        expect(row).toHaveAttribute("aria-pressed", "false");
        fireEvent.click(row);
        expect(setAssignStateMutate).toHaveBeenCalledWith("form_factor");
        expect(announce.mock.calls.at(-1)?.[0]).toBe("Marked as form factor");
      } finally {
        setAnnounceImpl(null);
      }
    });

    it("reflects an existing form-factor declaration and clears it to unindexed (null)", () => {
      state.assignment = { exposure_id: 7, state: "form_factor", members: [] };
      const announce = vi.fn();
      setAnnounceImpl(announce);
      try {
        renderAt(42);
        const row = screen.getByTestId("form-factor-row");
        expect(row).toHaveAttribute("aria-pressed", "true");
        fireEvent.click(row);
        expect(setAssignStateMutate).toHaveBeenCalledWith("null");
        expect(announce.mock.calls.at(-1)?.[0]).toBe("Form factor cleared");
      } finally {
        setAnnounceImpl(null);
      }
    });

    it("hides the form-factor row once a phase is in the call (mutually exclusive with indexing)", () => {
      state.assignment = { exposure_id: 7, state: "indexed", members: [1] };
      renderAt(42);
      expect(screen.queryByTestId("form-factor-row")).toBeNull();
    });
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

  it("arming '+ Peak' then Enter on an AUTO peak disables it (toggle exclude), not remove", () => {
    const announce = vi.fn();
    setAnnounceImpl(announce);
    try {
      renderAt(42); // default fixture peak id 1 is source:"auto"
      fireEvent.click(screen.getByText("+ Peak"));
      const mark = screen.getByRole("button", { name: /Auto peak at q = 0\.2000/ });
      fireEvent.keyDown(mark, { key: "Enter" });
      // Auto peaks belong to the indexer — a click disables (excludes) them,
      // it does NOT remove them (the old plain-click=remove failed on auto).
      expect(setPeakExclMutate).toHaveBeenCalledWith({ peakId: 1, excluded: true });
      expect(removePeakMutate).not.toHaveBeenCalled();
      expect(announce.mock.calls.at(-1)?.[0]).toBe("Auto peak disabled");
    } finally {
      setAnnounceImpl(null);
    }
  });

  it("Enter on an already-disabled AUTO peak restores it", () => {
    const announce = vi.fn();
    setAnnounceImpl(announce);
    state.peaks = [
      { id: 1, exposure_id: 7, q: 0.2, intensity: 40, prominence: 10, sharpness: 2, source: "auto", excluded: true },
    ];
    try {
      renderAt(42);
      fireEvent.click(screen.getByText("+ Peak"));
      const mark = screen.getByRole("button", { name: /Auto peak at q = 0\.2000 \(excluded\)/ });
      fireEvent.keyDown(mark, { key: "Enter" });
      expect(setPeakExclMutate).toHaveBeenCalledWith({ peakId: 1, excluded: false });
      expect(announce.mock.calls.at(-1)?.[0]).toBe("Auto peak restored");
    } finally {
      setAnnounceImpl(null);
    }
  });

  it("arming '+ Peak' then Enter on a MANUAL peak removes it and announces", () => {
    const announce = vi.fn();
    setAnnounceImpl(announce);
    state.peaks = [
      { id: 5, exposure_id: 7, q: 0.2, intensity: 40, prominence: null, sharpness: null, source: "manual", excluded: false },
    ];
    try {
      renderAt(42);
      fireEvent.click(screen.getByText("+ Peak"));
      const mark = screen.getByRole("button", { name: /Manual peak at q = 0\.2000/ });
      fireEvent.keyDown(mark, { key: "Enter" });
      expect(removePeakMutate).toHaveBeenCalledWith(5);
      expect(setPeakExclMutate).not.toHaveBeenCalled();
      expect(announce.mock.calls.at(-1)?.[0]).toBe("Peak removed");
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
    // The genuine empty state is NOT a loading state: a sample with zero usable
    // exposures must settle to the "no exposures" panel, not hold the skeleton.
    expect(screen.getByTestId("focus-skeleton-gate")).toHaveAttribute(
      "data-loading",
      "false",
    );
  });

  // FO-NAV-SKELETON: switching samples ([ ]) clears activeExposureId before the
  // new sample's exposures round-trip. The OLD skeleton gate only checked
  // traceQ/peaksQ behind an `activeExposureId !== undefined` guard, so during
  // that window it read false and flashed empty panels before the boneyard
  // appeared. The gate must hold the skeleton across the whole exposure-
  // resolution window.
  it("holds the skeleton while the sample's exposures are still loading (no empty flash)", () => {
    state.activeSampleId = 42;
    state.activeExposureId = undefined; // cleared by the sample switch
    state.corpus = [corpus()];
    state.exposures = undefined; // still fetching → which exposure is unknown
    state.trace = undefined;
    state.peaks = [];
    state.loading = false; // isolate the resolvingExposure branch from corpus
    renderAt(42);
    expect(screen.getByTestId("focus-skeleton-gate")).toHaveAttribute(
      "data-loading",
      "true",
    );
    // and it must NOT have fallen through to the genuine "no exposures" panel
    expect(screen.queryByText(/no exposures/i)).toBeNull();
  });

  // The caching win: navigating onto a sample whose exposures AND trace are
  // already cached must resolve to real content in the same render — even though
  // the STORED activeExposureId is momentarily undefined (the page resolves the
  // representative at render time via resolveActiveExposure). No skeleton.
  it("shows cached content instantly when the stored exposure id is undefined but exposures+trace are warm", () => {
    seedFull();
    state.activeExposureId = undefined; // store not yet re-seeded after a switch
    // exposures (id 7, selected) + trace + peaks are all cached from seedFull
    renderAt(42);
    expect(screen.getByTestId("focus-skeleton-gate")).toHaveAttribute(
      "data-loading",
      "false",
    );
    // and the workspace actually resolved to the cached exposure's content
    expect(screen.getByTestId("trace-plate")).toBeInTheDocument();
  });

  it("settles the skeleton once the active exposure resolves and its trace is in hand", () => {
    seedFull(); // activeExposureId 7, trace + peaks present, loading false
    renderAt(42);
    expect(screen.getByTestId("focus-skeleton-gate")).toHaveAttribute(
      "data-loading",
      "false",
    );
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
