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
const deleteIndexMutate = vi.fn();
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
  // Trace/peaks loading INDEPENDENT of structure (corpus/exposures): the page
  // reveals on structure and only the trace plot waits on this.
  traceLoading: false,
  activeSampleId: undefined as number | undefined,
  activeExposureId: undefined as number | undefined,
  // sibling order for the [ ] sample-step shortcut (useExperimentSiblings mock)
  sibPrev: undefined as { id: number } | undefined,
  sibNext: undefined as { id: number } | undefined,
  // Dock sample readout (5c): seeded sibling list + index, default empty so
  // tests that don't care render no readout.
  sibSiblings: [] as { id: number }[],
  sibIndex: 0,
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
  useTrace: () => ({ data: state.trace, isLoading: state.traceLoading || state.loading }),
  usePeaks: () => ({ data: state.peaks, isLoading: state.traceLoading || state.loading }),
  useIndices: () => ({ data: state.indices, isLoading: state.loading }),
  useAssignment: () => ({ data: state.assignment, isLoading: state.loading }),
  useAddPeak: () => ({ mutate: addPeakMutate }),
  useRemovePeak: () => ({ mutate: removePeakMutate }),
  useSetPeakExcluded: () => ({ mutate: setPeakExclMutate }),
  useAddAssignmentPhase: () => ({ mutate: addAssignMutate }),
  useRemoveAssignmentPhase: () => ({ mutate: removeAssignMutate }),
  useDeleteIndex: () => ({ mutate: deleteIndexMutate }),
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
    siblings: state.sibSiblings,
    index: state.sibIndex,
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
import { InteractionDock } from "../../src/print/interaction/InteractionDock";
import { useKeyboardLayer } from "../../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../../src/print/interaction/registry";
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
  state.traceLoading = false;
  state.sibPrev = undefined;
  state.sibNext = undefined;
  state.sibSiblings = [];
  state.sibIndex = 0;
}

function LocationProbe() {
  const loc = useLocation();
  return <div data-testid="loc" data-path={loc.pathname} />;
}

function TestShell({ children }: { children: React.ReactNode }): JSX.Element {
  useKeyboardLayer();
  return <>{children}<InteractionDock /></>;
}

function getScope(): HTMLElement {
  const el = document.querySelector("[data-interaction-scope]") as HTMLElement | null;
  if (!el) throw new Error("No [data-interaction-scope] found");
  return el;
}

function focusTreeAt(sampleId: number | string) {
  return (
    <MemoryRouter initialEntries={[`/sample/${sampleId}`]}>
      <LocationProbe />
      <Routes>
        <Route path="/sample/:sampleId" element={<TestShell><FocusPage /></TestShell>} />
        <Route path="/experiments/:id/corpus" element={<div data-testid="sheet">corpus</div>} />
        <Route path="/experiments" element={<div data-testid="experiments-home">home</div>} />
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
  useInteraction.getState().clearPage();
});

describe("FocusPage", () => {
  it("renders the trace plate, assignment rail, detector panel and combs panel", () => {
    renderAt(42);
    expect(screen.getByTestId("trace-plate")).toBeInTheDocument();
    expect(screen.getByTestId("assignment-rail")).toBeInTheDocument();
    expect(screen.getByTestId("detector-panel")).toBeInTheDocument();
    expect(screen.getByTestId("combs-panel")).toBeInTheDocument();
  });

  it("dock follows the §7 grammar: labeled Sample readout + right-anchored Loupe with L chip (5c)", () => {
    // Two experiment-siblings so the readout reports a real total.
    state.sibSiblings = [{ id: 42 }, { id: 43 }];
    state.sibIndex = 0;
    renderAt(42);
    expect(screen.getByTestId("dock-sample-count").textContent).toBe("1 / 2");
    // Dock key-caps route through keyComboLabel now (review round 1): single
    // letters render UPPERCASE, matching the help legend (and combos get glyphs:
    // Mod+Enter → ⌘↵, Space → space).
    expect(within(screen.getByTestId("dock-action-openLoupe")).getByTestId("kbkey").textContent).toBe("L");
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

  // ── keyboard: the two-axis model (shared shortcut library, rev-2 axes T2.5) ──
  // ↑/↓ = sample axis (prevSample/nextSample), ←/→ = candidate detail axis
  // (prevDetail/nextDetail). The old exposure axis ([]/[ and old ←/→ for frames
  // on Focus) is retired — exposure stepping is filmstrip-only via ThumbnailGallery.
  describe("keyboard — two-axis model", () => {
    it("↓ / ↑ step the sample (agree with the stepper via useExperimentSiblings)", () => {
      // useStepperOnly uses siblings.map(id) for its ordered list; prev/next alone
      // don't feed it. Provide a full ordered list so onNext can resolve ids[i+1].
      state.sibSiblings = [{ id: 41 }, { id: 42 }, { id: 43 }];
      state.sibNext = { id: 43 };
      state.sibPrev = { id: 41 };
      renderAt(42);
      fireEvent.keyDown(getScope(), { key: "ArrowDown" });
      expect(screen.getByTestId("loc")).toHaveAttribute("data-path", "/sample/43");
    });

    it("↑ steps to the previous sample", () => {
      state.sibSiblings = [{ id: 41 }, { id: 42 }];
      state.sibPrev = { id: 41 };
      renderAt(42);
      fireEvent.keyDown(getScope(), { key: "ArrowUp" });
      expect(screen.getByTestId("loc")).toHaveAttribute("data-path", "/sample/41");
    });

    // FO-NAV-STATE: React Router does NOT remount FocusPage on a same-route
    // ↑/↓ sample step, so page-owned interaction state survives a sample switch.
    // The "+ Peak" arm is INTENTIONALLY sticky (like scale/combView): once armed,
    // stepping keeps you armed so you can edit peaks sample-to-sample. But the
    // zoom window and candidate preview must still reset — a stale x-domain
    // renders the next trace off-screen, and a stale preview eats the first Esc.
    it("keeps the '+ Peak' arm but resets the candidate preview when the active sample changes (FO-NAV-STATE)", () => {
      seedFull();
      state.corpus = [corpus(), corpus({ id: 43, name: "JC043" })];
      const view = renderAt(42);
      // Arm "+ Peak" AND preview a candidate.
      fireEvent.click(screen.getAllByText("+ Peak")[0]);
      expect(screen.getAllByText("+ Peak")[0]).toHaveAttribute("aria-pressed", "true");
      fireEvent.keyDown(getScope(), { key: "ArrowLeft" }); // clamp-at-first → pn3m explicit
      expect(
        screen.getByRole("button", { name: /Pn3m, in assignment/ }),
      ).toHaveAttribute("data-previewed", "true");
      // A ↑/↓ step changes activeSampleId WITHOUT remounting FocusPage. In
      // production Zustand re-renders on that change; the mocked store is
      // non-reactive, so drive the sample change + re-render the SAME tree
      // (FocusPage is reused, not remounted — its useState survives). The reset
      // effect keyed on activeSampleId clears the candidate preview (asserted
      // below; the zoom also resets but is awkward to assert via public DOM), and
      // deliberately leaves the arm ON.
      act(() => {
        state.activeSampleId = 43;
      });
      view.rerender(focusTreeAt(42));
      expect(screen.getAllByText("+ Peak")[0]).toHaveAttribute("aria-pressed", "true");
      expect(
        screen.getByRole("button", { name: /Pn3m, in assignment/ }),
      ).not.toHaveAttribute("data-previewed");
    });

    it("→ / ← step the previewed candidate (detail axis), clamped, with a visible marker", () => {
      renderAt(42);
      const pn3m = screen.getByRole("button", { name: /Pn3m, in assignment/ });
      const lam = screen.getByRole("button", { name: /^Lamellar$/ });
      // cursor starts at Pn3m (ids[0]); ArrowLeft clamps at index 0 → pn3m explicit.
      expect(pn3m).not.toHaveAttribute("data-previewed");
      fireEvent.keyDown(getScope(), { key: "ArrowLeft" }); // clamp at first → pn3m
      expect(pn3m).toHaveAttribute("data-previewed", "true");
      fireEvent.keyDown(getScope(), { key: "ArrowRight" }); // first → second (lamellar)
      expect(lam).toHaveAttribute("data-previewed", "true");
      expect(pn3m).not.toHaveAttribute("data-previewed");
      fireEvent.keyDown(getScope(), { key: "ArrowRight" }); // clamp at last
      expect(lam).toHaveAttribute("data-previewed", "true");
      fireEvent.keyDown(getScope(), { key: "ArrowLeft" }); // back to first
      expect(pn3m).toHaveAttribute("data-previewed", "true");
    });

    it("exposure keys (old [] / old ←→ exposure axis) are gone — arrows step candidate, not exposure", () => {
      // ArrowRight must NOT change activeExposureId; it steps the candidate preview.
      state.exposures = [exp({ id: 7 }), exp({ id: 8, filename: "JC042-002.dat" })];
      state.activeExposureId = 7;
      renderAt(42);
      fireEvent.keyDown(document.body, { key: "ArrowRight" });
      // exposure is unchanged — the candidate list moved, not the exposure
      expect(state.activeExposureId).toBe(7);
      // and no [ / ] binding exists either
      fireEvent.keyDown(document.body, { key: "]" });
      expect(state.activeExposureId).toBe(7);
    });

    it("Escape is a ladder: clear the candidate preview first, THEN back to the sheet", () => {
      renderAt(42);
      fireEvent.keyDown(getScope(), { key: "ArrowLeft" }); // clamp-at-first → pn3m explicit
      const pn3m = () => screen.getByRole("button", { name: /Pn3m, in assignment/ });
      expect(pn3m()).toHaveAttribute("data-previewed", "true");
      // first Escape clears the preview, does NOT navigate
      fireEvent.keyDown(document.body, { key: "Escape" });
      expect(pn3m()).not.toHaveAttribute("data-previewed");
      expect(screen.getByTestId("loc")).toHaveAttribute("data-path", "/sample/42");
      // second Escape (no preview) backs out to the experiment corpus
      // (sample 42 has experiment_id: 1, no ?from=series → /experiments/1/corpus)
      fireEvent.keyDown(document.body, { key: "Escape" });
      expect(screen.getByTestId("loc")).toHaveAttribute("data-path", "/experiments/1/corpus");
    });
  });

  it("toggling a candidate row fires useAddAssignmentPhase().mutate(indexId)", () => {
    renderAt(42);
    // Lamellar (id 2) is NOT in the active set → its CandidateRow toggles it on.
    const candidate = screen.getByRole("button", { name: /Lamellar/ });
    fireEvent.click(candidate);
    expect(addAssignMutate).toHaveBeenCalledWith(2);
  });

  it("offers Discard only on speculative candidates, never on auto ones", () => {
    // The route 403s on kind != 'speculative', so an auto index must not even
    // show the affordance. Pn3m (id 1) is auto; Lamellar (id 2) is made
    // speculative here to stand in for a committed custom index.
    state.indices = [ix({ id: 1, phase: "Pn3m" }),
                     ix({ id: 2, phase: "Lamellar", kind: "speculative", score: null })];
    renderAt(42);
    expect(screen.queryByRole("button", { name: /Discard the Pn3m index/ })).toBeNull();
    expect(screen.getByRole("button", { name: /Discard the Lamellar index/ })).toBeTruthy();
  });

  it("Discard fires useDeleteIndex().mutate(indexId) without toggling the assignment", () => {
    // The discard control is a SIBLING of CandidateRow, not a child — clicking
    // it must not also fire the row's add/remove toggle.
    state.indices = [ix({ id: 1, phase: "Pn3m" }),
                     ix({ id: 2, phase: "Lamellar", kind: "speculative", score: null })];
    renderAt(42);
    fireEvent.click(screen.getByRole("button", { name: /Discard the Lamellar index/ }));
    expect(deleteIndexMutate).toHaveBeenCalledWith(2);
    expect(addAssignMutate).not.toHaveBeenCalled();
    expect(removeAssignMutate).not.toHaveBeenCalled();
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
    const addPeakBtn = screen.getAllByText("+ Peak")[0];
    fireEvent.click(addPeakBtn);
    expect(addPeakBtn).toHaveAttribute("aria-pressed", "true");
    // The plot's add path: click empty plot space → interaction.onAddPeak(q).
    const svg = container.querySelector('svg[data-testid="trace-plate-plot"]')!;
    fireEvent.click(svg, { clientX: 300, clientY: 150 });
    expect(addPeakMutate).toHaveBeenCalledTimes(1);
    expect(typeof addPeakMutate.mock.calls[0]![0]).toBe("number");
  });

  // Guard for the sticky-arm trade-off: on a no-exposure sibling activeExposureId
  // is undefined, so useAddPeak falls back to the sentinel exposure 0. A carried
  // (or manual) arm must NOT POST a doomed /exposures/0/peaks nor announce a false
  // "Peak added" — onAddPeak no-ops when there is no exposure.
  it("an armed click on a no-exposure sample is an inert no-op (no sentinel POST, no false announce)", () => {
    const announce = vi.fn();
    setAnnounceImpl(announce);
    state.activeSampleId = 42;
    state.activeExposureId = undefined; // no representative exposure
    state.corpus = [corpus()];
    state.exposures = []; // loaded + none usable → noUsableExposure, not a skeleton
    state.trace = undefined;
    state.peaks = [];
    state.loading = false;
    const { container } = renderAt(42);
    // The interactive plate still renders (work column is ungated); arm it.
    const addPeakBtn = screen.getAllByText("+ Peak")[0];
    fireEvent.click(addPeakBtn);
    expect(addPeakBtn).toHaveAttribute("aria-pressed", "true");
    const svg = container.querySelector('svg[data-testid="trace-plate-plot"]')!;
    fireEvent.click(svg, { clientX: 300, clientY: 150 });
    expect(addPeakMutate).not.toHaveBeenCalled();
    expect(announce).not.toHaveBeenCalledWith("Peak added");
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
    // A claimless index (peaks: [] — e.g. a custom index committed before
    // insert_custom_index! claimed its landed peaks). With every predicted_q
    // within tol of an observed peak it emits NO coloured or ghost ring, so the
    // caption must not name it — only the ring-emitting sibling.
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

  it("custom-index modal opens seeded at the phase's def, with the widened floor as the slider min", () => {
    // Seeding from the widened floor (10 Å) would open an off-scale 1 nm Pn3m comb;
    // the initial value must come from def (197), while the widened floor reaches the
    // slider bound (paramMin) — the FO-QWINDOW-BOUNDS wiring.
    renderAt(42);
    fireEvent.click(screen.getByTestId("custom-index-trigger"));
    const num = within(screen.getByTestId("custom-index-modal")).getByRole("spinbutton");
    expect(num).toHaveValue(197);
    expect(num).toHaveAttribute("min", "10");
  });

  it("clicking an observed peak in the preview snaps the lattice onto it (FO-CIX-SNAP)", async () => {
    // Regression: onSelectObserved was dropped at this call site, silently
    // removing click-to-snap (component tests + Storybook stayed green). Guard
    // the wiring at the page level — the hit-targets render AND clicking one
    // sets the lattice to the consumer's snap formula.
    const { latticeForFirstOrderOnPeak } = await import("../../src/lib/customIndex");
    const { CUSTOM_SYMS } = await import("../../src/print/pages/focusAdapters");
    renderAt(42);
    fireEvent.click(screen.getByTestId("custom-index-trigger"));
    const modal = screen.getByTestId("custom-index-modal");
    const hit = modal.querySelector("[data-observed-hit]") as SVGElement | null;
    expect(hit).not.toBeNull(); // click-to-snap targets present (absent when onSelectObserved is unwired)
    const q = Number(hit!.getAttribute("data-observed-hit"));
    fireEvent.click(hit!);
    const num = within(modal).getByRole("spinbutton") as HTMLInputElement;
    expect(Number(num.value)).toBeCloseTo(latticeForFirstOrderOnPeak(CUSTOM_SYMS[0]!.name, q), 5);
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
      fireEvent.click(screen.getAllByText("+ Peak")[0]);
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
      fireEvent.click(screen.getAllByText("+ Peak")[0]);
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
      fireEvent.click(screen.getAllByText("+ Peak")[0]);
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
      fireEvent.click(screen.getAllByText("+ Peak")[0]);
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

  it("commits the NORMALIZED ratios the modal drew, not raw q's (Hexagonal)", () => {
    // FocusPage is the sole producer of the `ratios` wire value. If this
    // regressed to `refls.map(r => r.q)` the values are still positive and
    // non-empty, so route validation passes, compute_snap then matches zero
    // positions, and the commit returns 200 with an index claiming NO peaks —
    // exactly the bug this PR exists to fix, silently and with no 4xx.
    //
    // Hexagonal on purpose: its drawn set [1,3,4,7,9,12] is 6 of the backend
    // series' 13 entries, so this gives the "send ratios, not a count" contract
    // behavioural teeth at the producing layer. (Pre-#304 the two series also
    // disagreed in VALUE at position 6, where the core carried a √11 that is
    // not a permitted 2D hexagonal reflection; that entry is gone.)
    renderAt(42);
    fireEvent.click(screen.getByTestId("custom-index-trigger"));
    const modal = screen.getByTestId("custom-index-modal");
    // SegmentedControl defaults to role="group", so segments are buttons.
    fireEvent.click(within(modal).getByRole("button", { name: "Hexagonal" }));
    fireEvent.click(within(modal).getByRole("button", { name: /^Add/ }));

    expect(commitCustomMutate).toHaveBeenCalledTimes(1);
    const [phase, , ratios] = commitCustomMutate.mock.calls[0]!;
    expect(phase).toBe("Hexagonal");
    // q ∝ √M for hex, normalized against the first reflection (M = 1).
    const expected = [1, 3, 4, 7, 9, 12].map((m) => Math.sqrt(m));
    expect(ratios).toHaveLength(expected.length);
    ratios.forEach((r: number, i: number) =>
      expect(r).toBeCloseTo(expected[i]!, 12));
    // First entry is exactly 1 by construction — a raw-q regression would put
    // the basis here instead, which is O(0.1).
    expect(ratios[0]).toBe(1);
  });

  it("Escape sequence: an open custom-index modal closes first; only the NEXT Escape disarms '+ Peak' (F7)", () => {
    renderAt(42);
    const addPeakBtn = screen.getAllByText("+ Peak")[0];
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

  it("not-found renders an EmptyState whose action leads back to experiments (FO-ERR, T3.2)", () => {
    state.corpus = [];
    state.activeSampleId = undefined;
    renderAt(999);
    const block = screen.getByTestId("focus-not-found");
    expect(within(block).getByTestId("empty-state")).toBeInTheDocument();
    expect(
      within(block).getByRole("heading", { name: "Sample not found" }),
    ).toBeInTheDocument();
    fireEvent.click(
      within(block).getByRole("button", { name: "Back to the experiments" }),
    );
    expect(screen.getByTestId("experiments-home")).toBeInTheDocument();
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

  // Progressive load: the page reveals on STRUCTURE (corpus + exposure). The
  // trace/peaks fetch no longer gates the whole page — only the trace plot waits,
  // behind its own plate skeleton, so the detector + rail SHELL paint a
  // round-trip earlier. (Detector depends only on the exposure.)
  //
  // Models the REAL production load window: trace/peaks in flight means their
  // data is absent, so peaks/indices/assignment are empty here (not seedFull's
  // populated arrays — that loaded+data-present combination can't happen).
  it("reveals the page while trace/peaks still load; trace plot AND rail skeleton, no misleading copy", () => {
    seedFull();
    state.traceLoading = true;
    state.trace = undefined;
    state.peaks = [];
    state.indices = [];
    state.assignment = undefined;
    renderAt(42);
    // Page is NOT held behind the whole-page skeleton…
    expect(screen.getByTestId("focus-skeleton-gate")).toHaveAttribute(
      "data-loading",
      "false",
    );
    // …the plate shell renders, but its plot is the skeleton, not the live plot.
    expect(screen.getByTestId("trace-plate")).toBeInTheDocument();
    expect(screen.getByTestId("trace-plate-skeleton")).toBeInTheDocument();
    expect(screen.queryByTestId("trace-plate-plot")).toBeNull();
    // …and the rail shows a neutral loading body, NOT the peak-derived empty copy
    // that would tell the user to mark peaks on a trace that's still skeletoning.
    expect(screen.getAllByTestId("rail-body-skeleton").length).toBeGreaterThan(0);
    expect(screen.queryByText(/No peaks marked/i)).toBeNull();
    expect(screen.queryByText(/Candidates appear once peaks are marked/i)).toBeNull();
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
