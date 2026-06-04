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
const reanalyzeMutate = vi.fn();
const deleteIndexMutate = vi.fn();

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
  useReanalyzeExposure: () => ({ mutate: reanalyzeMutate }),
  useDeleteIndex: () => ({ mutate: deleteIndexMutate }),
}));

// Route shim: seed the Zustand active sample / exposure from the mock state.
vi.mock("../../src/hooks/useSyncActiveSampleFromRoute", () => ({
  useSyncActiveSampleFromRoute: () => {},
}));
vi.mock("../../src/hooks/useAutoPickExposure", () => ({
  useAutoPickExposure: () => {},
}));

// Zustand store: return the mock active ids.
vi.mock("../../src/state", () => ({
  useAppState: (sel: (s: Record<string, unknown>) => unknown) =>
    sel({
      activeSampleId: state.activeSampleId,
      activeExposureId: state.activeExposureId,
    }),
}));

// Pending-peak-ops gate (lives under lib/queue/hooks): no pending ops in tests.
vi.mock("../../src/lib/queue/hooks", () => ({
  useExposureHasPendingPeakOps: () => false,
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

function renderAt(sampleId: number) {
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

  it("'+ custom index…' opens the CustomIndexModal", () => {
    renderAt(42);
    fireEvent.click(screen.getByTestId("custom-index-trigger"));
    expect(screen.getByTestId("custom-index-modal")).toBeInTheDocument();
  });

  it("not-found: no corpus sample for the id", () => {
    state.corpus = [];
    state.activeSampleId = undefined;
    renderAt(999);
    expect(screen.getByTestId("focus-not-found")).toBeInTheDocument();
  });

  it("no-exposure: the sample has no exposures", () => {
    state.exposures = [];
    state.activeExposureId = undefined;
    state.trace = undefined;
    renderAt(42);
    expect(screen.getByText(/no exposures/i)).toBeInTheDocument();
  });

  it("shows the stale-index banner + reanalyzes when an index hash is stale", () => {
    state.indices = [ix({ id: 1, inputs_hash: "OLD" })];
    state.assignment = { exposure_id: 7, state: "indexed", members: [1] };
    renderAt(42);
    const banner = screen.getByRole("alert");
    expect(within(banner).getByText(/stale/)).toBeInTheDocument();
    fireEvent.click(within(banner).getByText("Re-analyze"));
    expect(reanalyzeMutate).toHaveBeenCalled();
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

  it("export-button is disabled (data gate) when there are no peaks", () => {
    state.peaks = [];
    renderAt(42);
    const plate = screen.getByTestId("trace-plate");
    const exportBtn = within(plate).getByTestId("export-button");
    // Every button inside the export widget should be disabled when the data gate is off.
    const copyBtn = within(exportBtn).getByTestId("export-copy");
    expect(copyBtn).toBeDisabled();
  });
});
