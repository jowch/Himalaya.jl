// test/FocusPage.keyboard.test.tsx
//
// T2.5: Focus rev-2 axis pins.
//
// After T2.5:
//   ←/→ = prevDetail/nextDetail → step the candidate preview index.
//   ↑/↓ = prevSample/nextSample → navigate to sibling sample.
//   The old exposure axis (prevExposure/nextExposure on [ ]/←/→) is GONE —
//   exposure stepping is filmstrip-only via ThumbnailGallery.
//   l key = openLoupe → navigate to /sample/:id/loupe.
//
// These tests complement test/print-pages/FocusPage.test.tsx, which covers
// the full page render. Here we focus purely on the keyboard contract.
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import type { CorpusSample, Exposure, IndexEntry, Trace, Assignment } from "../src/api";

// ── mutator spies ─────────────────────────────────────────────────────────────
vi.mock("../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.corpus, isLoading: false }),
  useSamples: () => ({ data: [], isLoading: false }),
  useExperiment: () => ({ data: { id: 1, name: "BL-19", q_units: "A-1" } }),
  useExposures: () => ({ data: state.exposures, isLoading: false }),
  useExposure: () => ({ data: undefined }),
  useTrace: () => ({ data: state.trace, isLoading: false }),
  usePeaks: () => ({ data: [], isLoading: false }),
  useIndices: () => ({ data: state.indices, isLoading: false }),
  useAssignment: () => ({ data: state.assignment, isLoading: false }),
  useAddPeak: () => ({ mutate: vi.fn() }),
  useRemovePeak: () => ({ mutate: vi.fn() }),
  useSetPeakExcluded: () => ({ mutate: vi.fn() }),
  useAddAssignmentPhase: () => ({ mutate: vi.fn() }),
  useRemoveAssignmentPhase: () => ({ mutate: vi.fn() }),
  useSetAssignmentState: () => ({ mutate: vi.fn() }),
  useCommitCustomIndex: () => ({ mutate: vi.fn() }),
}));

vi.mock("../src/hooks/useAutoPickExposure", async (importOriginal) => ({
  ...(await importOriginal<typeof import("../src/hooks/useAutoPickExposure")>()),
  useAutoPickExposure: () => {},
}));

const state = {
  corpus: [] as CorpusSample[],
  exposures: [] as Exposure[],
  trace: undefined as Trace | undefined,
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
  activeSampleId: undefined as number | undefined,
  activeExposureId: undefined as number | undefined,
  sibPrev: undefined as { id: number } | undefined,
  sibNext: undefined as { id: number } | undefined,
};

vi.mock("../src/state", () => ({
  useAppState: (sel: (s: Record<string, unknown>) => unknown) =>
    sel({
      activeSampleId: state.activeSampleId,
      activeExposureId: state.activeExposureId,
      setActiveSample: (id: number | undefined) => { state.activeSampleId = id; },
      setActiveExposure: (id: number | undefined) => { state.activeExposureId = id; },
    }),
}));

vi.mock("../src/hooks/useExperimentSiblings", () => ({
  useExperimentSiblings: () => ({
    activeSample: state.activeSampleId !== undefined ? { id: state.activeSampleId } : undefined,
    siblings: [],
    index: 0,
    prev: state.sibPrev,
    next: state.sibNext,
  }),
}));

vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

vi.mock("../src/print/detector", async (orig) => {
  const real = (await orig()) as Record<string, unknown>;
  return { ...real, DetectorImage: () => <div data-testid="mock-detector-image" /> };
});
vi.mock("../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

vi.mock("../src/print/components/CombsPanel", () => ({
  CombsPanel: () => <div data-testid="combs-panel" />,
}));

// navigate spy wired before the import so the module-level mock is active.
const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

import { FocusPage } from "../src/print/pages/FocusPage";

// ── helpers ───────────────────────────────────────────────────────────────────
function mkCorpusSample(id: number): CorpusSample {
  return {
    id, experiment_id: 1, name: `S${id}`, display_name: null,
    notes: null, q_units: "A-1", tags: [],
  } as CorpusSample;
}
function mkExposure(id: number, sampleId: number): Exposure {
  return {
    id, sample_id: sampleId, filename: `file_${id}.dat`, kind: "file" as const,
    selected: false, status: "accepted", image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
  };
}
function mkIndex(id: number): IndexEntry {
  return {
    id, exposure_id: 7, phase: id === 1 ? "Pn3m" : "Lamellar",
    basis: 0.15, score: 0.9, r_squared: 0.99, lattice_d: 197,
    ngc: -1.5, status: "candidate", kind: "auto", inputs_hash: "h",
    peaks: [], predicted_q: [],
  };
}

function LocationProbe(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="loc" data-path={loc.pathname} />;
}

function renderFocus(sampleId = 42) {
  return render(
    <MemoryRouter initialEntries={[`/sample/${sampleId}`]}>
      <LocationProbe />
      <Routes>
        <Route path="/sample/:sampleId" element={<FocusPage />} />
        <Route path="/samples" element={<div data-testid="sheet">sheet</div>} />
        <Route path="/sample/:sampleId/loupe" element={<div data-testid="loupe">loupe</div>} />
      </Routes>
    </MemoryRouter>,
  );
}

beforeEach(() => {
  vi.clearAllMocks();
  state.activeSampleId = 42;
  state.activeExposureId = 7;
  state.corpus = [mkCorpusSample(42)];
  state.exposures = [mkExposure(7, 42)];
  state.trace = { q: [0.1, 0.2], I: [10, 20], sigma: [1, 1] };
  state.indices = [mkIndex(1), mkIndex(2)];
  state.assignment = { exposure_id: 7, state: "indexed", members: [1] };
  state.sibPrev = undefined;
  state.sibNext = undefined;
  navigate.mockClear();
});

// ── T2.5 axis contract ────────────────────────────────────────────────────────

describe("Focus T2.5 keyboard: rev-2 axis contract", () => {
  it("←/→ step the candidate preview (prevDetail/nextDetail)", () => {
    renderFocus();
    const pn3m = screen.getByRole("button", { name: /Pn3m, in assignment/ });
    const lam = screen.getByRole("button", { name: /^Lamellar$/ });

    // Start: no preview.
    expect(pn3m).not.toHaveAttribute("data-previewed");

    // ArrowRight → first candidate (Pn3m).
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(pn3m).toHaveAttribute("data-previewed", "true");

    // ArrowRight again → second candidate (Lamellar).
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(lam).toHaveAttribute("data-previewed", "true");
    expect(pn3m).not.toHaveAttribute("data-previewed");

    // ArrowLeft → back to first.
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    expect(pn3m).toHaveAttribute("data-previewed", "true");
  });

  it("↑/↓ step the sample (prevSample/nextSample)", () => {
    state.sibNext = { id: 43 };
    state.sibPrev = { id: 41 };
    renderFocus();

    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(navigate).toHaveBeenCalledWith("/sample/43");

    navigate.mockClear();
    fireEvent.keyDown(window, { key: "ArrowUp" });
    expect(navigate).toHaveBeenCalledWith("/sample/41");
  });

  it("exposure keys are gone — ArrowRight does NOT change activeExposureId", () => {
    state.exposures = [mkExposure(7, 42), mkExposure(8, 42)];
    state.activeExposureId = 7;
    renderFocus();

    // ArrowRight now steps the candidate, NOT the exposure.
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(state.activeExposureId).toBe(7); // unchanged
  });

  it("l key (openLoupe) navigates to /sample/:id/loupe", () => {
    renderFocus();
    fireEvent.keyDown(window, { key: "l" });
    expect(navigate).toHaveBeenCalledWith("/sample/42/loupe");
  });

  it("[ and ] are no longer bound — neither changes the active sample", () => {
    state.sibNext = { id: 43 };
    state.sibPrev = { id: 41 };
    renderFocus();

    fireEvent.keyDown(window, { key: "]" });
    fireEvent.keyDown(window, { key: "[" });
    // [ and ] are not in the registry; navigate must not be called
    expect(navigate).not.toHaveBeenCalledWith("/sample/43");
    expect(navigate).not.toHaveBeenCalledWith("/sample/41");
  });
});
