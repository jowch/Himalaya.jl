// test/FocusPage.return.test.tsx
// T3.2: Focus Esc dismiss — returns to /series when from=series, else to
// /experiments/:expId/corpus
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import type { CorpusSample, Exposure, IndexEntry, Assignment, Trace } from "../src/api";

// ── module mocks (must precede import of FocusPage) ──────────────────────────

const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

const state = {
  activeSampleId: 42 as number | undefined,
  activeExposureId: 7 as number | undefined,
  corpus: [] as CorpusSample[],
  exposures: [] as Exposure[],
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
  trace: undefined as Trace | undefined,
};

vi.mock("../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.corpus, isLoading: false }),
  useSamples: () => ({ data: [], isLoading: false }),
  useExperiment: () => ({ data: { id: 7, name: "BL-19", q_units: "A-1" } }),
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
    prev: undefined,
    next: undefined,
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

import { FocusPage } from "../src/print/pages/FocusPage";

// ── helpers ───────────────────────────────────────────────────────────────────
function mkSample(id: number, experiment_id: number): CorpusSample {
  return {
    id, experiment_id, name: `S${id}`, display_name: null,
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

function renderFocus(route: string) {
  return render(
    <MemoryRouter initialEntries={[route]}>
      <Routes>
        <Route path="/sample/:sampleId" element={<FocusPage />} />
        <Route path="/series" element={<div data-testid="series-page">series</div>} />
        <Route path="/experiments/:id/corpus" element={<div data-testid="corpus-page">corpus</div>} />
      </Routes>
    </MemoryRouter>,
  );
}

beforeEach(() => {
  vi.clearAllMocks();
  state.activeSampleId = 42;
  state.activeExposureId = 7;
  state.corpus = [mkSample(42, 7)];
  state.exposures = [mkExposure(7, 42)];
  state.indices = [];
  state.assignment = undefined;
  state.trace = { q: [0.1, 0.2], I: [10, 20], sigma: [1, 1] };
  navigate.mockClear();
});

// ── T3.2 dismiss return target ────────────────────────────────────────────────
describe("Focus T3.2: Esc dismiss return target", () => {
  it("Escape returns to /series when opened with ?from=series", () => {
    renderFocus("/sample/42?from=series");
    fireEvent.keyDown(window, { key: "Escape" });
    expect(navigate).toHaveBeenLastCalledWith("/series");
  });

  it("Escape returns to /experiments/:expId/corpus when no from param", () => {
    renderFocus("/sample/42");
    fireEvent.keyDown(window, { key: "Escape" });
    expect(navigate).toHaveBeenLastCalledWith("/experiments/7/corpus");
  });
});
