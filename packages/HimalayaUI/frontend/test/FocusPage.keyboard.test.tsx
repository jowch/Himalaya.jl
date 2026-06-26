// test/FocusPage.keyboard.test.tsx
//
// T2.5 + T4.1: Focus keyboard axis contract, reconciled to the new interaction
// architecture.
//
// After T4.1:
//   ←/→  = candidate cursor (on scope container, not window)
//   ↑/↓  = sample stepper (on scope container)
//   l    = openLoupe → navigate to /sample/:id/loupe (via keyboard layer)
//   p    = addPeak toggle (via keyboard layer)
//   Esc  = escapeLadder (via keyboard layer)
//   [ ]  = NOT bound (unchanged)
//   Exposure axis is filmstrip-only (unchanged)
//
// Arrows fire on [data-interaction-scope]; other keys fire on window.
import { describe, it, expect, vi, beforeEach } from "vitest";
import type React from "react";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import type { CorpusSample, Exposure, IndexEntry, Trace, Assignment } from "../src/api";
import { FocusPage } from "../src/print/pages/FocusPage";
import { InteractionDock } from "../src/print/interaction/InteractionDock";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../src/print/interaction/registry";

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
  // siblings drives useStepperOnly — must be the full ordered list
  siblings: [] as { id: number }[],
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
  useExperimentSiblings: () => {
    const i = state.siblings.findIndex((s) => s.id === state.activeSampleId);
    return {
      activeSample: state.activeSampleId !== undefined ? { id: state.activeSampleId } : undefined,
      siblings: state.siblings,
      index: i,
      prev: state.siblings[i - 1],
      next: state.siblings[i + 1],
    };
  },
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

const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

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
        <Route path="/sample/:sampleId" element={<TestShell><FocusPage /></TestShell>} />
        <Route path="/samples" element={<div data-testid="sheet">sheet</div>} />
        <Route path="/sample/:sampleId/loupe" element={<div data-testid="loupe">loupe</div>} />
      </Routes>
    </MemoryRouter>,
  );
}

function getScope(): HTMLElement {
  const scope = document.querySelector("[data-interaction-scope]") as HTMLElement | null;
  if (!scope) throw new Error("No [data-interaction-scope] found");
  return scope;
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
  state.siblings = [{ id: 42 }]; // default: lone sample, no navigation
  navigate.mockClear();
  useInteraction.getState().clearPage();
});

// ── T2.5 / T4.1 axis contract ─────────────────────────────────────────────────

describe("Focus T4.1 keyboard: candidate + sample axis via scope container", () => {
  it("←/→ step the candidate preview (on scope, not window)", async () => {
    // cursor starts at ids[0]=1 (Pn3m) but previewWasExplicit=false → nothing shown.
    // ArrowLeft: clamps at index 0 (Pn3m stays), sets previewWasExplicit=true → Pn3m shown.
    // ArrowRight: moves to index 1 (Lamellar), previewWasExplicit already true → Lamellar shown.
    // ArrowLeft: moves back to index 0 (Pn3m) → Pn3m shown.
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    const pn3m = screen.getByRole("button", { name: /Pn3m, in assignment/ });
    const lam = screen.getByRole("button", { name: /^Lamellar$/ });

    // Start: no preview.
    expect(pn3m).not.toHaveAttribute("data-previewed", "true");

    // ArrowLeft on scope → cursor clamps at Pn3m (index 0), becomes explicit.
    act(() => { fireEvent.keyDown(scope, { key: "ArrowLeft" }); });
    expect(pn3m).toHaveAttribute("data-previewed", "true");

    // ArrowRight → second candidate (Lamellar, index 1).
    act(() => { fireEvent.keyDown(scope, { key: "ArrowRight" }); });
    expect(lam).toHaveAttribute("data-previewed", "true");
    expect(pn3m).not.toHaveAttribute("data-previewed", "true");

    // ArrowLeft → back to first (Pn3m).
    act(() => { fireEvent.keyDown(scope, { key: "ArrowLeft" }); });
    expect(pn3m).toHaveAttribute("data-previewed", "true");
  });

  it("↑/↓ step the sample via scope container", async () => {
    state.siblings = [{ id: 41 }, { id: 42 }, { id: 43 }];
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();

    act(() => { fireEvent.keyDown(scope, { key: "ArrowDown" }); });
    expect(navigate).toHaveBeenCalledWith("/sample/43");

    navigate.mockClear();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowUp" }); });
    expect(navigate).toHaveBeenCalledWith("/sample/41");
  });

  it("exposure keys are gone — ArrowRight does NOT change activeExposureId", async () => {
    state.exposures = [mkExposure(7, 42), mkExposure(8, 42)];
    state.activeExposureId = 7;
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();

    // ArrowRight on scope now steps the candidate, NOT the exposure.
    act(() => { fireEvent.keyDown(scope, { key: "ArrowRight" }); });
    expect(state.activeExposureId).toBe(7); // unchanged
  });

  it("l key (openLoupe) navigates to /sample/:id/loupe", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    fireEvent.keyDown(window, { key: "l" });
    expect(navigate).toHaveBeenCalledWith("/sample/42/loupe");
  });

  it("[ and ] are no longer bound — neither changes the active sample", async () => {
    state.siblings = [{ id: 41 }, { id: 42 }, { id: 43 }];
    renderFocus();
    await screen.findByTestId("focus-workspace");

    fireEvent.keyDown(window, { key: "]" });
    fireEvent.keyDown(window, { key: "[" });
    // [ and ] are not in the registry; navigate must not be called
    expect(navigate).not.toHaveBeenCalledWith("/sample/43");
    expect(navigate).not.toHaveBeenCalledWith("/sample/41");
  });
});
