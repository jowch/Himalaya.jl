// test/FocusPage.actions.test.tsx
//
// Phase 4.1 — Focus interaction-architecture migration tests.
// Mounts the full TestShell (keyboard layer + InteractionDock) so dock buttons
// AND keyboard-triggered actions flow through the REAL registry path:
//   usePageActions → registry → InteractionDock / useKeyboardLayer → action.run
// Arrow keys (←/→/↑/↓) fire on the scope container (not window).
import { describe, it, expect, vi, beforeEach } from "vitest";
import type React from "react";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { CorpusSample, Exposure, IndexEntry, Assignment, Trace } from "../src/api";
import { FocusPage } from "../src/print/pages/FocusPage";
import { InteractionDock } from "../src/print/interaction/InteractionDock";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../src/print/interaction/registry";
import { useListCursor } from "../src/print/interaction/useListCursor";
import { DockStepper } from "../src/print/ui";
import { runCursorContract } from "./interaction/cursorContract";

// ── mock state ───────────────────────────────────────────────────────────────
const focusState = {
  corpus: [] as CorpusSample[],
  exposures: [] as Exposure[],
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
  trace: undefined as Trace | undefined,
  activeSampleId: undefined as number | undefined,
  activeExposureId: undefined as number | undefined,
  siblings: [] as { id: number }[],
};

const addAssignmentMutate = vi.fn();
const removeAssignmentMutate = vi.fn();

vi.mock("../src/queries", () => ({
  useCorpusSamples: () => ({ data: focusState.corpus, isLoading: false }),
  useExperiment: () => ({ data: { id: 1, name: "BL-19", q_units: "A-1" }, isLoading: false }),
  useExposures: () => ({ data: focusState.exposures, isLoading: false }),
  useTrace: () => ({ data: focusState.trace, isLoading: false }),
  usePeaks: () => ({ data: [], isLoading: false }),
  useIndices: () => ({ data: focusState.indices, isLoading: false }),
  useAssignment: () => ({ data: focusState.assignment, isLoading: false }),
  useAddPeak: () => ({ mutate: vi.fn() }),
  useRemovePeak: () => ({ mutate: vi.fn() }),
  useSetPeakExcluded: () => ({ mutate: vi.fn() }),
  useAddAssignmentPhase: () => ({ mutate: addAssignmentMutate }),
  useRemoveAssignmentPhase: () => ({ mutate: removeAssignmentMutate }),
  useDeleteIndex: () => ({ mutate: vi.fn() }),
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
      activeSampleId: focusState.activeSampleId,
      activeExposureId: focusState.activeExposureId,
      setActiveSample: (id: number | undefined) => { focusState.activeSampleId = id; },
      setActiveExposure: (id: number | undefined) => { focusState.activeExposureId = id; },
    }),
}));

vi.mock("../src/hooks/useExperimentSiblings", () => ({
  useExperimentSiblings: () => ({
    activeSample: focusState.activeSampleId !== undefined ? { id: focusState.activeSampleId } : undefined,
    siblings: focusState.siblings,
    index: focusState.siblings.findIndex((s) => s.id === focusState.activeSampleId),
    prev: focusState.siblings[focusState.siblings.findIndex((s) => s.id === focusState.activeSampleId) - 1],
    next: focusState.siblings[focusState.siblings.findIndex((s) => s.id === focusState.activeSampleId) + 1],
  }),
}));

vi.mock("boneyard-js/react", () => ({
  Skeleton: ({ children }: { children: React.ReactNode }) => <>{children}</>,
}));

vi.mock("../src/print/detector", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));
vi.mock("../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));
vi.mock("../src/print/components/CombsPanel", () => ({
  CombsPanel: () => <div data-testid="combs-panel" />,
}));

const navigateSpy = vi.fn();
vi.mock("react-router-dom", async (importOriginal) => {
  const actual = await importOriginal<typeof import("react-router-dom")>();
  return { ...actual, useNavigate: () => navigateSpy };
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
function makeSample(id: number, experimentId = 1): CorpusSample {
  return { id, experiment_id: experimentId, name: `S${id}`, display_name: null, notes: null, q_units: "A-1", tags: [] } as CorpusSample;
}
function makeExposure(id: number, sampleId: number): Exposure {
  return {
    id, sample_id: sampleId, filename: `file_${id}.dat`, kind: "file",
    selected: false, status: "accepted", image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
  };
}
function makeIndex(id: number, phase: string): IndexEntry {
  return {
    id, exposure_id: 7, phase, basis: 0.15, score: 0.9, r_squared: 0.99,
    lattice_d: 197, ngc: -1.5, status: "candidate", kind: "auto",
    inputs_hash: "h", peaks: [], predicted_q: [],
  };
}

function renderFocus(sampleId = 42, route = `/sample/${sampleId}`) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[route]}>
        <Routes>
          <Route
            path="/sample/:sampleId"
            element={<TestShell><FocusPage /></TestShell>}
          />
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

beforeEach(() => {
  navigateSpy.mockClear();
  addAssignmentMutate.mockClear();
  removeAssignmentMutate.mockClear();
  focusState.corpus = [makeSample(42, 1)];
  focusState.exposures = [makeExposure(7, 42)];
  focusState.trace = { q: [0.1, 0.2], I: [10, 20], sigma: [1, 1] };
  focusState.indices = [makeIndex(1, "Pn3m"), makeIndex(2, "Im3m")];
  focusState.assignment = { exposure_id: 7, state: "indexed", members: [] };
  focusState.activeSampleId = 42;
  focusState.activeExposureId = 7;
  focusState.siblings = [{ id: 41 }, { id: 42 }, { id: 43 }];
  useInteraction.getState().clearPage();
});

// ── tests ─────────────────────────────────────────────────────────────────────
describe("FocusPage action declaration", () => {
  it("dock primary label reads 'Apply' (not 'Focus')", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    // Primary button should be labelled Apply
    const primary = screen.getByTestId("dock-primary");
    expect(primary).toHaveTextContent("Apply");
  });

  it("Enter does NOT toggle assignment when no explicit preview (previewWasExplicit=false)", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    fireEvent.keyDown(window, { key: "Enter" });
    expect(addAssignmentMutate).not.toHaveBeenCalled();
    expect(removeAssignmentMutate).not.toHaveBeenCalled();
  });

  it("Enter toggles assignment after ArrowLeft makes first candidate explicit", async () => {
    // cursor starts at ids[0]=1 (Pn3m); ArrowLeft clamps at 0 → stays at 1,
    // sets previewWasExplicit=true. Enter then applies id=1.
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowLeft" }); });
    fireEvent.keyDown(window, { key: "Enter" });
    expect(addAssignmentMutate).toHaveBeenCalledWith(1);
  });

  it("Enter removes assignment when explicit candidate is already in the call", async () => {
    // id=2 (Im3m) is already in the call; ArrowRight moves cursor from 1→2.
    focusState.assignment = { exposure_id: 7, state: "indexed", members: [2] };
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowRight" }); });
    fireEvent.keyDown(window, { key: "Enter" });
    expect(removeAssignmentMutate).toHaveBeenCalledWith(2);
    expect(addAssignmentMutate).not.toHaveBeenCalled();
  });

  it("ArrowLeft on scope previews first candidate (cursor clamps at index 0)", async () => {
    // cursor starts at ids[0]=1 (Pn3m). ArrowLeft → clamps → stays at 1, sets preview.
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    const pn3mButton = screen.getByRole("button", { name: /Pn3m/ });
    expect(pn3mButton).not.toHaveAttribute("data-previewed", "true");
    act(() => { fireEvent.keyDown(scope, { key: "ArrowLeft" }); });
    expect(pn3mButton).toHaveAttribute("data-previewed", "true");
  });

  it("ArrowRight on scope previews second candidate (cursor advances from 0→1)", async () => {
    // cursor starts at ids[0]=1 (Pn3m). ArrowRight → ids[1]=2 (Im3m).
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    const im3mButton = screen.getByRole("button", { name: /Im3m/ });
    expect(im3mButton).not.toHaveAttribute("data-previewed", "true");
    act(() => { fireEvent.keyDown(scope, { key: "ArrowRight" }); });
    expect(im3mButton).toHaveAttribute("data-previewed", "true");
  });

  it("ArrowRight then ArrowLeft returns to first candidate", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    const pn3mButton = screen.getByRole("button", { name: /Pn3m/ });
    const im3mButton = screen.getByRole("button", { name: /Im3m/ });
    act(() => { fireEvent.keyDown(scope, { key: "ArrowRight" }); });
    expect(im3mButton).toHaveAttribute("data-previewed", "true");
    act(() => { fireEvent.keyDown(scope, { key: "ArrowLeft" }); });
    expect(pn3mButton).toHaveAttribute("data-previewed", "true");
  });

  it("p toggles addPeak mode (TracePlate button becomes armed)", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    // Before p: no button is armed
    expect(document.querySelector('[data-armed="true"]')).toBeNull();
    fireEvent.keyDown(window, { key: "p" });
    // After p: the + Peak button in TracePlate should be armed
    expect(document.querySelector('[data-armed="true"]')).not.toBeNull();
  });

  it("Escape exits addPeak first (rung 1), stays on page", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    // Arm addPeak
    fireEvent.keyDown(window, { key: "p" });
    expect(document.querySelector('[data-armed="true"]')).not.toBeNull();
    // Escape should disarm (not navigate)
    fireEvent.keyDown(window, { key: "Escape" });
    expect(navigateSpy).not.toHaveBeenCalled();
    expect(document.querySelector('[data-armed="true"]')).toBeNull();
  });

  it("Escape (rung 2) clears sticky preview after candidate was focused", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    // Make first candidate explicit via ArrowLeft (clamps at index 0 → Pn3m)
    act(() => { fireEvent.keyDown(scope, { key: "ArrowLeft" }); });
    const pn3mButton = screen.getByRole("button", { name: /Pn3m/ });
    expect(pn3mButton).toHaveAttribute("data-previewed", "true");
    // Escape clears preview (rung 2)
    fireEvent.keyDown(window, { key: "Escape" });
    expect(pn3mButton).not.toHaveAttribute("data-previewed", "true");
    expect(navigateSpy).not.toHaveBeenCalled();
  });

  it("Escape (rung 3) navigates out when no addPeak and no explicit preview", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    fireEvent.keyDown(window, { key: "Escape" });
    expect(navigateSpy).toHaveBeenCalledWith("/experiments/1/corpus");
  });

  it("Escape ladder: p → Esc(disarm) → arrow(explicit) → Esc(clear preview) → Esc(navigate)", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    // Step 1: arm addPeak
    fireEvent.keyDown(window, { key: "p" });
    // Step 2: Escape disarms addPeak (rung 1)
    fireEvent.keyDown(window, { key: "Escape" });
    expect(navigateSpy).not.toHaveBeenCalled();
    // Step 3: Move candidate (makes preview explicit, rung 2 target)
    act(() => { fireEvent.keyDown(scope, { key: "ArrowLeft" }); });
    // Step 4: Escape clears preview (rung 2)
    fireEvent.keyDown(window, { key: "Escape" });
    expect(navigateSpy).not.toHaveBeenCalled();
    // Step 5: Escape navigates out (rung 3)
    fireEvent.keyDown(window, { key: "Escape" });
    expect(navigateSpy).toHaveBeenCalledWith("/experiments/1/corpus");
  });

  it("ArrowDown on scope navigates to next sibling", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowDown" }); });
    expect(navigateSpy).toHaveBeenCalledWith("/sample/43");
  });

  it("ArrowUp on scope navigates to previous sibling", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const scope = getScope();
    act(() => { fireEvent.keyDown(scope, { key: "ArrowUp" }); });
    expect(navigateSpy).toHaveBeenCalledWith("/sample/41");
  });

  it("dock-next-sample navigates to next sibling", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    fireEvent.click(screen.getByTestId("dock-next-sample"));
    expect(navigateSpy).toHaveBeenCalledWith("/sample/43");
  });

  it("dock-prev-sample navigates to previous sibling", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    fireEvent.click(screen.getByTestId("dock-prev-sample"));
    expect(navigateSpy).toHaveBeenCalledWith("/sample/41");
  });

  it("renders both the candidate stepper and sample stepper in the dock", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    expect(screen.getByTestId("dock-prev-candidate")).toBeInTheDocument();
    expect(screen.getByTestId("dock-next-candidate")).toBeInTheDocument();
    expect(screen.getByTestId("dock-prev-sample")).toBeInTheDocument();
    expect(screen.getByTestId("dock-next-sample")).toBeInTheDocument();
  });

  it("sample stepper comes before candidate stepper in the dock DOM order", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    const sampleCount = screen.getByTestId("dock-sample-count");
    const candidateCount = screen.getByTestId("dock-candidate-count");
    // Sample (extraStepper) must appear before Candidate (cursor stepper)
    expect(sampleCount.compareDocumentPosition(candidateCount) & Node.DOCUMENT_POSITION_FOLLOWING).toBeTruthy();
  });

  it("l navigates to loupe", async () => {
    renderFocus();
    await screen.findByTestId("focus-workspace");
    fireEvent.keyDown(window, { key: "l" });
    expect(navigateSpy).toHaveBeenCalledWith("/sample/42/loupe");
  });

  it("from=series: Escape returns to /series", async () => {
    renderFocus(42, "/sample/42?from=series");
    await screen.findByTestId("focus-workspace");
    fireEvent.keyDown(window, { key: "Escape" });
    expect(navigateSpy).toHaveBeenCalledWith("/series");
  });
});

// ── cursor contract (candidate cursor) ───────────────────────────────────────
runCursorContract("Focus candidate cursor", () => {
  const onActivate = vi.fn();
  const IDS = [10, 20, 30];
  return {
    ui: (capture) => {
      function Probe(): JSX.Element {
        const cursor = useListCursor({
          ids: IDS,
          onActivate,
          stepperLabel: "Candidate",
          stepperTestIdBase: "candidate",
          axis: "horizontal",
        });
        capture({ cursor, ids: IDS, onActivate });
        return (
          <div data-interaction-scope="">
            {IDS.map((id) => (
              <div key={id} {...cursor.rowProps(id)}>
                candidate {id}
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
