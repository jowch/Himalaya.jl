// test/ExperimentCorpusPage.actions.test.tsx
//
// Corpus page action-declaration tests. These mount the full TestShell
// (keyboard layer + InteractionDock) so dock-primary AND keyboard-triggered
// actions flow through the REAL window listener (useKeyboardLayer) → registry →
// declared action.run — proving mouse/keyboard parity through one path.
// Navigate is mocked as a spy so we can assert on programmatic navigation
// without actual routing.
import { describe, it, expect, vi, beforeEach } from "vitest";
import type React from "react";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { CorpusSample, Exposure } from "../src/api";
import { ExperimentCorpusPage } from "../src/print/pages/ExperimentCorpusPage";
import { InteractionDock } from "../src/print/interaction/InteractionDock";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../src/print/interaction/registry";
import { useListCursor } from "../src/print/interaction/useListCursor";
import { DockStepper } from "../src/print/ui";
import { runCursorContract } from "./interaction/cursorContract";

// ── mock data plane ──────────────────────────────────────────────────────────
const actionsState = {
  samples: [] as CorpusSample[],
  byId: new Map<number, Exposure[]>(),
};
// Stable batch spy so keyboard-triggered cull verbs are assertable (the inline
// `vi.fn()` in the factory would mint a fresh fn per call and lose the calls).
const batchMutate = vi.fn();

vi.mock("../src/queries", () => ({
  useExperiment: () => ({ data: undefined, isLoading: false }),
  useLoads: () => ({ data: [], isLoading: false }),
  useTriggerScan: () => ({ mutate: vi.fn() }),
  useCorpusSamples: () => ({ data: actionsState.samples, isLoading: false, isError: false }),
  useCorpusExposures: (filtered: CorpusSample[]) => {
    const byId = new Map<number, Exposure[]>();
    for (const s of filtered) {
      const exps = actionsState.byId.get(s.id);
      if (exps !== undefined) byId.set(s.id, exps);
    }
    return { byId, isLoading: false };
  },
  useSetExposureStatusBatch: () => ({ mutate: batchMutate }),
}));

vi.mock("../src/state", () => ({
  useAppState: (sel: (s: Record<string, unknown>) => unknown) =>
    sel({ ingestInFlight: null, username: "tester" }),
}));

vi.mock("../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

// Mock navigate so programmatic navigation is assertable without real routing.
const navigateSpy = vi.fn();
vi.mock("react-router-dom", async (importOriginal) => {
  const actual = await importOriginal<typeof import("react-router-dom")>();
  return { ...actual, useNavigate: () => navigateSpy };
});

// ── shell (keyboard layer + dock) ────────────────────────────────────────────
function TestShell({ children }: { children: React.ReactNode }): JSX.Element {
  useKeyboardLayer();
  return (
    <>
      {children}
      <InteractionDock />
    </>
  );
}

// ── helpers ──────────────────────────────────────────────────────────────────
function makeSample(id: number, name: string): CorpusSample {
  return { id, experiment_id: 1, name, notes: null, q_units: "A-1", tags: [], phase: null } as CorpusSample;
}
function makeExposure(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 100, sample_id: 10, filename: "f.dat", kind: "file",
    selected: false, status: "accepted", image_path: "/x.tif", image_version: "v1",
    tags: [], sources: [], trace_hash: "h", analysis_inputs_hash: "a", ...over,
  };
}

function renderCorpus(
  samples: Array<{ id: number; name: string }>,
  opts: { expId?: number; byId?: Map<number, Exposure[]> } = {},
) {
  actionsState.samples = samples.map((s) => makeSample(s.id, s.name));
  actionsState.byId = opts.byId ?? new Map();
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[`/experiments/${opts.expId ?? 1}/corpus`]}>
        <Routes>
          <Route
            path="/experiments/:id/corpus"
            element={<TestShell><ExperimentCorpusPage /></TestShell>}
          />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

beforeEach(() => {
  navigateSpy.mockClear();
  batchMutate.mockClear();
  actionsState.samples = [];
  actionsState.byId = new Map();
  useInteraction.getState().clearPage();
});

// ── tests ────────────────────────────────────────────────────────────────────
describe("Corpus action declaration", () => {
  it("renders the shell dock primary 'Focus' and navigates on click", async () => {
    renderCorpus([{ id: 10, name: "A" }]);
    await screen.findByText("A");
    fireEvent.click(screen.getByTestId("dock-primary"));
    expect(navigateSpy).toHaveBeenCalledWith("/sample/10");
  });

  it("the sample stepper in the shell dock drives the cursor", async () => {
    renderCorpus([{ id: 10, name: "A" }, { id: 20, name: "B" }]);
    await screen.findByText("A");
    fireEvent.click(screen.getByTestId("dock-next-sample"));
    expect(screen.getByText("B").closest('[role="row"]')).toHaveAttribute("data-cursored", "true");
  });

  it("cull (x) dock button is enabled only in selection mode", async () => {
    renderCorpus([{ id: 10, name: "A" }]);
    const row = (await screen.findByText("A")).closest('[role="row"]')!;
    // browse mode: no Cull dock button (mode-gated action hidden when not enabled)
    expect(screen.queryByTestId("dock-action-cull")).toBeNull();
    fireEvent.click(within(row).getByRole("checkbox")); // → selection mode
    expect(screen.getByTestId("dock-action-cull")).toBeEnabled();
  });

  it("Enter (through the window keyboard layer) navigates to Focus", async () => {
    renderCorpus([{ id: 10, name: "A" }]);
    await screen.findByText("A");
    // Enter flows: window keydown → useKeyboardLayer → core("openFocus") →
    // sampleCursor.activate(). No row-level handler.
    fireEvent.keyDown(window, { key: "Enter" });
    expect(navigateSpy).toHaveBeenCalledWith("/sample/10");
  });

  it("l (openLoupe) through the keyboard layer navigates to the loupe", async () => {
    renderCorpus([{ id: 10, name: "A" }]);
    await screen.findByText("A");
    fireEvent.keyDown(window, { key: "l" });
    expect(navigateSpy).toHaveBeenCalledWith("/sample/10/loupe");
  });

  it("x (cull/drop) fires the batch verdict on the selected sample's frames", async () => {
    const byId = new Map<number, Exposure[]>([
      [10, [makeExposure({ id: 100, sample_id: 10 }), makeExposure({ id: 101, sample_id: 10 })]],
    ]);
    renderCorpus([{ id: 10, name: "A" }], { byId });
    const row = (await screen.findByText("A")).closest('[role="row"]')!;
    fireEvent.click(within(row).getByRole("checkbox")); // select → selection mode
    fireEvent.keyDown(window, { key: "x" });
    expect(batchMutate).toHaveBeenCalledWith({ sampleId: 10, exposureId: 100, status: "rejected" });
    expect(batchMutate).toHaveBeenCalledWith({ sampleId: 10, exposureId: 101, status: "rejected" });
  });

  it("k (keep) fires the accepted verdict on the selected sample's frames", async () => {
    const byId = new Map<number, Exposure[]>([[10, [makeExposure({ id: 100, sample_id: 10 })]]]);
    renderCorpus([{ id: 10, name: "A" }], { byId });
    const row = (await screen.findByText("A")).closest('[role="row"]')!;
    fireEvent.click(within(row).getByRole("checkbox"));
    fireEvent.keyDown(window, { key: "k" });
    expect(batchMutate).toHaveBeenCalledWith({ sampleId: 10, exposureId: 100, status: "accepted" });
  });

  it("r (restore) fires the null verdict on the selected sample's frames", async () => {
    const byId = new Map<number, Exposure[]>([[10, [makeExposure({ id: 100, sample_id: 10 })]]]);
    renderCorpus([{ id: 10, name: "A" }], { byId });
    const row = (await screen.findByText("A")).closest('[role="row"]')!;
    fireEvent.click(within(row).getByRole("checkbox"));
    fireEvent.keyDown(window, { key: "r" });
    expect(batchMutate).toHaveBeenCalledWith({ sampleId: 10, exposureId: 100, status: null });
  });

  it("x is inert in browse mode (mode-gated action not enabled)", async () => {
    const byId = new Map<number, Exposure[]>([[10, [makeExposure({ id: 100, sample_id: 10 })]]]);
    renderCorpus([{ id: 10, name: "A" }], { byId });
    await screen.findByText("A");
    fireEvent.keyDown(window, { key: "x" }); // no selection → browse → inert
    expect(batchMutate).not.toHaveBeenCalled();
  });
});

// ── shared cursor-contract harness (spec §10) ────────────────────────────────
// Proves the Corpus sample cursor satisfies click == arrow == stepper parity
// plus activate parity and selection orthogonality — the SAME contract the base
// useListCursor passes. Mounts the real hook over a few sample ids.
runCursorContract("Corpus sample cursor", () => {
  const onActivate = vi.fn();
  const IDS = [10, 20, 30];
  return {
    ui: (capture) => {
      function Probe(): JSX.Element {
        const cursor = useListCursor({
          ids: IDS,
          onActivate,
          stepperLabel: "Sample",
          stepperTestIdBase: "sample",
          axis: "vertical",
        });
        capture({ cursor, ids: IDS, onActivate });
        return (
          <div aria-multiselectable data-interaction-scope="">
            {IDS.map((id) => (
              <div key={id} {...cursor.rowProps(id)}>
                row {id}
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
