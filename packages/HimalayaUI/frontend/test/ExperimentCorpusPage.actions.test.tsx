// test/ExperimentCorpusPage.actions.test.tsx
//
// Corpus page action-declaration tests. These mount the full TestShell
// (keyboard layer + InteractionDock) so dock-primary and keyboard-triggered
// actions work correctly. Navigate is mocked as a spy so we can assert on
// programmatic navigation without actual routing.
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

// ── mock data plane ──────────────────────────────────────────────────────────
const actionsState = {
  samples: [] as CorpusSample[],
  byId: new Map<number, Exposure[]>(),
};

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
  useSetExposureStatusBatch: () => ({ mutate: vi.fn() }),
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

function renderCorpus(samples: Array<{ id: number; name: string }>, expId = 1) {
  actionsState.samples = samples.map((s) => makeSample(s.id, s.name));
  actionsState.byId = new Map();
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[`/experiments/${expId}/corpus`]}>
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

  it("cull (x) dock button only appears in selection mode", async () => {
    renderCorpus([{ id: 10, name: "A" }]);
    const row = (await screen.findByText("A")).closest('[role="row"]')!;
    // browse mode: no Cull dock button (mode-gated action hidden when not enabled)
    expect(screen.queryByTestId("dock-action-cull")).toBeNull();
    fireEvent.click(within(row).getByRole("checkbox"));
    expect(screen.getByTestId("dock-action-cull")).toBeInTheDocument();
  });

  it("Enter on a focused row navigates to Focus", async () => {
    renderCorpus([{ id: 10, name: "A" }]);
    const row = (await screen.findByText("A")).closest('[role="row"]')! as HTMLElement;
    fireEvent.keyDown(row, { key: "Enter" });
    expect(navigateSpy).toHaveBeenCalledWith("/sample/10");
  });
});
