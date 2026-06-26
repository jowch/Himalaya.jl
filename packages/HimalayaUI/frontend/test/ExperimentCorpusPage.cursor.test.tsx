// test/ExperimentCorpusPage.cursor.test.tsx
//
// ID-based cursor tests for ExperimentCorpusPage: roving tabindex, click-to-park,
// and Arrow key navigation. Arrow nav now flows through the shell window keyboard
// layer (scope-exempt), so the harness mounts it via KbLayer; the event still
// bubbles from the grid to the window listener.
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { CorpusSample, Exposure } from "../src/api";
import { ExperimentCorpusPage } from "../src/print/pages/ExperimentCorpusPage";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";

function KbLayer({ children }: { children: React.ReactNode }): JSX.Element {
  useKeyboardLayer();
  return <>{children}</>;
}

// ── mock data plane ──────────────────────────────────────────────────────────
const cursorState = {
  samples: [] as CorpusSample[],
  byId: new Map<number, Exposure[]>(),
};

vi.mock("../src/queries", () => ({
  useExperiment: () => ({ data: undefined, isLoading: false }),
  useLoads: () => ({ data: [], isLoading: false }),
  useTriggerScan: () => ({ mutate: vi.fn() }),
  useCorpusSamples: () => ({ data: cursorState.samples, isLoading: false, isError: false }),
  useCorpusExposures: (filtered: CorpusSample[]) => {
    const byId = new Map<number, Exposure[]>();
    for (const s of filtered) {
      const exps = cursorState.byId.get(s.id);
      if (exps !== undefined) byId.set(s.id, exps);
    }
    return { byId, isLoading: false };
  },
  useSetExposureStatusBatch: () => ({ mutate: vi.fn() }),
  useSetExposureStatus: () => ({ mutate: vi.fn() }),
  useSelectExposure: () => ({ mutate: vi.fn() }),
}));

vi.mock("../src/state", () => ({
  useAppState: (sel: (s: Record<string, unknown>) => unknown) =>
    sel({ ingestInFlight: null, username: "tester" }),
}));

vi.mock("../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

// ── helpers ──────────────────────────────────────────────────────────────────
function makeSample(id: number, name: string): CorpusSample {
  return { id, experiment_id: 1, name, notes: null, q_units: "A-1", tags: [], phase: null } as CorpusSample;
}

function renderCorpus(samples: Array<{ id: number; name: string }>, expId = 1) {
  cursorState.samples = samples.map((s) => makeSample(s.id, s.name));
  cursorState.byId = new Map();
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[`/experiments/${expId}/corpus`]}>
        <KbLayer>
          <Routes>
            <Route path="/experiments/:id/corpus" element={<ExperimentCorpusPage />} />
            <Route path="/sample/:sampleId" element={<div data-testid="focus-route" />} />
          </Routes>
        </KbLayer>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

beforeEach(() => { cursorState.samples = []; cursorState.byId = new Map(); });

// ── tests ────────────────────────────────────────────────────────────────────
describe("Corpus cursor (ID-based, click-to-park)", () => {
  it("clicking a sample row parks the cursor there (data-cursored)", async () => {
    renderCorpus([{ id: 10, name: "A" }, { id: 20, name: "B" }, { id: 30, name: "C" }]);
    const rowB = await screen.findByText("B");
    fireEvent.click(rowB.closest('[role="row"]')!);
    expect(rowB.closest('[role="row"]')).toHaveAttribute("data-cursored", "true");
  });

  it("exactly one row is tabbable at a time (roving)", async () => {
    const { container } = renderCorpus([{ id: 10, name: "A" }, { id: 20, name: "B" }]);
    await screen.findByText("A");
    expect(container.querySelectorAll('[role="row"][tabindex="0"]').length).toBe(1);
  });

  it("ArrowDown on the grid moves the cursor to the next sample", async () => {
    renderCorpus([{ id: 10, name: "A" }, { id: 20, name: "B" }]);
    await screen.findByText("A");
    const scope = screen.getByTestId("experiment-corpus");
    const grid = within(scope).getByTestId("corpus-grid");
    fireEvent.keyDown(grid, { key: "ArrowDown" });
    expect(within(scope).getByText("B").closest('[role="row"]')).toHaveAttribute("data-cursored", "true");
  });
});
