import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { CorpusSample, Exposure } from "../src/api";

// ── mutator spy + nav target probe ──────────────────────────────────────────
const batchMutate = vi.fn();

// ── mock data plane (mutated per test) ──────────────────────────────────────
const state = {
  samples: [] as CorpusSample[],
  byId: new Map<number, Exposure[]>(),
};

// Mock the queries hook layer. ExperimentCorpusPage reads:
//   useExperiment, useLoads, useTriggerScan, useCorpusSamples,
//   useCorpusExposures, useSetExposureStatusBatch.
vi.mock("../src/queries", () => ({
  useExperiment: () => ({ data: undefined, isLoading: false }),
  useLoads: () => ({ data: [], isLoading: false }), // empty → no review banner
  useTriggerScan: () => ({ mutate: vi.fn() }),
  useCorpusSamples: () => ({ data: state.samples, isLoading: false, isError: false }),
  useCorpusExposures: (filtered: CorpusSample[]) => {
    // Only the rows for the filtered (scoped) samples — mirrors the real hook.
    const byId = new Map<number, Exposure[]>();
    for (const s of filtered) {
      const exps = state.byId.get(s.id);
      if (exps !== undefined) byId.set(s.id, exps);
    }
    return { byId, isLoading: false };
  },
  useSetExposureStatusBatch: () => ({ mutate: batchMutate }),
  useSetExposureStatus: () => ({ mutate: vi.fn() }),
  useSelectExposure: () => ({ mutate: vi.fn() }),
}));

// Zustand store: ingestInFlight stays null so no takeover state pre-empts the
// sheet/dock. The page calls useAppState((s) => s.ingestInFlight?.[expId]).
vi.mock("../src/state", () => ({
  useAppState: (sel: (s: Record<string, unknown>) => unknown) =>
    sel({ ingestInFlight: null, username: "tester" }),
}));

// DetectorImage touches fetch / createImageBitmap (absent in JSDOM).
vi.mock("../src/print/detector/DetectorImage", () => ({
  DetectorImage: () => <div data-testid="mock-detector-image" />,
}));

import { ExperimentCorpusPage } from "../src/print/pages/ExperimentCorpusPage";
import { useKeyboardLayer } from "../src/print/interaction/useKeyboardLayer";

/** Mounts the shell window keyboard layer, the way AppShell does — arrow nav now
 *  flows through it (scope-exempt), not a container onKeyDown. */
function KbLayer({ children }: { children: React.ReactNode }): JSX.Element {
  useKeyboardLayer();
  return <>{children}</>;
}

function Landed(): JSX.Element {
  const loc = useLocation();
  return <div data-testid="landed">{loc.pathname + loc.search}</div>;
}

// ── fixtures ────────────────────────────────────────────────────────────────
function corpus(over: Partial<CorpusSample> = {}): CorpusSample {
  return {
    id: 1, experiment_id: 1, name: "JC001",
    notes: null, q_units: "A-1", tags: [], phase: null, ...over,
  } as CorpusSample;
}
function exp(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 100, sample_id: 1, filename: "JC001-001.dat", kind: "file",
    selected: false, status: "accepted", image_path: "/x.tif", image_version: "v1",
    tags: [], sources: [], trace_hash: "h", analysis_inputs_hash: "a", ...over,
  };
}

/** Render the corpus page at /experiments/:id/corpus. The route param is the
 *  experiment scope; the mocked corpus is filtered by experiment_id to it. */
function renderAt(expId = 1) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[`/experiments/${expId}/corpus`]}>
        <KbLayer>
          <Routes>
            <Route path="/experiments/:id/corpus" element={<ExperimentCorpusPage />} />
            <Route path="/experiments/:id/workspace/:sampleId" element={<Landed />} />
          </Routes>
        </KbLayer>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

beforeEach(() => {
  batchMutate.mockClear();
});

// ── grid Arrow navigation (cursor movement only; no keyboard layer here) ─────
// Arrow nav flows through the shell window keyboard layer (mounted via KbLayer),
// SCOPE-EXEMPT: ↑/↓/←/→ drive the cursor / frame axis from anywhere, so the page
// never scrolls when focus sits outside the grid. The last case fires on
// document.body to pin that global reach.
describe("ExperimentCorpusPage — Arrow navigation", () => {
  // Three scoped samples (experiment_id=1) so the cursor can step 0→1→2; S0
  // carries two frames so the ←/→ frame axis has a target.
  const S0 = corpus({ id: 10, experiment_id: 1, name: "Sample 10" });
  const S1 = corpus({ id: 11, experiment_id: 1, name: "Sample 11" });
  const S2 = corpus({ id: 12, experiment_id: 1, name: "Sample 12" });

  function seedKeyboard(): void {
    state.samples = [S0, S1, S2];
    state.byId = new Map<number, Exposure[]>([
      [10, [exp({ id: 100, sample_id: 10 }), exp({ id: 101, sample_id: 10 })]],
      [11, []],
      [12, []],
    ]);
  }

  beforeEach(seedKeyboard);

  function cursoredRowName(): string | null {
    const row = document.querySelector('[data-cursored="true"]');
    return row?.textContent?.includes("Sample 12") ? "Sample 12"
      : row?.textContent?.includes("Sample 11") ? "Sample 11"
      : row?.textContent?.includes("Sample 10") ? "Sample 10" : null;
  }

  it("ArrowDown on the grid moves the sample cursor forward (by id)", () => {
    renderAt(1);
    const grid = screen.getByTestId("corpus-grid");
    fireEvent.keyDown(grid, { key: "ArrowDown" });
    fireEvent.keyDown(grid, { key: "ArrowDown" });
    expect(cursoredRowName()).toBe("Sample 12");
  });

  it("ArrowDown clamps at the last sample (not circular)", () => {
    renderAt(1);
    const grid = screen.getByTestId("corpus-grid");
    fireEvent.keyDown(grid, { key: "ArrowDown" }); // → S1
    fireEvent.keyDown(grid, { key: "ArrowDown" }); // → S2
    fireEvent.keyDown(grid, { key: "ArrowDown" }); // stays S2
    expect(cursoredRowName()).toBe("Sample 12");
  });

  it("ArrowUp moves the cursor back", () => {
    renderAt(1);
    const grid = screen.getByTestId("corpus-grid");
    fireEvent.keyDown(grid, { key: "ArrowDown" }); // → S1
    fireEvent.keyDown(grid, { key: "ArrowDown" }); // → S2
    fireEvent.keyDown(grid, { key: "ArrowUp" });   // → S1
    expect(cursoredRowName()).toBe("Sample 11");
  });

  it("ArrowRight walks the active sample's frame axis (cursored thumbnail moves)", () => {
    renderAt(1);
    const grid = screen.getByTestId("corpus-grid");
    // Cursor on S0 (id 10), which has two frames. Frame 0 is cursored at start.
    const thumbsBefore = screen.getAllByTestId("thumbnail");
    expect(thumbsBefore[0]!.getAttribute("data-state")).toContain("cursored");
    fireEvent.keyDown(grid, { key: "ArrowRight" }); // frame 0 → 1
    const thumbsAfter = screen.getAllByTestId("thumbnail");
    expect(thumbsAfter[1]!.getAttribute("data-state")).toContain("cursored");
  });

  it("arrows drive the cursor even when focus is OUTSIDE the grid (on body)", () => {
    // The regression guard for the global-arrow change: with focus on
    // document.body (never the scope), ArrowDown still moves the surface cursor
    // instead of scrolling the page.
    renderAt(1);
    fireEvent.keyDown(document.body, { key: "ArrowDown" });
    expect(cursoredRowName()).toBe("Sample 11");
  });
});

describe("ExperimentCorpusPage — opening a frame", () => {
  it("double-clicking a thumbnail opens the workspace on that frame, not the representative", () => {
    state.samples = [corpus({ id: 10, experiment_id: 1, name: "Sample 10" })];
    state.byId = new Map<number, Exposure[]>([
      [10, [exp({ id: 100, sample_id: 10, selected: true }), exp({ id: 101, sample_id: 10 })]],
    ]);
    renderAt(1);
    fireEvent.doubleClick(screen.getAllByTestId("thumbnail")[1]!);
    expect(screen.getByTestId("landed")).toHaveTextContent("/experiments/1/workspace/10?exposure=101");
  });
});
