// test/ExperimentCorpusPage.interaction.test.tsx
//
// Interaction coverage for ExperimentCorpusPage (the /experiments/:id/corpus
// home), ported from the retired SamplesPage unit suite when the /samples
// contact-sheet page was deleted. ExperimentCorpusPage carries the SAME
// cull/keep/compose/keyboard behaviour verbatim, scoped to one experiment via
// the :id route param.
//
// Ported from:
//   - test/SamplesPage.keyboard.test.tsx     → "keyboard cursor + cull verbs"
//   - test/SamplesPage.samplePicker.test.tsx → "sample-grain compose bar"
//   - test/SamplesPage.staleSelect.test.tsx  → "selection clears on scope change"
//
// Harness: mocks the `queries` hook layer (as the source suites did) so the
// batch mutator is a directly-assertable spy and exposures are present
// synchronously (no async React-Query settle before fireEvent). The route
// param drives `expId`; scoped samples are pre-filtered by experiment_id, so
// the mocked corpus lists exactly this experiment's rows.
//
// INTENTIONALLY DROPPED (SamplesPage-only chrome, not present in
// ExperimentCorpusPage — the corpus is already experiment-scoped):
//   - the contact-sheet "N / M samples screened" header / useScreenedProgress.
//   - the ?experiment search-param filter + unknown-experiment EmptyState +
//     "experiment N" title fallback (resolved by the :id route param instead).
//   - URL sort key/dir round-trip (no sortable columns on this surface).
//   - the page-title heading ("The corpus" / experiment name).
// These have no analogue here, so their assertions are not ported.
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter, Routes, Route, useNavigate } from "react-router-dom";
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
        <Routes>
          <Route path="/experiments/:id/corpus" element={<ExperimentCorpusPage />} />
          <Route path="/sample/:sampleId/loupe" element={<div data-testid="loupe-route" />} />
          <Route path="/sample/:sampleId" element={<div data-testid="focus-route" />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

beforeEach(() => {
  batchMutate.mockClear();
});

// ── grid Arrow navigation (cursor movement only; no keyboard layer here) ─────
// This harness mounts the page WITHOUT the shell keyboard layer, so it proves
// the grid container's own ↑/↓/←/→ onKeyDown drives the cursor / frame axis.
// Enter→navigate flows through the window keyboard layer and is proven in
// ExperimentCorpusPage.actions.test.tsx (which mounts useKeyboardLayer).
describe("ExperimentCorpusPage — grid Arrow navigation", () => {
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
});

// ── sample-grain compose segment (dock-as-action-bar, items 2/5) ─────────────
// The floating ComposeBar/CullBar were folded into the Dock: a sample-grain
// `dock-compose` segment appears on check; a frame-grain `dock-selection-count`
// readout appears on exposure select. Both are absent (not just hidden/inert)
// when their selection is empty.
describe("ExperimentCorpusPage — sample-grain picker (dock compose segment)", () => {
  function seedPicker(): void {
    // Two scoped samples, no exposures (compose segment is sample-grain only).
    state.samples = [
      corpus({ id: 10, experiment_id: 1, name: "A" }),
      corpus({ id: 11, experiment_id: 1, name: "B" }),
    ];
    state.byId = new Map();
  }
  beforeEach(seedPicker);

  it("the compose segment is absent when no samples are checked", () => {
    renderAt(1);
    expect(screen.queryByTestId("dock-compose")).toBeNull();
  });

  it("the compose segment appears (count 1) when a sample is checked", async () => {
    const user = userEvent.setup();
    renderAt(1);
    const checkboxes = screen.getAllByRole("checkbox");
    await user.click(checkboxes[0]!);
    expect(screen.getByTestId("dock-compose")).toHaveTextContent("1");
  });

  it("Clear resets the sample selection and removes the compose segment", async () => {
    const user = userEvent.setup();
    renderAt(1);
    await user.click(screen.getAllByRole("checkbox")[0]!);
    await user.click(
      within(screen.getByTestId("dock-compose")).getByRole("button", { name: /clear/i }),
    );
    expect(screen.queryByTestId("dock-compose")).toBeNull();
  });

  it("the frame-grain selection readout stays independent of sample-grain checks", async () => {
    const user = userEvent.setup();
    renderAt(1);
    await user.click(screen.getAllByRole("checkbox")[0]!);
    // No exposure selected → no frame-grain selection readout in the dock.
    expect(screen.queryByTestId("dock-selection-count")).toBeNull();
  });
});

// ── stale-selection reset on scope change (ported from staleSelect.test) ─────
describe("ExperimentCorpusPage — selection clears on experiment-scope change", () => {
  beforeEach(() => {
    // Two samples in DISTINCT experiments. The corpus is filtered by the route's
    // :id, so re-rendering at a different /experiments/:id/corpus re-scopes the
    // visible list — exactly the SamplesPage ?beamtime swap, now a route param.
    state.samples = [
      corpus({ id: 10, experiment_id: 1, name: "A" }),
      corpus({ id: 11, experiment_id: 2, name: "B" }),
    ];
    state.byId = new Map();
  });

  // A sibling control that re-scopes the page IN PLACE (same router, new :id)
  // the way the shell's experiment switcher does — a real navigation, so the
  // page's expId changes without remounting and the reset effect fires.
  function ScopeSwitch(): JSX.Element {
    const navigate = useNavigate();
    return (
      <button data-testid="scope-switch" onClick={() => navigate("/experiments/2/corpus")}>
        scope to exp 2
      </button>
    );
  }

  it("a sample-grain check does NOT survive a change of experiment scope", async () => {
    const user = userEvent.setup();
    const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
    render(
      <QueryClientProvider client={qc}>
        <MemoryRouter initialEntries={["/experiments/1/corpus"]}>
          <ScopeSwitch />
          <Routes>
            <Route path="/experiments/:id/corpus" element={<ExperimentCorpusPage />} />
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );

    // Check the only experiment-1 sample → the dock compose segment carries one pick.
    await user.click(screen.getAllByRole("checkbox")[0]!);
    expect(screen.getByTestId("dock-compose")).toHaveTextContent("1");

    // Re-scope to experiment 2 in place. The page's expId-change effect must
    // clear both selection grains, so the carry does not linger pointing at an
    // off-scope sample (no row remains to uncheck it).
    await user.click(screen.getByTestId("scope-switch"));
    expect(screen.queryByTestId("dock-compose")).toBeNull();
  });
});
