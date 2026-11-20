import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import type { CorpusSample, Exposure, IndexEntry, Assignment } from "../../src/api";
import { WorkspacePage } from "../../src/workspace/WorkspacePage";

const spies = {
  add: vi.fn(),
  remove: vi.fn(),
  select: vi.fn(),
  status: vi.fn(),
};

const state = {
  corpus: [] as CorpusSample[],
  exposures: [] as Exposure[],
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
};

vi.mock("../../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.corpus, isLoading: false }),
  useExperiment: () => ({ data: { id: 1, name: "BL-19" } }),
  useExposures: () => ({ data: state.exposures }),
  useTrace: () => ({ data: undefined }),
  usePeaks: () => ({ data: [] }),
  useIndices: () => ({ data: state.indices }),
  useAssignment: () => ({ data: state.assignment }),
  useAddAssignmentPhase: () => ({ mutate: spies.add }),
  useRemoveAssignmentPhase: () => ({ mutate: spies.remove }),
  useSelectExposure: () => ({ mutate: spies.select }),
  useSetExposureStatus: () => ({ mutate: spies.status }),
}));

function sample(id: number, name: string, extra: Partial<CorpusSample> = {}): CorpusSample {
  return { id, experiment_id: 1, name, notes: null, tags: [], q_units: "A-1", ...extra } as CorpusSample;
}

function exposure(id: number, extra: Partial<Exposure> = {}): Exposure {
  return {
    id, sample_id: 10, filename: `f${id}.dat`, kind: "file", selected: false, status: null,
    image_path: null, image_version: "", tags: [], sources: [], ...extra,
  } as Exposure;
}

function index(id: number, phase: string, score: number): IndexEntry {
  return {
    id, exposure_id: 101, phase, basis: 0.1, score, r_squared: 0.99, lattice_d: 100,
    ngc: null, status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [0.1, 0.14],
  } as IndexEntry;
}

function Where(): JSX.Element {
  return <div data-testid="where">{useLocation().pathname}</div>;
}

function renderAt(path: string): void {
  render(
    <MemoryRouter initialEntries={[path]}>
      <Routes>
        <Route path="/experiments/:id/workspace/:sampleId?" element={<><WorkspacePage /><Where /></>} />
      </Routes>
    </MemoryRouter>,
  );
}

beforeEach(() => {
  Element.prototype.scrollIntoView = vi.fn();
  Object.values(spies).forEach((s) => s.mockReset());
  state.corpus = [
    sample(10, "A1", { assignment_state: "indexed", phase: "Pn3m" }),
    sample(11, "A2"),
    sample(12, "A3", { assignment_state: "form_factor" }),
    sample(13, "A4"),
    { ...sample(99, "other"), experiment_id: 2 },
  ];
  state.exposures = [exposure(100), exposure(101, { selected: true })];
  state.indices = [index(1, "Im3m", 0.4), index(2, "Pn3m", 0.9)];
  state.assignment = { exposure_id: 101, state: "indexed", members: [2] };
});

describe("WorkspacePage", () => {
  it("lists only this experiment's samples and marks unreviewed ones", () => {
    renderAt("/experiments/1/workspace/10");
    const list = screen.getByTestId("workspace-samples");
    expect(within(list).queryByText("other")).toBeNull();
    expect(within(list).getAllByLabelText("Unreviewed")).toHaveLength(2);
    expect(screen.getByText("2 of 4 reviewed")).toBeInTheDocument();
  });

  it("J jumps to the next unreviewed sample and Shift+J goes back", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "j" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
    fireEvent.keyDown(window, { key: "J", shiftKey: true });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/13");
  });

  it("arrow down steps to the next sample in the list", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
  });

  it("ranks candidates by score and number keys toggle them", () => {
    renderAt("/experiments/1/workspace/10");
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows[0]).toHaveTextContent("Pn3m");
    expect(rows[0]).toHaveAttribute("aria-pressed", "true");
    fireEvent.keyDown(window, { key: "1" });
    expect(spies.remove).toHaveBeenCalledWith(2);
    fireEvent.keyDown(window, { key: "2" });
    expect(spies.add).toHaveBeenCalledWith(1);
    expect(spies.select).not.toHaveBeenCalled();
  });

  it("selects the exposure before adding a phase to an unselected one", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    fireEvent.keyDown(window, { key: "2" });
    expect(spies.select).toHaveBeenCalledWith(100);
    expect(spies.add).toHaveBeenCalledWith(1);
  });

  it("X drops the current exposure", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "x" });
    expect(spies.status).toHaveBeenCalledWith({ exposureId: 101, status: "rejected" });
  });

  it("ignores shortcuts while typing in the filter", () => {
    renderAt("/experiments/1/workspace/10");
    const input = screen.getByLabelText("Filter samples");
    fireEvent.change(input, { target: { value: "A" } });
    fireEvent.keyDown(input, { key: "1" });
    expect(spies.remove).not.toHaveBeenCalled();
  });

  it("filters the sample list by name", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.change(screen.getByLabelText("Filter samples"), { target: { value: "a3" } });
    const list = screen.getByTestId("workspace-samples");
    expect(within(list).getAllByRole("button")).toHaveLength(1);
    expect(within(list).getByText("A3")).toBeInTheDocument();
  });
});
