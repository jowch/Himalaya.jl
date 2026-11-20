import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import type { CorpusSample, Exposure, IndexEntry, Assignment } from "../../src/api";
import { WorkspacePage } from "../../src/workspace/WorkspacePage";

type Call = { hook: string; boundTo: number; arg: unknown; callbacks?: { onSuccess?: () => void } };
const calls: Call[] = [];

const state = {
  corpus: [] as CorpusSample[],
  exposures: [] as Exposure[],
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
  username: "jc" as string | undefined,
};

function recorder(hook: string) {
  return (boundTo: number) => ({
    mutate: (arg: unknown, callbacks?: { onSuccess?: () => void }) => {
      calls.push({ hook, boundTo, arg, ...(callbacks ? { callbacks } : {}) });
    },
  });
}

vi.mock("../../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.corpus, isLoading: false }),
  useExperiment: () => ({ data: { id: 1, name: "BL-19" } }),
  useExposures: () => ({ data: state.exposures }),
  useTrace: () => ({ data: undefined }),
  usePeaks: () => ({ data: [] }),
  useIndices: () => ({ data: state.indices }),
  useAssignment: () => ({ data: state.assignment }),
  useAddAssignmentPhase: recorder("add"),
  useRemoveAssignmentPhase: recorder("remove"),
  useSelectExposure: recorder("select"),
  useSetExposureStatus: recorder("status"),
}));

vi.mock("../../src/state", () => ({
  useAppState: (sel: (s: { username: string | undefined }) => unknown) => sel({ username: state.username }),
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

function page(path: string): JSX.Element {
  return (
    <MemoryRouter initialEntries={[path]}>
      <Routes>
        <Route path="/experiments/:id/workspace/:sampleId?" element={<><WorkspacePage /><Where /></>} />
      </Routes>
    </MemoryRouter>
  );
}

function renderAt(path: string) {
  return render(page(path));
}

function currentFrame(): string | null {
  const strip = screen.getByTestId("workspace-exposures");
  return strip.querySelector("[aria-current='true']")?.getAttribute("aria-label") ?? null;
}

beforeEach(() => {
  Element.prototype.scrollIntoView = vi.fn();
  calls.length = 0;
  state.username = "jc";
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
    expect(within(list).getAllByRole("img", { name: "Unreviewed" })).toHaveLength(2);
    expect(screen.getByTestId("workspace-progress")).toHaveTextContent("2 of 4 reviewed");
  });

  it("updates the list and the count when the corpus status changes", () => {
    const { rerender } = renderAt("/experiments/1/workspace/11");
    state.corpus = state.corpus.map((s) => (s.id === 11 ? { ...s, assignment_state: "indexed", phase: "Im3m" } : s));
    rerender(page("/experiments/1/workspace/11"));
    const row = within(screen.getByTestId("workspace-samples")).getByRole("button", { current: true });
    expect(within(row).getByTestId("workspace-status")).toHaveTextContent("Im3m");
    expect(screen.getByTestId("workspace-progress")).toHaveTextContent("3 of 4 reviewed");
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

  it("opens on the representative frame", () => {
    renderAt("/experiments/1/workspace/10");
    expect(currentFrame()).toBe("Frame 2");
  });

  it("ranks candidates by score and number keys toggle them on the shown exposure", () => {
    renderAt("/experiments/1/workspace/10");
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows[0]).toHaveTextContent("Pn3m");
    expect(rows[0]).toHaveAttribute("aria-pressed", "true");
    fireEvent.keyDown(window, { key: "1" });
    fireEvent.keyDown(window, { key: "2" });
    expect(calls).toEqual([
      { hook: "remove", boundTo: 101, arg: 2 },
      { hook: "add", boundTo: 101, arg: 1 },
    ]);
  });

  it("on an unselected frame, adds the phase only after the select succeeds", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    fireEvent.keyDown(window, { key: "2" });
    expect(calls).toHaveLength(1);
    expect(calls[0]).toMatchObject({ hook: "select", boundTo: 10, arg: 100 });
    calls[0]!.callbacks!.onSuccess!();
    expect(calls[1]).toEqual({ hook: "add", boundTo: 100, arg: 1 });
  });

  it("X drops the shown frame and keeps it in view, so a second X restores it", () => {
    const { rerender } = renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "x" });
    expect(calls[0]).toEqual({ hook: "status", boundTo: 10, arg: { exposureId: 101, status: "rejected" } });
    state.exposures = [exposure(100), exposure(101, { selected: true, status: "rejected" })];
    rerender(page("/experiments/1/workspace/10"));
    expect(currentFrame()).toBe("Frame 2, dropped");
    fireEvent.keyDown(window, { key: "x" });
    expect(calls[1]).toEqual({ hook: "status", boundTo: 10, arg: { exposureId: 101, status: null } });
  });

  it("refuses a phase call on a dropped frame and says why", () => {
    state.exposures = [exposure(100), exposure(101, { selected: true, status: "rejected" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(currentFrame()).toBe("Frame 2, dropped");
    fireEvent.keyDown(window, { key: "2" });
    expect(calls).toEqual([]);
    expect(screen.getByText("This frame is dropped. Restore it to call a phase.")).toBeInTheDocument();
  });

  it("ignores shortcuts while typing in the filter", () => {
    renderAt("/experiments/1/workspace/10");
    const input = screen.getByLabelText("Filter samples");
    fireEvent.change(input, { target: { value: "A" } });
    fireEvent.keyDown(input, { key: "1" });
    expect(calls).toEqual([]);
  });

  it("ignores shortcuts from inside a dialog", () => {
    renderAt("/experiments/1/workspace/10");
    const dialog = document.createElement("div");
    dialog.setAttribute("role", "dialog");
    const button = document.createElement("button");
    dialog.appendChild(button);
    document.body.appendChild(dialog);
    fireEvent.keyDown(button, { key: "x" });
    fireEvent.keyDown(button, { key: "j" });
    expect(calls).toEqual([]);
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/10");
    dialog.remove();
  });

  it("does not write anything before a user is set", () => {
    state.username = undefined;
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "x" });
    fireEvent.keyDown(window, { key: "1" });
    expect(calls).toEqual([]);
  });

  it("filters the sample list by name", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.change(screen.getByLabelText("Filter samples"), { target: { value: "a3" } });
    const list = screen.getByTestId("workspace-samples");
    expect(within(list).getAllByRole("button")).toHaveLength(1);
    expect(within(list).getByText("A3")).toBeInTheDocument();
  });

  it("says when the sample in the link is not in this experiment", () => {
    renderAt("/experiments/1/workspace/99");
    expect(screen.getByText("This sample is not in this experiment.")).toBeInTheDocument();
  });

  it("says when the link is not valid", () => {
    renderAt("/experiments/abc/workspace");
    expect(screen.getByText("This workspace link is not valid.")).toBeInTheDocument();
  });
});
