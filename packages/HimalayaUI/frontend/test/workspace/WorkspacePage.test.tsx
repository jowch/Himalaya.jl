import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import type { CorpusSample, Exposure, IndexEntry, Assignment, Peak, Trace } from "../../src/api";
import { WorkspacePage } from "../../src/workspace/WorkspacePage";

type Call = { hook: string; boundTo: number; arg: unknown; callbacks?: { onSuccess?: () => void } };
const calls: Call[] = [];

const state = {
  corpus: [] as CorpusSample[],
  exposures: [] as Exposure[],
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
  username: "jc" as string | undefined,
  trace: undefined as Trace | undefined,
  peaks: [] as Peak[],
  compareTraces: new Map<number, Trace>(),
  compareExposures: new Map<number, Exposure>(),
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
  useTrace: () => ({ data: state.trace }),
  usePeaks: () => ({ data: state.peaks }),
  useIndices: () => ({ data: state.indices }),
  useAssignment: () => ({ data: state.assignment }),
  useAddAssignmentPhase: recorder("add"),
  useRemoveAssignmentPhase: recorder("remove"),
  useSelectExposure: recorder("select"),
  useSetExposureStatus: recorder("status"),
  useAddPeak: recorder("addPeak"),
  useRemovePeak: recorder("removePeak"),
  useSetPeakExcluded: recorder("excludePeak"),
}));

vi.mock("../../src/workspace/useCompareTraces", () => ({
  useCompareTraces: (ids: readonly number[]) =>
    new Map(ids.flatMap((id) => (state.compareTraces.has(id)
      ? [[id, { trace: state.compareTraces.get(id)!, exposure: state.compareExposures.get(id) ?? exposure(500 + id) }] as const]
      : []))),
}));

function makeTrace(level: number): Trace {
  const q = Array.from({ length: 21 }, (_, i) => 0.05 + i * 0.01);
  return { q, I: q.map((x) => level * (1 + 10 * Math.exp(-((x - 0.1) ** 2) / 1e-5))), sigma: q.map(() => 0) };
}

function peak(id: number, q: number, source: "auto" | "manual", excluded = false): Peak {
  return { id, exposure_id: 101, q, intensity: null, prominence: null, sharpness: null, source, excluded };
}

function clickPlotAt(q: number, detail = 1): void {
  const plotWidth = 640 - 60 - 24;
  const px = 60 + ((q - 0.05) / 0.2) * plotWidth;
  fireEvent.click(screen.getByTestId("workspace-trace"), { clientX: px, clientY: 100, detail });
}

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
  state.trace = undefined;
  state.peaks = [];
  state.compareTraces = new Map();
  state.compareExposures = new Map();
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
  it("a click near an auto peak disables it, near a manual peak removes it, and elsewhere adds one", () => {
    state.trace = makeTrace(1);
    state.peaks = [peak(7, 0.1, "auto"), peak(8, 0.2, "manual")];
    renderAt("/experiments/1/workspace/10");
    clickPlotAt(0.1);
    clickPlotAt(0.2);
    clickPlotAt(0.15);
    expect(calls.slice(0, 2)).toEqual([
      { hook: "excludePeak", boundTo: 101, arg: { peakId: 7, excluded: true } },
      { hook: "removePeak", boundTo: 101, arg: 8 },
    ]);
    expect(calls[2]).toMatchObject({ hook: "addPeak", boundTo: 101 });
    expect(calls[2]!.arg as number).toBeCloseTo(0.15, 2);
  });

  it("a double-click adds one peak, not two", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    clickPlotAt(0.15, 1);
    clickPlotAt(0.15, 2);
    expect(calls.filter((c) => c.hook === "addPeak")).toHaveLength(1);
  });

  it("a click on a disabled auto peak restores it", () => {
    state.trace = makeTrace(1);
    state.peaks = [peak(7, 0.1, "auto", true)];
    renderAt("/experiments/1/workspace/10");
    clickPlotAt(0.1);
    expect(calls).toEqual([{ hook: "excludePeak", boundTo: 101, arg: { peakId: 7, excluded: false } }]);
  });

  it("does not edit peaks before a user is set", () => {
    state.trace = makeTrace(1);
    state.username = undefined;
    renderAt("/experiments/1/workspace/10");
    clickPlotAt(0.15);
    expect(calls).toEqual([]);
  });

  it("does not edit peaks on a dropped frame", () => {
    state.trace = makeTrace(1);
    state.exposures = [exposure(100), exposure(101, { selected: true, status: "rejected" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(currentFrame()).toBe("Frame 2, dropped");
    clickPlotAt(0.15);
    expect(calls).toEqual([]);
  });

  it("L switches the intensity axis between log and linear", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    const ticks = (): string =>
      Array.from(screen.getByTestId("workspace-trace").querySelectorAll('[data-role="axis-left"] text'))
        .map((t) => t.textContent)
        .join(" ");
    const logTicks = ticks();
    expect(screen.getByRole("button", { name: /Log/ })).toHaveAttribute("aria-pressed", "true");
    fireEvent.keyDown(window, { key: "l" });
    expect(screen.getByRole("button", { name: /Linear/ })).toHaveAttribute("aria-pressed", "false");
    expect(ticks()).not.toBe(logTicks);
  });

  it("C opens a picker, Enter adds the sample scaled to the current trace and subtracts it", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    const picker = screen.getByRole("dialog", { name: "Compare with a sample" });
    expect(within(picker).queryByText("A1")).toBeNull();
    fireEvent.change(within(picker).getByLabelText("Find a sample to compare"), { target: { value: "a2" } });
    fireEvent.keyDown(within(picker).getByLabelText("Find a sample to compare"), { key: "Enter" });
    expect(screen.queryByRole("dialog")).toBeNull();
    const chip = screen.getByRole("button", { name: "Remove A2" }).parentElement!;
    expect(chip).toHaveTextContent("×0.500");
    expect(screen.getByRole("button", { name: "A2" })).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByText("minus A2")).toBeInTheDocument();
    expect(screen.getByTestId("workspace-compares").querySelector("[data-compare='11']")?.getAttribute("d")).toMatch(/^M/);
  });

  it("the picker only adds: samples already compared are not offered", () => {
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    fireEvent.keyDown(window, { key: "c" });
    expect(within(screen.getByRole("dialog")).queryByText("A2")).toBeNull();
    fireEvent.keyDown(screen.getByLabelText("Find a sample to compare"), { key: "Escape" });
    expect(screen.getByRole("button", { name: "Remove A2" })).toBeInTheDocument();
  });

  it("Escape and an outside click close the picker without adding anything", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.keyDown(screen.getByLabelText("Find a sample to compare"), { key: "Escape" });
    expect(screen.queryByRole("dialog")).toBeNull();
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.mouseDown(screen.getByTestId("workspace-trace"));
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(screen.queryByRole("button", { name: /^Remove/ })).toBeNull();
  });

  it("subtracts with the scale that matches the high-q background, ignoring the sample's own peaks there", () => {
    const q = Array.from({ length: 101 }, (_, i) => 0.05 + i * 0.002);
    const bg = (x: number): number => 10 + 1 / x;
    const peakQ = [0.232, 0.244];
    state.trace = { q, I: q.map((x) => 2 * bg(x) + peakQ.reduce((a, p) => a + 500 * Math.exp(-((x - p) ** 2) / 4e-6), 0)), sigma: q.map(() => 0) };
    state.peaks = peakQ.map((p, i) => peak(7 + i, p, "auto"));
    state.compareTraces = new Map([[11, { q, I: q.map(bg), sigma: q.map(() => 0) }]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    expect(screen.getByRole("button", { name: "Remove A2" }).parentElement!).toHaveTextContent("×2.00");
  });

  it("scales by the transmitted beam when both frames carry monitors, even where high-q matching would disagree", () => {
    const monitors = (transmitted_counts: number) => ({ transmitted_counts, exposure_time: 15, counting_time: 15.8 });
    state.exposures = [exposure(100), exposure(101, { selected: true, ...monitors(6000) })];
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    state.compareExposures = new Map([[11, exposure(111, { sample_id: 11, ...monitors(8000) })]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    const scale = screen.getByRole("button", { name: "Remove A2" }).parentElement!.querySelector("[data-scale-source]")!;
    expect(scale).toHaveAttribute("data-scale-source", "monitor");
    expect(scale).toHaveTextContent(/^×0\.750$/);
  });

  it("says when a comparison had no monitors and was matched at high q instead", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    const scale = screen.getByRole("button", { name: "Remove A2" }).parentElement!.querySelector("[data-scale-source]")!;
    expect(scale).toHaveAttribute("data-scale-source", "matched");
    expect(scale).toHaveTextContent("×0.500 matched");
  });

  it("clicking the subtracted chip stops subtracting, R hides comparisons, and × removes one", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    fireEvent.click(screen.getByRole("button", { name: "A2" }));
    expect(screen.queryByTestId("workspace-subtraction")).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "A2" }));
    expect(screen.getByTestId("workspace-subtraction")).toBeInTheDocument();
    fireEvent.keyDown(window, { key: "r" });
    expect(screen.queryByTestId("workspace-subtraction")).toBeNull();
    expect(screen.getByTestId("workspace-compares").children).toHaveLength(0);
    fireEvent.keyDown(window, { key: "r" });
    fireEvent.click(screen.getByRole("button", { name: "Remove A2" }));
    expect(screen.queryByTestId("workspace-subtraction")).toBeNull();
    expect(screen.getByTestId("workspace-compares").children).toHaveLength(0);
  });

  it("re-adding a removed comparison leaves the current subtraction alone", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)], [12, makeTrace(3)]]);
    renderAt("/experiments/1/workspace/10");
    const pick = (id: number): void => {
      fireEvent.keyDown(window, { key: "c" });
      fireEvent.click(screen.getByRole("dialog").querySelector(`[data-sample-id='${id}']`)!);
    };
    pick(11);
    pick(12);
    fireEvent.click(screen.getByRole("button", { name: "A3" }));
    fireEvent.click(screen.getByRole("button", { name: "Remove A2" }));
    pick(11);
    expect(screen.getByRole("button", { name: "A2" })).toHaveAttribute("aria-pressed", "false");
    expect(screen.getByRole("button", { name: "A3" })).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByText("minus A3")).toBeInTheDocument();
  });

  it("keeps comparisons when moving to another sample, but never compares a sample with itself", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    const overlay = (): Element | null => screen.getByTestId("workspace-compares").querySelector("[data-compare='11']");
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
    expect(screen.queryByRole("button", { name: "Remove A2" })).toBeNull();
    expect(overlay()).toBeNull();
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.getByRole("button", { name: "Remove A2" })).toBeInTheDocument();
    expect(overlay()).not.toBeNull();
  });

  it("drops a kept zoom that falls outside the next sample's q range, and Z resets the zoom", () => {
    state.trace = makeTrace(1);
    const { rerender } = renderAt("/experiments/1/workspace/10");
    const plot = screen.getByTestId("workspace-trace");
    const xTicks = (): string =>
      Array.from(plot.ownerDocument.querySelectorAll('[data-testid="workspace-trace"] [data-role="axis-bottom"] text')).map((t) => t.textContent).join(" ");
    const full = xTicks();
    for (let i = 0; i < 5; i++) fireEvent.wheel(plot.parentElement!, { deltaY: -300, clientX: 500 });
    expect(xTicks()).not.toBe(full);
    fireEvent.keyDown(window, { key: "z" });
    expect(xTicks()).toBe(full);
    for (let i = 0; i < 5; i++) fireEvent.wheel(plot.parentElement!, { deltaY: -300, clientX: 540 });
    const q = Array.from({ length: 21 }, (_, i) => 0.01 + i * 0.002);
    state.trace = { q, I: q.map(() => 5), sigma: q.map(() => 0) };
    rerender(page("/experiments/1/workspace/10"));
    expect(screen.getByTestId("workspace-trace").querySelector('[data-role="trace-line"] path')?.getAttribute("d")).toMatch(/^M/);
    expect(xTicks()).toContain("0.03");
  });
});
