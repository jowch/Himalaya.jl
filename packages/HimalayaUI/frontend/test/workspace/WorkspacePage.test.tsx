import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { act, render, screen, fireEvent, within } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import type { CorpusSample, Exposure, IndexEntry, Assignment, Peak, SampleBridge, Trace } from "../../src/api";
import { WorkspacePage } from "../../src/workspace/WorkspacePage";
import { peakList } from "../fixtures/peaks";

type Callbacks = { onSuccess?: () => void; onError?: () => void };
type Call = { hook: string; boundTo: number; arg: unknown; callbacks?: Callbacks };
const calls: Call[] = [];

const traceRefetch = vi.fn();
const state = {
  corpus: [] as CorpusSample[],
  exposures: [] as Exposure[],
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
  username: "jc" as string | undefined,
  trace: undefined as Trace | undefined,
  traceError: false,
  traceFetching: false,
  peaks: [] as Peak[],
  compareTraces: new Map<number, Trace>(),
  bridges: [] as SampleBridge[],
  loads: [] as { samples: { flag: object | null }[] }[],
  ingest: undefined as { status: string; stage?: string; processed: number; total: number } | undefined,
  rescanError: false,
};
const rescans: boolean[] = [];

function recorder(hook: string) {
  return (boundTo: number) => ({
    mutate: (arg: unknown, callbacks?: Callbacks) => {
      calls.push({ hook, boundTo, arg, ...(callbacks ? { callbacks } : {}) });
    },
  });
}

vi.mock("../../src/queries", () => ({
  useCorpusSamples: () => ({ data: state.corpus, isLoading: false }),
  useExperiment: () => ({ data: { id: 1, name: "BL-19", last_scanned_at: "2026-11-01T00:00:00" } }),
  useLoads: () => ({ data: state.loads }),
  useTriggerScan: () => ({ mutate: (force = false) => rescans.push(force), isError: state.rescanError }),
  useExposures: () => ({ data: state.exposures }),
  useTrace: () => ({ data: state.trace, isError: state.traceError, isFetching: state.traceFetching, refetch: traceRefetch }),
  usePeaks: () => ({ data: state.peaks }),
  useIndices: () => ({ data: state.indices }),
  useAssignment: () => ({ data: state.assignment }),
  useBridges: () => ({ data: state.bridges }),
  useBridgeSamples: recorder("bridge"),
  useUnbridgeSample: recorder("unbridge"),
  useAddAssignmentPhase: recorder("add"),
  useRemoveAssignmentPhase: recorder("remove"),
  useSetAssignmentState: recorder("state"),
  useSelectExposure: recorder("select"),
  useSetExposureStatus: recorder("status"),
  useAddPeak: recorder("addPeak"),
  useRemovePeak: recorder("removePeak"),
  useSetPeakExcluded: recorder("excludePeak"),
  useCommitCustomIndex: (boundTo: number) => ({
    mutate: (phase: string, basis: number, ratios: number[]) => calls.push({ hook: "custom", boundTo, arg: { phase, basis, ratios } }),
  }),
}));

vi.mock("../../src/workspace/useCompareTraces", () => ({
  useCompareTraces: (ids: readonly number[]) =>
    new Map(ids.flatMap((id) => (state.compareTraces.has(id) ? [[id, state.compareTraces.get(id)!] as const] : []))),
}));

function makeTrace(level: number): Trace {
  const q = Array.from({ length: 21 }, (_, i) => 0.05 + i * 0.01);
  return { q, I: q.map((x) => level * (1 + 10 * Math.exp(-((x - 0.1) ** 2) / 1e-5))), sigma: q.map(() => 0) };
}

function peak(id: number, q: number, source: "auto" | "manual", excluded = false): Peak {
  return { id, exposure_id: 101, q, intensity: null, prominence: null, sharpness: null, source, excluded };
}

function currentView(): string {
  fireEvent.keyDown(window, { key: "k", ctrlKey: true });
  const input = screen.getByLabelText("Command or sample");
  fireEvent.change(input, { target: { value: "current view" } });
  const option = screen.getAllByRole("option").find((o) => /current view/.test(o.textContent ?? ""));
  const name = /Plot (.+), current view/.exec(option?.textContent ?? "")?.[1] ?? "";
  fireEvent.keyDown(input, { key: "Escape" });
  return name;
}

function pickView(label: string): void {
  fireEvent.keyDown(window, { key: "k", ctrlKey: true });
  const input = screen.getByLabelText("Command or sample");
  fireEvent.change(input, { target: { value: label } });
  fireEvent.keyDown(input, { key: "Enter" });
}

function clickPlotAt(q: number, detail = 1): void {
  const plotWidth = 640 - 52 - 12;
  const px = 52 + ((q - 0.05) / 0.2) * plotWidth;
  fireEvent.click(screen.getByTestId("workspace-trace"), { clientX: px, clientY: 100, detail });
}

const clearUsername = vi.fn();
vi.mock("../../src/state", () => ({
  useAppState: (sel: (s: { username: string | undefined; clearUsername: () => void; ingestInFlight: Record<number, unknown> }) => unknown) =>
    sel({ username: state.username, clearUsername, ingestInFlight: state.ingest ? { 1: state.ingest } : {} }),
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
        <Route path="*" element={<Where />} />
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

afterEach(() => {
  window.localStorage.removeItem("himalaya-ui:workspace-theme");
  window.localStorage.removeItem("himalaya-ui:workspace-plot-view");
});

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
  state.traceError = false;
  state.traceFetching = false;
  state.peaks = [];
  state.compareTraces = new Map();
  state.bridges = [];
  state.loads = [];
  state.ingest = undefined;
  state.rescanError = false;
  rescans.length = 0;
});

describe("WorkspacePage", () => {
  it("lists only this experiment's samples and marks unreviewed ones", () => {
    renderAt("/experiments/1/workspace/10");
    const list = screen.getByTestId("workspace-samples");
    expect(within(list).queryByText("other")).toBeNull();
    expect(within(list).getAllByRole("img", { name: "Unreviewed" })).toHaveLength(2);
    expect(screen.getByTestId("workspace-progress")).toHaveTextContent("2 of 4 reviewed");
  });

  it("switches user from the palette, which names who is signed in", () => {
    clearUsername.mockClear();
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    const input = screen.getByLabelText("Command or sample");
    fireEvent.change(input, { target: { value: "switch user" } });
    expect(screen.getAllByRole("option").some((o) => /Switch user, signed in as jc/.test(o.textContent ?? ""))).toBe(true);
    fireEvent.keyDown(input, { key: "Enter" });
    expect(clearUsername).toHaveBeenCalledTimes(1);
  });

  it("says when the trace failed to load and retries on request", () => {
    state.traceError = true;
    traceRefetch.mockClear();
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByText("Couldn't load this sample")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "Retry" }));
    expect(traceRefetch).toHaveBeenCalledTimes(1);
  });

  it("disables Retry while the retry is in flight", () => {
    state.traceError = true;
    state.traceFetching = true;
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByRole("button", { name: "Retrying…" })).toBeDisabled();
  });

  it("holding Space opens the detector view once and does not close it", () => {
    state.exposures = [exposure(100), exposure(101, { selected: true, image_path: "/x.tif" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: " " });
    fireEvent.keyDown(screen.getByRole("dialog", { name: "Detector image" }), { key: " ", repeat: true });
    expect(screen.getByRole("dialog", { name: "Detector image" })).toBeInTheDocument();
  });

  it("Space enlarges the shown frame's detector image, and Space or Esc closes it", () => {
    state.exposures = [exposure(100), exposure(101, { selected: true, image_path: "/x.tif", image_version: "7" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: " " });
    const dialog = screen.getByRole("dialog", { name: "Detector image" });
    expect(dialog.querySelector("img:not([aria-hidden])")!.getAttribute("src")).toBe("/api/exposures/101/image?v=7");
    expect(screen.getByTestId("workspace-detector-caption")).toHaveTextContent("Frame 2 of 2 · f101.dat");
    fireEvent.keyDown(dialog, { key: " " });
    expect(screen.queryByRole("dialog", { name: "Detector image" })).toBeNull();
    fireEvent.keyDown(window, { key: " " });
    fireEvent.keyDown(screen.getByRole("dialog", { name: "Detector image" }), { key: "Escape" });
    expect(screen.queryByRole("dialog", { name: "Detector image" })).toBeNull();
  });

  it("flips through frames inside the enlarged view and drops one without closing it", () => {
    state.exposures = [
      exposure(100, { image_path: "/a.tif" }),
      exposure(101, { selected: true, image_path: "/b.tif" }),
      exposure(102, { image_path: "/c.tif" }),
    ];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: " " });
    const dialog = screen.getByRole("dialog", { name: "Detector image" });
    expect(dialog.querySelectorAll("img")).toHaveLength(3);
    const caption = screen.getByTestId("workspace-detector-caption");
    expect(caption).toHaveTextContent("Frame 2 of 3");
    fireEvent.keyDown(dialog, { key: "ArrowRight" });
    expect(caption).toHaveTextContent("Frame 3 of 3 · f102.dat");
    expect(dialog.querySelector("img:not([aria-hidden])")!.getAttribute("data-exposure-id")).toBe("102");
    fireEvent.keyDown(dialog, { key: "ArrowLeft" });
    fireEvent.keyDown(dialog, { key: "ArrowLeft" });
    expect(caption).toHaveTextContent("Frame 1 of 3");
    fireEvent.keyDown(dialog, { key: "x" });
    expect(calls.filter((c) => c.hook === "status").map((c) => c.arg)).toEqual([{ exposureId: 100, status: "rejected" }]);
    expect(screen.getByRole("dialog", { name: "Detector image" })).toBeInTheDocument();
    expect(currentFrame()).toBe("Frame 1");
  });

  it("does not drop from the enlarged view before a user is set", () => {
    state.username = undefined;
    state.exposures = [exposure(100, { selected: true, image_path: "/a.tif" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: " " });
    fireEvent.keyDown(screen.getByRole("dialog", { name: "Detector image" }), { key: "x" });
    expect(calls.filter((c) => c.hook === "status")).toHaveLength(0);
  });

  it("Space does nothing on a frame without a detector image", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: " " });
    expect(screen.queryByRole("dialog", { name: "Detector image" })).toBeNull();
  });

  it("links the Himalaya brand to all experiments", () => {
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByRole("link", { name: "Himalaya" })).toHaveAttribute("href", "/experiments");
  });

  it("opens an experiment without a sample on its first unreviewed sample", () => {
    renderAt("/experiments/1/workspace");
    expect(screen.getByRole("heading", { level: 2, name: "A2" })).toBeInTheDocument();
  });

  it("puts the landing sample in the URL, so calling its phase does not move on", () => {
    const { rerender } = renderAt("/experiments/1/workspace");
    expect(screen.getByTestId("where")).toHaveTextContent("/experiments/1/workspace/11");
    state.corpus = state.corpus.map((s) => (s.id === 11 ? { ...s, assignment_state: "indexed", phase: "Im3m" } : s));
    rerender(page("/experiments/1/workspace"));
    expect(screen.getByRole("heading", { level: 2, name: "A2" })).toBeInTheDocument();
  });

  it("says when the last rescan failed, and Details opens the failure", () => {
    state.ingest = { status: "failed", processed: 0, total: 0 };
    renderAt("/experiments/1/workspace/10");
    const failed = screen.getByTestId("workspace-rescan-failed");
    expect(failed).toHaveTextContent("The last rescan failed.");
    fireEvent.click(within(failed).getByRole("button", { name: "Details" }));
    expect(screen.getByTestId("where")).toHaveTextContent(/^\/experiments\/1$/);
  });

  it("says the rescan failed when the request to start it fails", () => {
    state.rescanError = true;
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-rescan-failed")).toBeInTheDocument();
  });

  it("shows nothing above the list while grouping is settled and no rescan runs", () => {
    state.loads = [{ samples: [{ flag: null }] }];
    renderAt("/experiments/1/workspace/10");
    expect(screen.queryByTestId("workspace-notice")).toBeNull();
  });

  it("says how many samples need a grouping check, and Review opens the grouping review", () => {
    state.loads = [{ samples: [{ flag: { kind: "merge" } }, { flag: null }, { flag: { kind: "split" } }] }];
    renderAt("/experiments/1/workspace/10");
    const notice = screen.getByTestId("workspace-notice");
    expect(notice).toHaveTextContent("2 samples need a grouping check.");
    fireEvent.click(within(notice).getByRole("button", { name: "Review" }));
    expect(screen.getByTestId("where")).toHaveTextContent("/experiments/1/grouping");
  });

  it("shows rescan progress and hides the Rescan command while it runs", () => {
    state.ingest = { status: "analyzing", processed: 5, total: 10 };
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-rescan")).toHaveTextContent("Rescanning: 5 / 10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "rescan" } });
    expect(screen.queryAllByRole("option")).toHaveLength(0);
  });

  it("rescans, configures and reviews grouping from the palette", () => {
    renderAt("/experiments/1/workspace/10");
    const run = (label: string): void => {
      fireEvent.keyDown(window, { key: "k", ctrlKey: true });
      const input = screen.getByLabelText("Command or sample");
      fireEvent.change(input, { target: { value: label } });
      fireEvent.keyDown(input, { key: "Enter" });
    };
    run("rescan");
    expect(rescans).toEqual([false]);
    run("configure experiment");
    expect(screen.getByTestId("where")).toHaveTextContent("/experiments/1/config");
  });

  it("updates the list and the count when the corpus status changes", () => {
    const { rerender } = renderAt("/experiments/1/workspace/11");
    state.corpus = state.corpus.map((s) => (s.id === 11 ? { ...s, assignment_state: "indexed", phase: "Im3m" } : s));
    rerender(page("/experiments/1/workspace/11"));
    const row = within(screen.getByTestId("workspace-samples")).getByRole("button", { current: true });
    expect(within(row).getByTestId("workspace-status")).toHaveTextContent("Im3m");
    expect(screen.getByTestId("workspace-progress")).toHaveTextContent("3 of 4 reviewed");
  });

  it("Z jumps to the next unreviewed sample, like Gradescope's next ungraded", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "z" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
  });

  it("J, L and R are no longer shortcuts", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(screen.getByRole("dialog").querySelector("[data-sample-id='11']")!);
    expect(screen.getByRole("button", { name: "A2" })).toBeInTheDocument();
    fireEvent.keyDown(window, { key: "j" });
    fireEvent.keyDown(window, { key: "l" });
    fireEvent.keyDown(window, { key: "r" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/10");
    expect(currentView()).toBe("Log");
    expect(screen.getByRole("button", { name: "A2" })).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: "A2, hidden" })).toBeNull();
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

  it("every exposure card has its own Drop button, which drops that frame and shows it", () => {
    renderAt("/experiments/1/workspace/10");
    const cards = within(screen.getByTestId("workspace-exposures")).getAllByRole("listitem");
    fireEvent.click(within(cards[0]!).getByRole("button", { name: "Drop" }));
    expect(calls[0]).toEqual({ hook: "status", boundTo: 10, arg: { exposureId: 100, status: "rejected" } });
    expect(currentFrame()).toBe("Frame 1");
  });

  it("a dropped card offers Restore and the header counts the kept frames", () => {
    state.exposures = [exposure(100, { status: "rejected" }), exposure(101, { selected: true })];
    renderAt("/experiments/1/workspace/10");
    const strip = screen.getByTestId("workspace-exposures");
    expect(within(within(strip).getAllByRole("listitem")[0]!).getByRole("button", { name: "Restore" })).toBeInTheDocument();
    expect(screen.getByText("1 of 2 kept")).toBeInTheDocument();
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
    fireEvent.keyDown(button, { key: "z" });
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
    state.peaks = peakList(101, { auto: [0.1, 0.12], manual: [0.2] });
    renderAt("/experiments/1/workspace/10");
    clickPlotAt(0.1);
    clickPlotAt(0.2);
    clickPlotAt(0.15);
    expect(calls.slice(0, 2)).toEqual([
      { hook: "excludePeak", boundTo: 101, arg: { peakId: 1, excluded: true } },
      { hook: "removePeak", boundTo: 101, arg: 1 },
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

  it("a click on a disabled auto peak restores it, even when a manual peak shares its id", () => {
    state.trace = makeTrace(1);
    state.peaks = peakList(101, { auto: [{ q: 0.2, excluded: true }], manual: [0.1] });
    renderAt("/experiments/1/workspace/10");
    clickPlotAt(0.2);
    expect(calls).toEqual([{ hook: "excludePeak", boundTo: 101, arg: { peakId: 1, excluded: false } }]);
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

  it("the palette steps through log, linear, log-log and q²I", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    const ticks = (): string =>
      Array.from(screen.getByTestId("workspace-trace").querySelectorAll('[data-role="axis-left"] text'))
        .map((t) => t.textContent)
        .join(" ");
    const logTicks = ticks();
    expect(currentView()).toBe("Log");
    pickView("Plot Linear");
    expect(currentView()).toBe("Linear");
    expect(ticks()).not.toBe(logTicks);
    pickView("Plot Log-log");
    expect(currentView()).toBe("Log-log");
    pickView("Plot q²I");
    expect(currentView()).toBe("q²I (Kratky)");
    expect(screen.getByTestId("workspace-trace").textContent).toContain("q²·I");
    pickView("Plot Log");
    expect(currentView()).toBe("Log");
    expect(ticks()).toBe(logTicks);
  });

  it("the difference strip follows the plot view", () => {
    state.trace = makeTrace(1);
    const flat = makeTrace(1);
    state.compareTraces = new Map([[11, { ...flat, I: flat.q.map(() => 1) }]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    const strip = (): string => screen.getByTestId("workspace-subtraction").querySelector("[data-role='difference-line']")!.getAttribute("d")!;
    const log = strip();
    pickView("Plot Log-log");
    const loglog = strip();
    expect(loglog).not.toBe(log);
    pickView("Plot q²I");
    expect(strip()).not.toBe(log);
    expect(strip()).not.toBe(loglog);
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

  it("scales by the transmitted beam when both traces carry monitors, even where high-q matching would disagree", () => {
    const frame = (exposure_id: number, transmitted_counts: number) => ({ exposure_id, transmitted_counts, exposure_time: 15, counting_time: 15.8 });
    state.trace = { ...makeTrace(1), frames: [frame(101, 6000)] };
    state.compareTraces = new Map([[11, { ...makeTrace(2), frames: [frame(111, 8000)] }]]);
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

  it("clicking the subtracted chip stops subtracting, the palette hides comparisons, and × removes one", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    fireEvent.click(screen.getByRole("button", { name: "A2" }));
    expect(screen.queryByTestId("workspace-subtraction")).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "A2" }));
    expect(screen.getByTestId("workspace-subtraction")).toBeInTheDocument();
    pickView("Hide comparisons");
    expect(screen.queryByTestId("workspace-subtraction")).toBeNull();
    expect(screen.getByTestId("workspace-compares").children).toHaveLength(0);
    expect(screen.getByRole("button", { name: "A2, hidden" })).toBeInTheDocument();
    pickView("Show comparisons");
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
    fireEvent.keyDown(window, { key: "0", ctrlKey: true });
    expect(xTicks()).toBe(full);
    for (let i = 0; i < 5; i++) fireEvent.wheel(plot.parentElement!, { deltaY: -300, clientX: 540 });
    const q = Array.from({ length: 21 }, (_, i) => 0.01 + i * 0.002);
    state.trace = { q, I: q.map(() => 5), sigma: q.map(() => 0) };
    rerender(page("/experiments/1/workspace/10"));
    expect(screen.getByTestId("workspace-trace").querySelector('[data-role="trace-line"] path')?.getAttribute("d")).toMatch(/^M/);
    expect(xTicks()).toContain("0.03");
  });
  it("? opens the shortcut list, the header button opens it too, and Escape closes it", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "?", shiftKey: true });
    const help = screen.getByRole("dialog", { name: "Keyboard shortcuts" });
    expect(within(help).getByText("Reset the zoom")).toBeInTheDocument();
    fireEvent.keyDown(help, { key: "Escape" });
    expect(screen.queryByRole("dialog")).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "Keyboard shortcuts" }));
    expect(screen.getByRole("dialog", { name: "Keyboard shortcuts" })).toBeInTheDocument();
  });

  it("shortcuts do nothing while the shortcut list is open", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "?", shiftKey: true });
    const help = screen.getByRole("dialog", { name: "Keyboard shortcuts" });
    expect(document.activeElement).toBe(help);
    fireEvent.keyDown(help, { key: "x" });
    fireEvent.keyDown(help, { key: "z" });
    fireEvent.keyDown(help, { key: "ArrowDown" });
    expect(calls).toEqual([]);
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/10");
  });

  it("⌘K over the compare picker replaces it, so only one dialog is open", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    expect(screen.getByRole("dialog", { name: "Compare with a sample" })).toBeInTheDocument();
    fireEvent.keyDown(screen.getByLabelText("Find a sample to compare"), { key: "k", metaKey: true });
    expect(screen.getAllByRole("dialog").map((d) => d.getAttribute("aria-label"))).toEqual(["Commands"]);
  });

  it("⌘ Shift K is left alone", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "K", metaKey: true, shiftKey: true });
    expect(screen.queryByRole("dialog")).toBeNull();
  });

  it("the palette shows each command's key as a keycap and tracks the cursor for screen readers", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", metaKey: true });
    const input = screen.getByRole("combobox", { name: "Command or sample" });
    fireEvent.change(input, { target: { value: "reset zoom" } });
    const option = screen.getByRole("option", { name: /Reset zoom/ });
    expect(option.querySelector("kbd")?.textContent).toBe("Ctrl 0");
    expect(input.getAttribute("aria-activedescendant")).toBe(option.id);
    fireEvent.change(input, { target: { value: "go" } });
    fireEvent.keyDown(input, { key: "ArrowDown" });
    const second = screen.getAllByRole("option")[1]!;
    expect(second).toHaveAttribute("aria-selected", "true");
    expect(input.getAttribute("aria-activedescendant")).toBe(second.id);
    expect(second.querySelector("kbd")).toBeNull();
  });

  it("⌘K opens the command palette and Enter goes to the matching sample", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", metaKey: true });
    const input = screen.getByLabelText("Command or sample");
    expect(document.activeElement).toBe(input);
    fireEvent.change(input, { target: { value: "go a3" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/12");
  });

  it("Ctrl+K toggles the palette, and a command runs the same action as its key", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "k", ctrlKey: true });
    expect(screen.queryByRole("dialog")).toBeNull();
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "kratky" } });
    fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Enter" });
    expect(currentView()).toBe("q²I (Kratky)");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "plot" } });
    expect(screen.getByRole("option", { name: /q²I \(Kratky\), current view/ })).toBeInTheDocument();
    fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Escape" });
  });

  it("the palette switches the workspace to dark and back, remembers it, and can follow the system again", () => {
    const { unmount } = renderAt("/experiments/1/workspace/10");
    const root = screen.getByTestId("workspace");
    expect(root).toHaveAttribute("data-workspace-theme", "system");
    fireEvent.keyDown(window, { key: "k", metaKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "dark theme" } });
    fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Enter" });
    expect(root).toHaveAttribute("data-workspace-theme", "dark");
    unmount();
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace")).toHaveAttribute("data-workspace-theme", "dark");
    fireEvent.keyDown(window, { key: "k", metaKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "light theme" } });
    fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Enter" });
    expect(screen.getByTestId("workspace")).toHaveAttribute("data-workspace-theme", "light");
    fireEvent.keyDown(window, { key: "k", metaKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "system theme" } });
    fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Enter" });
    expect(screen.getByTestId("workspace")).toHaveAttribute("data-workspace-theme", "system");
    expect(window.localStorage.getItem("himalaya-ui:workspace-theme")).toBeNull();
    fireEvent.keyDown(window, { key: "k", metaKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "system theme" } });
    expect(screen.queryByRole("option")).toBeNull();
  });

  it("I opens a custom index on the top candidate; arrows nudge the lattice and Enter adds it", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    const panel = screen.getByTestId("workspace-custom-index");
    expect(within(panel).getByRole("button", { name: "Pn3m" })).toHaveAttribute("aria-pressed", "true");
    fireEvent.keyDown(window, { key: "ArrowRight" });
    fireEvent.keyDown(window, { key: "ArrowRight", shiftKey: true });
    expect(within(panel).getByText(/Add Pn3m 105\.5 Å/)).toBeInTheDocument();
    expect(currentFrame()).toBe("Frame 2");
    fireEvent.keyDown(window, { key: "Enter" });
    expect(screen.queryByTestId("workspace-custom-index")).toBeNull();
    expect(calls).toHaveLength(1);
    const c = calls[0]!;
    expect(c).toMatchObject({ hook: "custom", boundTo: 101 });
    const arg = c.arg as { phase: string; basis: number; ratios: number[] };
    expect(arg.phase).toBe("Pn3m");
    expect(arg.basis).toBeCloseTo((2 * Math.PI * Math.SQRT2) / 105.5, 6);
    expect(arg.ratios[0]).toBe(1);
    expect(arg.ratios[1]).toBeCloseTo(Math.sqrt(3 / 2), 6);
  });

  it("Escape cancels a custom index without writing", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.queryByTestId("workspace-custom-index")).toBeNull();
    expect(calls).toEqual([]);
  });

  it("a custom index on an unselected frame selects it first", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    fireEvent.keyDown(window, { key: "i" });
    fireEvent.click(within(screen.getByTestId("workspace-custom-index")).getByRole("button", { name: /^Add/ }));
    expect(calls[0]).toMatchObject({ hook: "select", boundTo: 10, arg: 100 });
    calls[0]!.callbacks!.onSuccess!();
    expect(calls[1]).toMatchObject({ hook: "custom", boundTo: 100 });
  });

  it("does not open a custom index on a dropped frame or before a user is set", () => {
    state.username = undefined;
    const { unmount } = renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    expect(screen.queryByTestId("workspace-custom-index")).toBeNull();
    unmount();
    state.username = "jc";
    state.exposures = [exposure(100), exposure(101, { selected: true, status: "rejected" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(currentFrame()).toBe("Frame 2, dropped");
    fireEvent.keyDown(window, { key: "i" });
    expect(screen.queryByTestId("workspace-custom-index")).toBeNull();
  });

  it("the indexing plot shows the assigned candidate's fit against the peaks", () => {
    state.trace = makeTrace(1);
    state.peaks = [peak(7, 0.1, "auto"), peak(8, 0.14, "auto"), peak(9, 0.2, "auto", true)];
    state.indices = [
      index(1, "Im3m", 0.4),
      {
        ...index(2, "Pn3m", 0.9),
        peaks: [
          { peak_id: 7, ratio_position: 1, residual: 0, q_observed: 0.1 },
          { peak_id: 8, ratio_position: 2, residual: 0, q_observed: 0.14 },
        ],
      },
    ];
    renderAt("/experiments/1/workspace/10");
    const plot = screen.getByTestId("workspace-indexing");
    expect(within(plot).getByTestId("workspace-indexing-fit")).toHaveTextContent("Pn3m a = 100.0 Å R² 0.9900 2/2 on peaks");
    expect(plot.querySelectorAll("[data-observed]")).toHaveLength(2);
  });

  it("draws empty indexing axes before any candidate is hovered or ticked", () => {
    state.assignment = { exposure_id: 101, state: "indexed", members: [] };
    renderAt("/experiments/1/workspace/10");
    const plot = screen.getByTestId("workspace-indexing");
    expect(plot.querySelector("svg")).not.toBeNull();
    expect(within(plot).queryByTestId("workspace-indexing-fit")).toBeNull();
  });

  it("titles the phase list Call and the exposure details Exposure", () => {
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByRole("heading", { name: "Call" })).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: "Exposure" })).toBeInTheDocument();
  });

  it("lists the shown exposure's time, I₀ and I signal under the .prp names", () => {
    state.exposures = [exposure(100), exposure(101, { selected: true, exposure_time: 10, incident_monitor: 1.5007, transmitted_counts: 7798.4 })];
    renderAt("/experiments/1/workspace/10");
    const section = screen.getByRole("heading", { name: "Exposure" }).closest("section")!;
    const text = section.textContent ?? "";
    expect(text).toContain("Exposure10.0 s");
    expect(text).toContain("I₀1.501");
    expect(text).toContain("I signal7,798");
  });

  it("the indexing plot follows the custom index while it is open", () => {
    state.trace = makeTrace(1);
    const a = (2 * Math.PI * Math.SQRT2) / 0.1;
    state.peaks = [peak(7, 0.1, "auto"), peak(8, 0.1 * Math.sqrt(1.5), "auto")];
    state.indices = [{ ...index(2, "Pn3m", 0.9), lattice_d: a }];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    expect(screen.getByTestId("workspace-custom-landed")).toHaveTextContent("2 of 6");
    expect(screen.getByTestId("workspace-indexing-fit")).toHaveTextContent("2/6 on peaks");
    const combIds = Array.from(screen.getByTestId("workspace-combs").querySelectorAll("[data-comb]")).map((l) => l.getAttribute("data-comb"));
    expect(new Set(combIds)).toEqual(new Set(["-1"]));
  });

  it("an added custom index is listed as Custom and ticked, though it has no score", () => {
    state.indices = [index(1, "Im3m", 0.4), { ...index(5, "Hexagonal", 0), score: null, kind: "speculative", r_squared: null }];
    state.assignment = { exposure_id: 101, state: "indexed", members: [5] };
    renderAt("/experiments/1/workspace/10");
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows[0]).toHaveTextContent("Custom");
    expect(rows[0]).toHaveAttribute("aria-pressed", "true");
    expect(rows[1]).toHaveTextContent("0.40");
  });

  it("Enter on a button in the custom panel presses that button and does not add", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    const panel = screen.getByTestId("workspace-custom-index");
    const im3m = within(panel).getByRole("button", { name: "Im3m" });
    fireEvent.keyDown(im3m, { key: "Enter" });
    fireEvent.keyDown(within(panel).getByRole("button", { name: /^Cancel/ }), { key: "Enter" });
    expect(calls).toEqual([]);
    expect(screen.getByTestId("workspace-custom-index")).toBeInTheDocument();
  });

  it("arrows nudge the lattice after a symmetry click, while focus is on that button", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    const panel = screen.getByTestId("workspace-custom-index");
    const im3m = within(panel).getByRole("button", { name: "Im3m" });
    fireEvent.click(im3m);
    fireEvent.keyDown(im3m, { key: "ArrowRight" });
    expect(within(panel).getByRole("button", { name: /^Add/ })).toHaveTextContent("Add Im3m 100.5 Å");
  });

  it("a lattice typed outside the range is flagged and blocks Add until it is clamped on blur", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    const panel = screen.getByTestId("workspace-custom-index");
    const field = within(panel).getByLabelText("Lattice a value");
    const add = within(panel).getByRole("button", { name: /^Add/ });
    fireEvent.focus(field);
    fireEvent.change(field, { target: { value: "5" } });
    expect(field).toHaveAttribute("aria-invalid", "true");
    expect(add).toBeDisabled();
    expect(within(panel).getByTestId("workspace-custom-invalid")).toHaveTextContent("Enter a value from 10.0 to");
    fireEvent.keyDown(field, { key: "Enter" });
    expect(calls).toEqual([]);
    expect(field).toHaveValue(10);
    expect(field).toHaveAttribute("aria-invalid", "false");
    expect(add).toHaveTextContent("Add Pn3m 10.0 Å");
    expect(within(panel).getByTestId("workspace-custom-landed")).toHaveTextContent("No reflections fall inside this trace's q range.");
    fireEvent.focus(field);
    fireEvent.change(field, { target: { value: "" } });
    expect(add).toBeDisabled();
    fireEvent.blur(field);
    expect(field).toHaveValue(10);
    fireEvent.focus(field);
    fireEvent.change(field, { target: { value: "123.4" } });
    expect(add).toHaveTextContent("Add Pn3m 123.4 Å");
    fireEvent.blur(field);
    fireEvent.click(add);
    expect((calls[0]!.arg as { basis: number }).basis).toBeCloseTo((2 * Math.PI * Math.SQRT2) / 123.4, 6);
  });

  it("while a custom index is open, the single-key shortcuts do nothing", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    for (const key of ["x", "1", "2", "c", "j", "l", "r", "z", "ArrowDown", "ArrowUp"]) fireEvent.keyDown(window, { key });
    expect(calls).toEqual([]);
    expect(screen.getByTestId("where")).toHaveTextContent("/experiments/1/workspace/10");
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(screen.getByTestId("workspace-custom-index")).toBeInTheDocument();
    fireEvent.keyDown(window, { key: "?" });
    expect(screen.getByRole("dialog")).toBeInTheDocument();
  });

  it("while a custom index is open, a plot click puts the first order on the nearest peak and edits no peaks", () => {
    state.trace = makeTrace(1);
    state.peaks = peakList(101, { auto: [0.2], manual: [0.1] });
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    clickPlotAt(0.2);
    const add = within(screen.getByTestId("workspace-custom-index")).getByRole("button", { name: /^Add/ });
    expect(add).toHaveTextContent(`Add Pn3m ${((2 * Math.PI * Math.SQRT2) / 0.2).toFixed(1)} Å`);
    clickPlotAt(0.1);
    expect(add).toHaveTextContent(`Add Pn3m ${((2 * Math.PI * Math.SQRT2) / 0.1).toFixed(1)} Å`);
    expect(calls).toEqual([]);
  });

  it("moving to another frame closes the custom index", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    fireEvent.click(screen.getByRole("button", { name: "Frame 1" }));
    expect(currentFrame()).toBe("Frame 1");
    expect(screen.queryByTestId("workspace-custom-index")).toBeNull();
  });

  it("says so and keeps the custom index open when the frame cannot be selected", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    fireEvent.keyDown(window, { key: "i" });
    fireEvent.click(within(screen.getByTestId("workspace-custom-index")).getByRole("button", { name: /^Add/ }));
    act(() => calls[0]!.callbacks!.onError!());
    expect(screen.getByRole("alert")).toHaveTextContent("could not be selected");
    expect(screen.getByTestId("workspace-custom-index")).toBeInTheDocument();
    expect(calls.filter((c) => c.hook === "custom")).toEqual([]);
  });

  it("hovering a candidate shows its fit in place of the ticked ones", () => {
    state.trace = makeTrace(1);
    state.peaks = [peak(7, 0.1, "auto")];
    renderAt("/experiments/1/workspace/10");
    expect(screen.getAllByTestId("workspace-indexing-fit").map((n) => n.textContent)).toEqual([expect.stringMatching(/^Pn3m/)]);
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    fireEvent.mouseEnter(rows[1]!);
    expect(screen.getAllByTestId("workspace-indexing-fit").map((n) => n.textContent)).toEqual([expect.stringMatching(/^Im3m/)]);
    fireEvent.mouseLeave(screen.getByTestId("workspace-candidates"));
    expect(screen.getAllByTestId("workspace-indexing-fit").map((n) => n.textContent)).toEqual([expect.stringMatching(/^Pn3m/)]);
  });
  it("the plot draws a candidate's reflections only while it is hovered, not because it is ticked", () => {
    state.trace = makeTrace(1);
    state.peaks = [peak(7, 0.1, "auto")];
    renderAt("/experiments/1/workspace/10");
    const combIds = (): string[] => Array.from(screen.getByTestId("workspace-combs").querySelectorAll("[data-comb]")).map((l) => l.getAttribute("data-comb") ?? "");
    expect(combIds()).toEqual([]);
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    fireEvent.mouseEnter(rows[1]!);
    expect(combIds().length).toBeGreaterThan(0);
    expect(new Set(combIds())).toEqual(new Set(["1"]));
    fireEvent.mouseLeave(screen.getByTestId("workspace-candidates"));
    expect(combIds()).toEqual([]);
  });
  it("the indexing plot shows the peaks the stored index claimed, even where a fresh match would differ", () => {
    state.trace = makeTrace(1);
    state.peaks = [peak(7, 0.1, "auto"), peak(8, 0.14, "auto")];
    state.indices = [{ ...index(2, "Pn3m", 0.9), peaks: [{ peak_id: 7, ratio_position: 1, residual: 0, q_observed: 0.1 }] }];
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-indexing-fit")).toHaveTextContent("1/2 on peaks");
  });
  it("Enter on a focused button outside the custom panel does not add the index", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    const row = screen.getByRole("button", { name: /A2/ });
    fireEvent.keyDown(row, { key: "Enter" });
    expect(calls).toEqual([]);
    expect(screen.getByTestId("workspace-custom-index")).toBeInTheDocument();
  });

  it("the lattice field keeps what is typed, keystroke by keystroke", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    const field = within(screen.getByTestId("workspace-custom-index")).getByLabelText("Lattice a value") as HTMLInputElement;
    fireEvent.focus(field);
    for (const typed of ["2", "25", "250"]) {
      fireEvent.change(field, { target: { value: typed } });
      expect(field.value).toBe(typed);
    }
    fireEvent.blur(field);
    expect(field.value).toBe("250.0");
  });
  it("opens on the frame a link names, and falls back to the representative for a frame of another sample", () => {
    const { unmount } = renderAt("/experiments/1/workspace/10?exposure=100");
    expect(currentFrame()).toBe("Frame 1");
    unmount();
    renderAt("/experiments/1/workspace/10?exposure=999");
    expect(currentFrame()).toBe("Frame 2");
  });

  it("0 marks the frame as no phase (form factor), and 0 again clears it", () => {
    state.indices = [];
    renderAt("/experiments/1/workspace/10");
    const row = screen.getByTestId("workspace-form-factor");
    expect(row).toHaveAttribute("aria-pressed", "false");
    fireEvent.keyDown(window, { key: "0" });
    expect(calls).toEqual([{ hook: "state", boundTo: 101, arg: "form_factor" }]);
    calls.length = 0;
    state.assignment = { exposure_id: 101, state: "form_factor", members: [] };
    fireEvent.click(screen.getByRole("button", { name: "Frame 1" }));
    fireEvent.click(screen.getByRole("button", { name: "Frame 2" }));
    expect(screen.getByTestId("workspace-form-factor")).toHaveAttribute("aria-pressed", "true");
    fireEvent.click(screen.getByTestId("workspace-form-factor"));
    expect(calls).toEqual([{ hook: "state", boundTo: 101, arg: "null" }]);
  });

  it("on an unselected frame, marks form factor only after the select succeeds", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    fireEvent.keyDown(window, { key: "0" });
    expect(calls).toHaveLength(1);
    expect(calls[0]).toMatchObject({ hook: "select", boundTo: 10, arg: 100 });
    calls[0]!.callbacks!.onSuccess!();
    expect(calls[1]).toEqual({ hook: "state", boundTo: 100, arg: "form_factor" });
  });

  it("form factor is not offered on a dropped frame or before a user is set", () => {
    state.username = undefined;
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-form-factor")).toBeDisabled();
    fireEvent.keyDown(window, { key: "0" });
    expect(calls).toEqual([]);
  });
  it("draws bridged samples on the plot as soon as the sample opens, and marks them in the list", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[12, makeTrace(3)]]);
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-compares").querySelector("[data-compare='12']")?.getAttribute("d")).toMatch(/^M/);
    expect(screen.getByRole("button", { name: "Remove bridge to A3" })).toBeInTheDocument();
    const list = screen.getByTestId("workspace-samples");
    expect(list.querySelector("[data-sample-id='10'] [data-bridged]")).toHaveTextContent("1 bridge");
    expect(list.querySelector("[data-sample-id='12'] [data-bridged]")).not.toBeNull();
    expect(list.querySelector("[data-sample-id='11'] [data-bridged]")).toBeNull();
    expect(screen.queryByText(/^minus /)).toBeNull();
  });

  it("the palette opens a bridge picker that leaves out this sample and its bridges, and Enter bridges to the pick", () => {
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    pickView("Bridge A1 to a sample");
    const picker = screen.getByRole("dialog", { name: "Bridge A1 to a sample" });
    expect(within(picker).queryByText("A1")).toBeNull();
    expect(within(picker).queryByText("A3")).toBeNull();
    const input = within(picker).getByLabelText("Find a sample to bridge to");
    fireEvent.change(input, { target: { value: "a4" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(calls).toContainEqual({ hook: "bridge", boundTo: 1, arg: { sampleId: 10, toSampleId: 13 } });
  });

  it("the palette takes the open sample out of its bridges, and × on a bridged chip takes that sample out", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[12, makeTrace(3)]]);
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    pickView("Take A1 out of its bridges");
    expect(calls).toContainEqual({ hook: "unbridge", boundTo: 1, arg: { sampleId: 10 } });
    fireEvent.click(screen.getByRole("button", { name: "Remove bridge to A3" }));
    expect(calls).toContainEqual({ hook: "unbridge", boundTo: 1, arg: { sampleId: 12 } });
  });

  it("a bridged sample can be the subtraction reference", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[12, makeTrace(3)]]);
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    fireEvent.click(screen.getByRole("button", { name: "A3" }));
    expect(screen.getByTestId("workspace-subtraction")).toBeInTheDocument();
    expect(screen.getByText("minus A3")).toBeInTheDocument();
  });

  it("bridging while comparisons are hidden shows them again", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[12, makeTrace(3)]]);
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    pickView("Hide comparisons");
    expect(screen.getByTestId("workspace-compares").children).toHaveLength(0);
    pickView("Bridge A1 to a sample");
    const input = within(screen.getByRole("dialog")).getByLabelText("Find a sample to bridge to");
    fireEvent.change(input, { target: { value: "a4" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(screen.getByTestId("workspace-compares").children).toHaveLength(1);
  });

  it("B is not a shortcut, the palette offers no unbridging without bridges, and no bridging before a user is set", () => {
    const options = (): string[] => {
      fireEvent.keyDown(window, { key: "k", ctrlKey: true });
      fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "bridge" } });
      const labels = screen.queryAllByRole("option").map((o) => o.textContent ?? "");
      fireEvent.keyDown(screen.getByRole("dialog"), { key: "Escape" });
      return labels;
    };
    const first = renderAt("/experiments/1/workspace/11");
    fireEvent.keyDown(window, { key: "b" });
    fireEvent.keyDown(window, { key: "B", shiftKey: true });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(calls.filter((c) => c.hook === "unbridge")).toHaveLength(0);
    expect(options()).toEqual(["Bridge A2 to a sample"]);
    first.unmount();
    state.username = undefined;
    renderAt("/experiments/1/workspace/11");
    expect(options()).toEqual([]);
  });
});
