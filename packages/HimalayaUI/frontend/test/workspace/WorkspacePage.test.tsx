import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { act, render, screen, fireEvent, within, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route, useLocation } from "react-router-dom";
import type { CorpusSample, Exposure, IndexEntry, Assignment, Peak, ProblemExposure, SampleBridge, SampleCall, SampleCallMember, SampleLabel, LabelColorEntry, LabelColor, Trace } from "../../src/api";
import { ApiError } from "../../src/api";
import { WorkspacePage } from "../../src/workspace/WorkspacePage";
import type { PrefetchTarget } from "../../src/workspace/prefetch";
import { peakList } from "../fixtures/peaks";
import { makeExposure } from "../helpers/exposure";
import { useSyncExternalStore } from "react";
import type { CompareChange, CompareState } from "../../src/api";
import { applyCompareChange, EMPTY_COMPARES } from "../../src/lib/compares";
import { COMPARES_FULL, MAX_COMPARES } from "../../src/workspace/compare";
import { COMPARE_COLORS } from "../../src/ui/WorkspacePlot";
import { phaseColor } from "../../src/phases";

type Callbacks = { onSuccess?: () => void; onError?: () => void };
type Call = { hook: string; boundTo: number; arg: unknown; callbacks?: Callbacks };
const calls: Call[] = [];

const traceRefetch = vi.fn();
const corpusRequests: number[] = [];
const state = {
  experimentLoaded: true,
  experiments: [] as { id: number; name: string }[],
  allCorpus: [] as CorpusSample[],
  allCorpusEnabled: [] as boolean[],
  setups: [] as ({ id: number; name: string; folder: string; included: boolean } & Record<string, unknown>)[],
  detectorSize: undefined as { width: number; height: number } | undefined,
  sizeRequests: [] as string[],
  corpus: [] as CorpusSample[],
  exposures: [] as Exposure[],
  mateExposures: new Map<number, Exposure[]>(),
  indices: [] as IndexEntry[],
  assignment: undefined as Assignment | undefined,
  sampleCall: undefined as SampleCall | undefined,
  frameAssignments: new Map<number, Assignment>(),
  username: "jc" as string | undefined,
  trace: undefined as Trace | undefined,
  traceError: false,
  traceErr: null as Error | null,
  traceFetching: false,
  peaks: [] as Peak[],
  compareTraces: new Map<number, Trace>(),
  bridges: [] as SampleBridge[],
  labels: [] as SampleLabel[],
  problems: { ingest_error: null as string | null, exposures: [] as ProblemExposure[] },
  labelColors: [] as LabelColorEntry[],
  ingest: undefined as { status: string; stage?: string; processed: number; total: number } | undefined,
  rescanError: false,
};
const compareStore = { data: EMPTY_COMPARES as CompareState, listeners: new Set<() => void>() };
function setCompares(next: CompareState): void {
  compareStore.data = next;
  compareStore.listeners.forEach((l) => l());
}
const subscribeCompares = (l: () => void): (() => void) => {
  compareStore.listeners.add(l);
  return () => compareStore.listeners.delete(l);
};
const rescans: boolean[] = [];
const removals: number[] = [];
const opened: number[] = [];

function recorder(hook: string) {
  return (boundTo: number) => ({
    mutate: (arg: unknown, callbacks?: Callbacks) => {
      calls.push({ hook, boundTo, arg, ...(callbacks ? { callbacks } : {}) });
    },
  });
}

const toasts: string[] = [];
const compareFetches: number[][] = [];
vi.mock("../../src/lib/toast", () => ({ showToast: (msg: string) => toasts.push(msg) }));

let mateCache: (Exposure[] | undefined)[] = [];
const server = {
  indices: new Map<number, IndexEntry[]>(),
  assignments: new Map<number, Assignment>(),
  saving: false,
  saved: undefined as (() => void) | undefined,
  status: new Map<number, Exposure["status"]>(),
  colors: new Map<string, LabelColor | null>(),
  experiment: { name: "BL-19", description: null as string | null },
  callMembers: new Map<number, SampleCallMember[]>(),
};
const serverStatus = (sid: number, eid: number, status: Exposure["status"]) => {
  server.status.set(eid, status);
  calls.push({ hook: "op:setStatus", boundTo: 0, arg: [sid, eid, status] });
};
const serverColor = (label: string, color: LabelColor | null) => {
  server.colors.set(label.toLowerCase(), color);
  calls.push({ hook: "op:setLabelColor", boundTo: 0, arg: [label, color] });
};
const editOp = (name: string) => (...args: unknown[]) => calls.push({ hook: `op:${name}`, boundTo: 0, arg: args });
const editOps = {
  addPhase: editOp("addPhase"),
  removePhase: editOp("removePhase"),
  setCallState: editOp("setCallState"),
  setStatus: serverStatus,
  addPeak: editOp("addPeak"),
  removePeak: editOp("removePeak"),
  setPeakExcluded: editOp("setPeakExcluded"),
  commitCustom: editOp("commitCustom"),
  deleteIndex: editOp("deleteIndex"),
  addLabel: editOp("addLabel"),
  removeLabel: editOp("removeLabel"),
  restoreLabelSuggestion: editOp("restoreLabelSuggestion"),
  setLabelColor: serverColor,
  patchExperiment: editOp("patchExperiment"),
  bridge: editOp("bridge"),
  unbridge: editOp("unbridge"),
  patchCompares: (sampleId: number, change: CompareChange) => {
    calls.push({ hook: "op:patchCompares", boundTo: 0, arg: [sampleId, change] });
    setCompares(applyCompareChange(compareStore.data, sampleId, change, state.username ?? null));
  },
  read: {
    peaks: async () => state.peaks,
    indices: async (eid: number) => server.indices.get(eid) ?? state.indices,
    knownIndices: async (eid: number) => server.indices.get(eid) ?? state.indices,
    assignment: async (eid: number) =>
      server.assignments.get(eid)
      ?? (state.assignment?.exposure_id === eid ? state.assignment : state.frameAssignments.get(eid))
      ?? { exposure_id: eid, state: "null" as const, members: [] },
    labels: async () => state.labels,
    bridges: async () => state.bridges,
    call: async (sid: number) => ({ sample_id: sid, state: "indexed" as const, phase: null, call_phases: [], members: server.callMembers.get(sid) ?? (sid === 10 ? [{ index_id: 2, exposure_id: 101, phase: "Pn3m" }] : []) }),
    exposures: async () => state.exposures.map((e) => (server.status.has(e.id) ? { ...e, status: server.status.get(e.id)! } : e)),
    labelColors: async () => {
      const merged = new Map(state.labelColors.map((e) => [e.label.toLowerCase(), e]));
      for (const [k, color] of server.colors) {
        if (color === null) merged.delete(k);
        else merged.set(k, { label: merged.get(k)?.label ?? k, color });
      }
      return [...merged.values()];
    },
    experiment: async () => ({ id: 1, ...server.experiment }),
    compares: async () => compareStore.data,
  },
  saving: () => server.saving,
  whenSaved: () => (server.saving ? new Promise<void>((r) => { server.saved = r; }) : Promise.resolve()),
};

vi.mock("../../src/queries", () => ({
  useEditOps: () => editOps,
  useExperiments: () => ({ data: state.experiments }),
  useCorpusSamples: ({ enabled = true }: { enabled?: boolean } = {}) => {
    state.allCorpusEnabled.push(enabled);
    return { data: enabled ? state.allCorpus : undefined };
  },
  useExperimentCorpusSamples: (experimentId: number) => {
    corpusRequests.push(experimentId);
    return { data: state.corpus, isLoading: false };
  },
  useExperiment: () => ({
    data: state.experimentLoaded ? { id: 1, name: "BL-19", last_scanned_at: "2026-11-01T00:00:00", setups: state.setups } : undefined,
    isError: false,
  }),
  useProblems: () => ({ data: state.problems }),
  useTriggerScan: () => ({ mutate: (force = false) => rescans.push(force), isError: state.rescanError }),
  useCheckForNewFilesOnOpen: (experimentId: number) => { opened.push(experimentId); },
  useDeleteExperiment: (experimentId: number) => ({
    mutate: (_: unknown, cb?: { onSuccess?: () => void }) => { removals.push(experimentId); cb?.onSuccess?.(); },
    reset: () => {}, isPending: false, error: null,
  }),
  useExposures: () => ({ data: state.exposures }),
  useSamplesExposures: (ids: readonly number[]) => {
    const next = ids.map((id) => state.mateExposures.get(id));
    if (next.length !== mateCache.length || next.some((d, i) => d !== mateCache[i])) mateCache = next;
    return mateCache;
  },
  useTrace: () => ({
    data: state.trace, error: state.traceErr, isError: state.traceError, isFetching: state.traceFetching,
    fetchStatus: state.traceFetching ? "fetching" : "idle", refetch: traceRefetch,
  }),
  usePeaks: () => ({ data: state.peaks, fetchStatus: "idle" }),
  useIndices: () => ({ data: state.indices, fetchStatus: "idle" }),
  useAssignment: () => ({ data: state.assignment, fetchStatus: "idle" }),
  useSampleCall: () => ({ data: state.sampleCall }),
  useFrameAssignments: (ids: number[]) => ids.map((id) => ({ data: state.frameAssignments.get(id) })),
  useRemoveFramePhase: () => ({
    mutate: (exposureId: number, indexId: number) => calls.push({ hook: "removeFrame", boundTo: exposureId, arg: indexId }),
  }),
  useBridges: () => ({ data: state.bridges }),
  useDetectorSizes: (urls: readonly string[]) => {
    state.sizeRequests.push(...urls);
    return new Map(state.detectorSize ? urls.map((u) => [u, state.detectorSize!] as const) : []);
  },
  useCompares: () => ({ data: useSyncExternalStore(subscribeCompares, () => compareStore.data) }),
  usePatchCompares: () => ({
    mutate: (arg: { sampleId: number; change: CompareChange }) => {
      calls.push({ hook: "patchCompares", boundTo: arg.sampleId, arg: arg.change });
      setCompares(applyCompareChange(compareStore.data, arg.sampleId, arg.change, state.username ?? null));
    },
  }),
  useBridgeSamples: recorder("bridge"),
  useUnbridgeSample: recorder("unbridge"),
  useLabels: () => ({ data: state.labels }),
  useAddLabel: recorder("addLabel"),
  useRemoveLabel: recorder("removeLabel"),
  useLabelColors: () => ({ data: state.labelColors }),
  useSetLabelColor: () => ({
    mutate: (arg: { label: string; color: LabelColor | null }) => {
      server.colors.set(arg.label.toLowerCase(), arg.color);
      calls.push({ hook: "setLabelColor", boundTo: 0, arg });
    },
  }),
  useAddAssignmentPhase: recorder("add"),
  useRemoveAssignmentPhase: recorder("remove"),
  useSetSampleCallState: recorder("state"),
  useSetExposureStatus: (boundTo: number) => ({
    mutate: (arg: { exposureId: number; status: Exposure["status"] }, callbacks?: Callbacks) => {
      server.status.set(arg.exposureId, arg.status);
      calls.push({ hook: "status", boundTo, arg, ...(callbacks ? { callbacks } : {}) });
    },
  }),
  useAddPeak: recorder("addPeak"),
  useRemovePeak: recorder("removePeak"),
  useSetPeakExcluded: recorder("excludePeak"),
  useCommitCustomIndex: (boundTo: number) => ({
    mutate: (phase: string, basis: number, ratios: number[]) => calls.push({ hook: "custom", boundTo, arg: { phase, basis, ratios } }),
  }),
}));

const prefetchTargets: PrefetchTarget[] = [];
vi.mock("../../src/workspace/prefetch", () => ({
  usePrefetchNeighbours: (target: PrefetchTarget) => {
    prefetchTargets.push(target);
  },
}));

vi.mock("../../src/workspace/useCompareTraces", () => ({
  useCompareTraces: (ids: readonly number[]) => {
    compareFetches.push([...ids]);
    return new Map(ids.flatMap((id) => (state.compareTraces.has(id) ? [[id, state.compareTraces.get(id)!] as const] : [])));
  },
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
  const row = screen.getByRole("dialog").querySelector("[data-palette-item='view']");
  const name = row?.children[1]?.textContent ?? "";
  fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Escape" });
  return name;
}

function runFromPalette(...steps: string[]): void {
  fireEvent.keyDown(window, { key: "k", ctrlKey: true });
  for (const step of steps) {
    const input = screen.getByRole("combobox");
    fireEvent.change(input, { target: { value: step } });
    fireEvent.keyDown(input, { key: "Enter" });
  }
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
  return { id, experiment_id: 1, name, notes: null, q_units: "A-1", ...extra } as CorpusSample;
}

function exposure(id: number, extra: Partial<Exposure> = {}): Exposure {
  return makeExposure({ id, sample_id: 10, filename: `f${id}.dat`, ...extra });
}

function index(id: number, phase: string, score: number): IndexEntry {
  return {
    id, exposure_id: 101, phase, basis: 0.1, score, r_squared: 0.99, lattice_d: 100,
    status: "candidate", kind: "auto", inputs_hash: null, peaks: [], predicted_q: [0.1, 0.14],
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
  window.localStorage.removeItem("himalaya-ui:theme");
  window.localStorage.removeItem("himalaya-ui:workspace-plot-view");
});

beforeEach(() => {
  Element.prototype.scrollIntoView = vi.fn();
  calls.length = 0;
  toasts.length = 0;
  compareFetches.length = 0;
  prefetchTargets.length = 0;
  state.username = "jc";
  state.experimentLoaded = true;
  state.experiments = [{ id: 1, name: "BL-19" }, { id: 2, name: "Beamtime 12" }];
  state.allCorpus = [];
  state.allCorpusEnabled = [];
  state.setups = [];
  state.detectorSize = undefined;
  state.sizeRequests = [];
  state.corpus = [
    sample(10, "A1", { assignment_state: "indexed", phase: "Pn3m" }),
    sample(11, "A2"),
    sample(12, "A3", { assignment_state: "form_factor" }),
    sample(13, "A4"),
    { ...sample(99, "other"), experiment_id: 2 },
  ];
  state.exposures = [exposure(100), exposure(101)];
  state.mateExposures = new Map();
  state.indices = [index(1, "Im3m", 0.4), index(2, "Pn3m", 0.9)];
  state.assignment = { exposure_id: 101, state: "indexed", members: [2] };
  state.sampleCall = undefined;
  state.frameAssignments = new Map();
  state.trace = undefined;
  state.traceError = false;
  state.traceErr = null;
  state.traceFetching = false;
  state.peaks = [];
  state.compareTraces = new Map();
  state.bridges = [];
  compareStore.data = EMPTY_COMPARES;
  localStorage.removeItem("himalaya-ui:recent:comparison:1");
  localStorage.removeItem("himalaya-ui:recent:label");
  state.labels = [];
  state.problems = { ingest_error: null, exposures: [] };
  state.labelColors = [];
  state.ingest = undefined;
  state.rescanError = false;
  server.indices = new Map();
  server.assignments = new Map();
  server.callMembers = new Map();
  server.saving = false;
  server.saved = undefined;
  server.status = new Map();
  server.colors = new Map();
  server.experiment = { name: "BL-19", description: null };
  rescans.length = 0;
  removals.length = 0;
  opened.length = 0;
});

describe("WorkspacePage", () => {
  it("prefetches around the open sample in the listed order, once its own data has settled", () => {
    state.traceFetching = true;
    const { rerender } = renderAt("/experiments/1/workspace/11");
    const loading = prefetchTargets.at(-1)!;
    expect(loading.samples.map((s) => s.id)).toEqual([10, 11, 12, 13]);
    expect(loading.currentId).toBe(11);
    expect(loading.exposures.map((e) => e.id)).toEqual([100, 101]);
    expect(loading.shownId).toBe(101);
    expect(loading.settled).toBe(false);

    state.traceFetching = false;
    rerender(page("/experiments/1/workspace/11"));
    expect(prefetchTargets.at(-1)!.settled).toBe(true);
  });

  it("holds the sample list until the setups are known, so it paints once in its final order", () => {
    state.experimentLoaded = false;
    state.setups = [
      { id: 2, name: "setup2", folder: "setup2", included: true },
      { id: 1, name: "setup1", folder: "setup1", included: false },
    ];
    state.corpus = [
      { ...sample(10, "A1"), setup_id: 1 },
      { ...sample(11, "B1"), setup_id: 2 },
      { ...sample(12, "B2"), setup_id: 2 },
    ];
    const { rerender } = renderAt("/experiments/1/workspace");
    expect(within(screen.getByTestId("workspace-samples")).queryAllByRole("button")).toHaveLength(0);
    expect(screen.getByText("Loading samples…")).toBeInTheDocument();
    state.experimentLoaded = true;
    rerender(page("/experiments/1/workspace"));
    const rows = within(screen.getByTestId("workspace-samples")).getAllByRole("button");
    expect(rows.map((r) => r.getAttribute("data-sample-id"))).toEqual(["11", "12"]);
  });

  it("lists only this experiment's samples and marks unreviewed ones", () => {
    corpusRequests.length = 0;
    renderAt("/experiments/1/workspace/10");
    expect(new Set(corpusRequests)).toEqual(new Set([1]));
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
    expect(screen.getAllByRole("option")[0]).toHaveAccessibleName("Switch user, jc");
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
    state.exposures = [exposure(100), exposure(101, { image_path: "/x.tif" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: " " });
    fireEvent.keyDown(screen.getByRole("dialog", { name: "Detector image" }), { key: " ", repeat: true });
    expect(screen.getByRole("dialog", { name: "Detector image" })).toBeInTheDocument();
  });

  it("Space enlarges the shown frame's detector image, and Space or Esc closes it", () => {
    state.exposures = [exposure(100), exposure(101, { image_path: "/x.tif", image_version: "7" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: " " });
    const dialog = screen.getByRole("dialog", { name: "Detector image" });
    expect(dialog.querySelector("img:not([aria-hidden])")!.getAttribute("src")).toBe("/api/exposures/101/image?v=7");
    expect(screen.getByTestId("workspace-detector-caption")).toHaveTextContent("Image 2 of 2 · f101.dat");
    fireEvent.keyDown(dialog, { key: " " });
    expect(screen.queryByRole("dialog", { name: "Detector image" })).toBeNull();
    fireEvent.keyDown(window, { key: " " });
    fireEvent.keyDown(screen.getByRole("dialog", { name: "Detector image" }), { key: "Escape" });
    expect(screen.queryByRole("dialog", { name: "Detector image" })).toBeNull();
  });

  it("flips through frames inside the enlarged view, and dropping one moves to the nearest kept frame without closing it", () => {
    state.exposures = [
      exposure(100, { image_path: "/a.tif" }),
      exposure(101, { image_path: "/b.tif" }),
      exposure(102, { image_path: "/c.tif" }),
    ];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: " " });
    const dialog = screen.getByRole("dialog", { name: "Detector image" });
    expect(dialog.querySelectorAll("img")).toHaveLength(3);
    const caption = screen.getByTestId("workspace-detector-caption");
    expect(caption).toHaveTextContent("Image 3 of 3 · f102.dat");
    expect(dialog.querySelector("img:not([aria-hidden])")!.getAttribute("data-exposure-id")).toBe("102");
    fireEvent.keyDown(dialog, { key: "ArrowLeft" });
    expect(caption).toHaveTextContent("Image 2 of 3 · f101.dat");
    fireEvent.keyDown(dialog, { key: "ArrowLeft" });
    expect(caption).toHaveTextContent("Image 1 of 3");
    fireEvent.keyDown(dialog, { key: "x" });
    fireEvent.keyDown(dialog, { key: "x", repeat: true });
    expect(calls.filter((c) => c.hook === "status").map((c) => c.arg)).toEqual([{ exposureId: 100, status: "rejected" }]);
    expect(screen.getByRole("dialog", { name: "Detector image" })).toBeInTheDocument();
    expect(caption).toHaveTextContent("Image 2 of 3");
    expect(currentFrame()).toBe("Image 2");
  });

  it("does not drop from the enlarged view before a user is set", () => {
    state.username = undefined;
    state.exposures = [exposure(100, { image_path: "/a.tif" })];
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

  it("says when the last ingest failed, and Details lists the error", () => {
    state.ingest = { status: "failed", processed: 0, total: 0 };
    state.problems = { ingest_error: "SystemError: opening file", exposures: [] };
    renderAt("/experiments/1/workspace/10");
    const notice = screen.getByTestId("workspace-notice");
    expect(within(notice).getByTestId("workspace-ingest-failed")).toHaveTextContent("The last ingest failed.");
    fireEvent.click(within(notice).getByRole("button", { name: "Details" }));
    expect(screen.getByTestId("problems-ingest-error")).toHaveTextContent("SystemError: opening file");
  });

  it("opens the problem list from the palette", () => {
    state.problems = { ingest_error: null, exposures: [
      { exposure_id: 3, filename: "C_0_001", sample_id: 12, sample_name: "B1", problem: "failed", reason: "bad header" },
    ] };
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    const input = screen.getByLabelText("Command or sample");
    fireEvent.change(input, { target: { value: "show ingest problems" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(screen.getByTestId("workspace-problems")).toHaveTextContent("bad header");
  });

  it("forces the rescan after a failed ingest so it runs even with nothing new", () => {
    state.ingest = { status: "failed", processed: 0, total: 0 };
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    const input = screen.getByLabelText("Command or sample");
    fireEvent.change(input, { target: { value: "rescan" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(rescans).toEqual([true]);
  });

  it("removes the experiment from the palette after a confirmation, then goes home", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    const input = screen.getByLabelText("Command or sample");
    fireEvent.change(input, { target: { value: "remove this experiment" } });
    fireEvent.keyDown(input, { key: "Enter" });
    const dialog = screen.getByTestId("workspace-remove-experiment");
    expect(dialog).toHaveTextContent("Remove BL-19?");
    expect(removals).toEqual([]);
    fireEvent.click(within(dialog).getByRole("button", { name: "Remove experiment" }));
    expect(removals).toEqual([1]);
    expect(screen.getByTestId("where")).toHaveTextContent(/^\/experiments$/);
  });

  it("lists removing the experiment only once it is typed, and goes back to all experiments from the palette", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    expect(screen.queryByRole("option", { name: "Remove this experiment" })).toBeNull();
    const input = screen.getByLabelText("Command or sample");
    fireEvent.change(input, { target: { value: "back to all" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(screen.getByTestId("where")).toHaveTextContent(/^\/experiments$/);
  });

  it("offers no remove command while an ingest runs", () => {
    state.ingest = { status: "analyzing", stage: "analyzing", processed: 1, total: 2 };
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "remove this experiment" } });
    expect(screen.queryAllByRole("option")).toHaveLength(0);
  });

  it("says the ingest failed when the request to start it fails", () => {
    state.rescanError = true;
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-ingest-failed")).toBeInTheDocument();
  });

  it("shows nothing above the list when nothing is running and nothing is wrong", () => {
    renderAt("/experiments/1/workspace/10");
    expect(screen.queryByTestId("workspace-notice")).toBeNull();
  });

  it("counts frames waiting for a trace and failed analyses, and a problem opens its sample", () => {
    state.problems = {
      ingest_error: null,
      exposures: [
        { exposure_id: 1, filename: "A_0_001", sample_id: 11, sample_name: "A2", problem: "no_trace", reason: "no trace" },
        { exposure_id: 2, filename: "B_0_001", sample_id: 11, sample_name: "A2", problem: "no_trace", reason: "no trace" },
        { exposure_id: 3, filename: "C_0_001", sample_id: 12, sample_name: "B1", problem: "failed", reason: "bad header" },
      ],
    };
    renderAt("/experiments/1/workspace/10");
    const notice = screen.getByTestId("workspace-notice");
    expect(within(notice).getByTestId("workspace-ingest-status")).toHaveTextContent("2 images waiting for a trace · 1 analysis failed");
    fireEvent.click(within(notice).getByRole("button", { name: "Details" }));
    const rows = screen.getAllByTestId("problem-row");
    expect(rows).toHaveLength(3);
    expect(rows[2]).toHaveTextContent("bad header");
    fireEvent.click(within(rows[2]!).getByRole("button"));
    expect(screen.getByTestId("where")).toHaveTextContent("/experiments/1/workspace/12");
    expect(screen.queryByTestId("workspace-problems")).toBeNull();
  });

  it("shows ingest progress by stage and hides the Rescan command while it runs", () => {
    state.ingest = { status: "analyzing", stage: "analyzing", processed: 5, total: 1200 };
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-ingest-status")).toHaveTextContent("Analyzing 5 of 1,200");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "rescan" } });
    expect(screen.queryAllByRole("option")).toHaveLength(0);
  });

  it("says samples are on their way while the first ingest has found none yet", () => {
    state.corpus = [];
    state.ingest = { status: "scanning", stage: "discovery", processed: 0, total: 0 };
    renderAt("/experiments/1/workspace");
    expect(screen.getByText("No samples yet. They appear as the ingest reads the data folder.")).toBeInTheDocument();
    expect(screen.queryByText("No samples match.")).toBeNull();
  });

  it("says no samples match only when a filter is typed", () => {
    renderAt("/experiments/1/workspace/10");
    expect(screen.queryByText("No samples match.")).toBeNull();
    fireEvent.change(screen.getByRole("searchbox", { name: "Filter samples" }), { target: { value: "zzz" } });
    expect(screen.getByText("No samples match.")).toBeInTheDocument();
  });

  it("checks this experiment's folder for new files when it opens", () => {
    renderAt("/experiments/1/workspace/10");
    expect(opened.length).toBeGreaterThan(0);
    expect(new Set(opened)).toEqual(new Set([1]));
    expect(rescans).toEqual([]);
  });

  it("rescans from the palette, and offers no configuration or grouping review", () => {
    renderAt("/experiments/1/workspace/10");
    const run = (label: string): void => {
      fireEvent.keyDown(window, { key: "k", ctrlKey: true });
      const input = screen.getByLabelText("Command or sample");
      fireEvent.change(input, { target: { value: label } });
      fireEvent.keyDown(input, { key: "Enter" });
    };
    run("rescan");
    expect(rescans).toEqual([false]);
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "configure" } });
    expect(screen.queryAllByRole("option")).toHaveLength(0);
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "grouping" } });
    expect(screen.queryAllByRole("option")).toHaveLength(0);
  });

  it("updates the list and the count when the corpus status changes", () => {
    const { rerender } = renderAt("/experiments/1/workspace/11");
    state.corpus = state.corpus.map((s) => (s.id === 11 ? { ...s, assignment_state: "indexed", phase: "Im3m" } : s));
    rerender(page("/experiments/1/workspace/11"));
    const row = within(screen.getByTestId("workspace-samples")).getByRole("button", { current: true });
    expect(within(row).getByTestId("workspace-status")).toHaveTextContent("Im3m");
    expect(screen.getByTestId("workspace-progress")).toHaveTextContent("3 of 4 reviewed");
  });

  it("lists every phase of a coexistence call, as the mock-up does", () => {
    const { rerender } = renderAt("/experiments/1/workspace/11");
    state.corpus = state.corpus.map((s) => (s.id === 11 ? { ...s, assignment_state: "indexed", phase: "Im3m", call_phases: ["Im3m", "Pn3m"] } : s));
    rerender(page("/experiments/1/workspace/11"));
    const row = within(screen.getByTestId("workspace-samples")).getByRole("button", { current: true });
    expect(within(row).getByTestId("workspace-status")).toHaveTextContent("Im3m + Pn3m");
  });

  it("Z jumps to the next unreviewed sample, like Gradescope's next ungraded", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "z" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
  });

  it("J and R are no longer shortcuts", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(screen.getByRole("dialog").querySelector("[data-sample-id='11']")!);
    expect(screen.getByRole("button", { name: "A2" })).toBeInTheDocument();
    fireEvent.keyDown(window, { key: "j" });
    fireEvent.keyDown(window, { key: "r" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/10");
    expect(currentView()).toBe("Log");
    expect(screen.getByRole("button", { name: "A2" })).toBeInTheDocument();
    expect(screen.queryByRole("button", { name: "A2, hidden" })).toBeNull();
  });

  it("⌘K opens another experiment, and finds its samples once something is typed", () => {
    state.allCorpus = [...state.corpus, { ...sample(40, "lipid_A"), experiment_id: 2 }];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    expect(state.allCorpusEnabled.every((e) => !e)).toBe(true);
    const input = screen.getByRole("combobox");
    fireEvent.change(input, { target: { value: "BL-19" } });
    expect(screen.queryByRole("option", { name: /Open BL-19/ })).toBeNull();
    fireEvent.change(input, { target: { value: "lipid" } });
    expect(state.allCorpusEnabled.at(-1)).toBe(true);
    expect(screen.getByRole("option", { name: /Go to lipid_A/ })).toHaveTextContent("Beamtime 12");
    fireEvent.change(input, { target: { value: "Beamtime 12" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/2");
  });

  it("/ puts the cursor in the filter, where letters type, and ? opens the shortcuts only while the filter is empty", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "/" });
    const filter = screen.getByRole("searchbox", { name: "Filter samples" });
    expect(document.activeElement).toBe(filter);
    fireEvent.keyDown(filter, { key: "z" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/10");
    fireEvent.keyDown(filter, { key: "?" });
    const help = screen.getByRole("dialog", { name: "Keyboard shortcuts" });
    fireEvent.keyDown(help, { key: "Escape" });
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
    filter.focus();
    fireEvent.change(filter, { target: { value: "a" } });
    fireEvent.keyDown(filter, { key: "?" });
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
  });

  it("/ selects an earlier filter, so the next query replaces it", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "/" });
    const filter = screen.getByRole("searchbox", { name: "Filter samples" }) as HTMLInputElement;
    fireEvent.change(filter, { target: { value: "buffer" } });
    filter.blur();
    fireEvent.keyDown(window, { key: "/" });
    expect(document.activeElement).toBe(filter);
    expect([filter.selectionStart, filter.selectionEnd]).toEqual([0, "buffer".length]);
  });

  it("Esc leaves a page field that has no Esc handling of its own", () => {
    renderAt("/experiments/1/workspace/10");
    const field = document.createElement("input");
    document.body.appendChild(field);
    field.focus();
    fireEvent.keyDown(field, { key: "Escape" });
    expect(document.activeElement).not.toBe(field);
    field.remove();
  });

  it("? types into the empty palette field instead of opening the shortcuts over it", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", metaKey: true });
    const input = screen.getByLabelText("Command or sample");
    expect(input).toHaveValue("");
    fireEvent.keyDown(input, { key: "?" });
    expect(screen.queryByRole("dialog", { name: "Keyboard shortcuts" })).toBeNull();
  });

  it("⌘/Ctrl+0 in a field leaves the zoom alone", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    const plot = screen.getByTestId("workspace-trace");
    const xTicks = (): string =>
      Array.from(plot.ownerDocument.querySelectorAll('[data-testid="workspace-trace"] [data-role="axis-bottom"] text')).map((t) => t.textContent).join(" ");
    const full = xTicks();
    for (let i = 0; i < 5; i++) fireEvent.wheel(plot.parentElement!, { deltaY: -300, clientX: 500 });
    const zoomed = xTicks();
    expect(zoomed).not.toBe(full);
    const filter = screen.getByRole("searchbox", { name: "Filter samples" });
    filter.focus();
    fireEvent.keyDown(filter, { key: "0", ctrlKey: true });
    expect(xTicks()).toBe(zoomed);
    fireEvent.keyDown(window, { key: "0", ctrlKey: true });
    expect(xTicks()).toBe(full);
  });

  it("page keys stay off while a modal is open, even before focus lands in it", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "?" });
    expect(screen.getByRole("dialog", { name: "Keyboard shortcuts" })).toBeInTheDocument();
    fireEvent.keyDown(document.body, { key: "z" });
    fireEvent.keyDown(document.body, { key: "ArrowDown" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/10");
  });

  it("arrow down steps to the next sample in the list", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
  });

  it("says the beam was likely off when the server reads the shown frame as beam-off, and adds I₀'s sign as support", () => {
    state.trace = makeTrace(0.001);
    state.exposures = [exposure(100, { incident_monitor: 1.5 }), exposure(101, { incident_monitor: -0.1, beam_off: true })];
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-plot-notice")).toHaveTextContent(
      "Beam likely off: the detector recorded no counts, so the beam was off or the shutter closed. I₀ also reads below zero.",
    );
  });

  it("says the beam was likely off for a beam-off frame even when I₀ looks fine", () => {
    state.trace = makeTrace(0.001);
    state.exposures = [exposure(101, { incident_monitor: 1.5, beam_off: true })];
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-plot-notice").textContent).toBe(
      "Beam likely off: the detector recorded no counts, so the beam was off or the shutter closed.",
    );
  });

  it("leaves the beam-off call to the server, so a dark-looking trace the server reads as lit gets no caption", () => {
    state.trace = makeTrace(0.001);
    state.exposures = [exposure(101, { incident_monitor: -0.1, beam_off: false })];
    renderAt("/experiments/1/workspace/10");
    expect(screen.queryByTestId("workspace-plot-notice")).toBeNull();
  });

  it("never calls a frame with counts beam-off, whatever I₀ reads", () => {
    state.trace = makeTrace(20);
    state.exposures = [exposure(101, { incident_monitor: -0.1 })];
    renderAt("/experiments/1/workspace/10");
    expect(screen.queryByTestId("workspace-plot-notice")).toBeNull();
  });

  it("says a frame has no integrated trace instead of offering a retry, with I₀'s sign as a hint", () => {
    state.traceError = true;
    state.traceErr = new ApiError(404, ".dat file not found", null);
    state.exposures = [exposure(101, { incident_monitor: -0.1 })];
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-trace-empty")).toHaveTextContent(
      "No integrated trace for this image. I₀ reads below zero, so the beam was likely off.",
    );
    expect(screen.queryByRole("button", { name: "Retry" })).toBeNull();
  });

  it("says a frame has no integrated trace without guessing when I₀ is unknown", () => {
    state.traceError = true;
    state.traceErr = new ApiError(404, ".dat file not found", null);
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-trace-empty").textContent).toBe("No integrated trace for this image.");
  });

  it("says a frame is waiting for its trace when ingest reports it pending, matching the notice", () => {
    state.traceError = true;
    state.traceErr = new ApiError(404, ".dat file not found", null);
    state.exposures = [exposure(101, { incident_monitor: -0.1 })];
    state.problems = { ingest_error: null, exposures: [
      { exposure_id: 101, filename: "f101", sample_id: 10, sample_name: "S", problem: "no_trace", reason: "no trace" },
    ] };
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-trace-empty").textContent).toBe("Waiting for this image's trace.");
    expect(screen.getByTestId("workspace-ingest-status")).toHaveTextContent("1 image waiting for a trace");
  });

  it("keeps saying no integrated trace when only another frame is pending", () => {
    state.traceError = true;
    state.traceErr = new ApiError(404, ".dat file not found", null);
    state.exposures = [exposure(101)];
    state.problems = { ingest_error: null, exposures: [
      { exposure_id: 999, filename: "f999", sample_id: 11, sample_name: "T", problem: "no_trace", reason: "no trace" },
    ] };
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-trace-empty").textContent).toBe("No integrated trace for this image.");
  });

  it("opens on the representative frame", () => {
    renderAt("/experiments/1/workspace/10");
    expect(currentFrame()).toBe("Image 2");
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

  it("on an unselected frame, adds the phase to that frame without selecting it", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    fireEvent.keyDown(window, { key: "2" });
    expect(calls).toEqual([{ hook: "add", boundTo: 100, arg: 1 }]);
  });

  it("X drops the shown frame and shows the nearest kept frame, so a phase key calls there", () => {
    const { rerender } = renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "x" });
    expect(calls[0]).toEqual({ hook: "status", boundTo: 10, arg: { exposureId: 101, status: "rejected" } });
    state.exposures = [exposure(100), exposure(101, { status: "rejected" })];
    rerender(page("/experiments/1/workspace/10"));
    expect(currentFrame()).toBe("Image 1");
    fireEvent.keyDown(window, { key: "2" });
    expect(calls[1]).toEqual({ hook: "add", boundTo: 100, arg: 1 });
  });

  it("X on a dropped frame restores it and keeps it in view", () => {
    state.exposures = [exposure(100), exposure(101, { status: "rejected" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(currentFrame()).toBe("Image 2, dropped");
    fireEvent.keyDown(window, { key: "x" });
    expect(calls[0]).toEqual({ hook: "status", boundTo: 10, arg: { exposureId: 101, status: null } });
    expect(currentFrame()).toBe("Image 2, dropped");
  });

  it("X keeps the dropped frame in view when no other frame is kept", () => {
    state.exposures = [exposure(100, { status: "rejected" }), exposure(101)];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "x" });
    expect(calls[0]).toEqual({ hook: "status", boundTo: 10, arg: { exposureId: 101, status: "rejected" } });
    expect(currentFrame()).toBe("Image 2");
  });

  it("a card's Drop button drops that frame and leaves the view on the shown one", () => {
    state.exposures = [exposure(100), exposure(101), exposure(102)];
    renderAt("/experiments/1/workspace/10");
    expect(currentFrame()).toBe("Image 3");
    const cards = within(screen.getByTestId("workspace-exposures")).getAllByRole("listitem");
    fireEvent.click(within(cards[0]!).getByRole("button", { name: "Drop" }));
    expect(calls[0]).toEqual({ hook: "status", boundTo: 10, arg: { exposureId: 100, status: "rejected" } });
    expect(currentFrame()).toBe("Image 3");
  });

  it("the shown frame's own Drop button moves to the nearest kept frame", () => {
    state.exposures = [exposure(100), exposure(101), exposure(102)];
    renderAt("/experiments/1/workspace/10");
    const cards = within(screen.getByTestId("workspace-exposures")).getAllByRole("listitem");
    fireEvent.click(within(cards[2]!).getByRole("button", { name: "Drop" }));
    expect(calls[0]).toEqual({ hook: "status", boundTo: 10, arg: { exposureId: 102, status: "rejected" } });
    expect(currentFrame()).toBe("Image 2");
  });

  it("holding X drops only the shown frame, not a run of frames", () => {
    state.exposures = [exposure(100), exposure(101), exposure(102)];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "x" });
    fireEvent.keyDown(window, { key: "x", repeat: true });
    fireEvent.keyDown(window, { key: "x", repeat: true });
    expect(calls.filter((c) => c.hook === "status")).toHaveLength(1);
    expect(currentFrame()).toBe("Image 2");
  });

  describe("undo and redo", () => {
    const ops = () => calls.filter((c) => c.hook.startsWith("op:")).map((c) => [c.hook.slice(3), ...(c.arg as unknown[])]);
    const settle = () => act(async () => { for (let i = 0; i < 10; i++) await Promise.resolve(); });
    const press = async (init: KeyboardEventInit) => {
      fireEvent.keyDown(window, init);
      await settle();
    };
    const undo = () => press({ key: "z", ctrlKey: true });
    const redo = () => press({ key: "Z", ctrlKey: true, shiftKey: true });
    const onServer = (eid: number, a: Partial<Assignment>) =>
      server.assignments.set(eid, { exposure_id: eid, state: "indexed", members: [], ...a });

    it("ticks a Bonnet group whole with its number key and undoes it as one edit", async () => {
      const peaksAt = (qs: number[]) => qs.map((q, i) => ({ peak_id: i + 1, ratio_position: i + 1, residual: 0, q_observed: q }));
      state.indices = [
        { ...index(1, "Im3m", 0.9), lattice_d: 127.9, peaks: peaksAt([0.05, 0.07, 0.09, 0.11]) },
        { ...index(2, "Pn3m", 0.6), lattice_d: 100, peaks: peaksAt([0.06, 0.08, 0.1, 0.12]) },
      ];
      state.assignment = { exposure_id: 101, state: "null", members: [] };
      renderAt("/experiments/1/workspace/10");
      expect(screen.getByTestId("workspace-bonnet-group")).toHaveAttribute("aria-pressed", "false");
      fireEvent.keyDown(window, { key: "1" });
      expect(calls).toEqual([
        { hook: "add", boundTo: 101, arg: 1 },
        { hook: "add", boundTo: 101, arg: 2 },
      ]);
      onServer(101, { members: [1, 2] });
      await undo();
      expect(ops()).toEqual([["removePhase", 101, 2], ["removePhase", 101, 1]]);
      expect(toasts.at(-1)).toBe("Undone: ticked Im3m + Pn3m on A1");
    });

    it("adds only the missing member when the header of a partly ticked group is clicked, and undoes only that", async () => {
      const peaksAt = (qs: number[]) => qs.map((q, i) => ({ peak_id: i + 1, ratio_position: i + 1, residual: 0, q_observed: q }));
      state.indices = [
        { ...index(1, "Im3m", 0.9), lattice_d: 127.9, peaks: peaksAt([0.05, 0.07, 0.09, 0.11]) },
        { ...index(2, "Pn3m", 0.6), lattice_d: 100, peaks: peaksAt([0.06, 0.08, 0.1, 0.12]) },
      ];
      state.assignment = { exposure_id: 101, state: "indexed", members: [1] };
      renderAt("/experiments/1/workspace/10");
      const header = screen.getByTestId("workspace-bonnet-group");
      expect(header).toHaveAttribute("aria-pressed", "mixed");
      fireEvent.click(header);
      expect(calls).toEqual([{ hook: "add", boundTo: 101, arg: 2 }]);
      onServer(101, { members: [1, 2] });
      await undo();
      expect(ops()).toEqual([["removePhase", 101, 2]]);
    });

    it("unticks a whole Bonnet group when every member is in the call", () => {
      const peaksAt = (qs: number[]) => qs.map((q, i) => ({ peak_id: i + 1, ratio_position: i + 1, residual: 0, q_observed: q }));
      state.indices = [
        { ...index(1, "Im3m", 0.9), lattice_d: 127.9, peaks: peaksAt([0.05, 0.07, 0.09, 0.11]) },
        { ...index(2, "Pn3m", 0.6), lattice_d: 100, peaks: peaksAt([0.06, 0.08, 0.1, 0.12]) },
      ];
      state.assignment = { exposure_id: 101, state: "indexed", members: [1, 2] };
      renderAt("/experiments/1/workspace/10");
      fireEvent.click(screen.getByTestId("workspace-bonnet-group"));
      expect(calls).toEqual([
        { hook: "remove", boundTo: 101, arg: 1 },
        { hook: "remove", boundTo: 101, arg: 2 },
      ]);
    });

    it("Ctrl+Z after 0 clears no phase and ticks every phase it wiped, on this frame and others", async () => {
      state.sampleCall = {
        sample_id: 10, state: "indexed", phase: "Pn3m", call_phases: ["Pn3m", "Im3m"],
        members: [{ index_id: 2, exposure_id: 101, phase: "Pn3m" }, { index_id: 7, exposure_id: 100, phase: "Im3m" }],
      };
      server.indices.set(100, [{ ...index(7, "Im3m", 0.5), exposure_id: 100 }]);
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "0" });
      expect(calls[0]).toEqual({ hook: "state", boundTo: 10, arg: "form_factor" });
      onServer(101, { state: "form_factor" });
      onServer(100, { state: "form_factor" });
      await undo();
      expect(ops()).toEqual([["setCallState", 10, "null"], ["addPhase", 101, 2], ["addPhase", 100, 7]]);
      expect(toasts.at(-1)).toBe("Undone: marked no phase on A1");
    });

    it("Ctrl+Shift+Z redoes what was undone", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "1" });
      onServer(101, { members: [] });
      await undo();
      onServer(101, { members: [2] });
      await redo();
      expect(ops()).toEqual([["addPhase", 101, 2], ["removePhase", 101, 2]]);
      expect(toasts.at(-1)).toBe("Redone: unticked Pn3m on A1");
    });

    it("unticks a phase by its lattice after a peak edit renumbered the candidates", async () => {
      state.assignment = { exposure_id: 101, state: "null", members: [] };
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "1" });
      server.indices.set(101, [{ ...index(33, "Pn3m", 0.9), basis: 0.101 }, index(34, "Im3m", 0.4)]);
      onServer(101, { members: [33] });
      await undo();
      expect(ops()).toEqual([["removePhase", 101, 33]]);
      onServer(101, { members: [] });
      await redo();
      expect(ops()).toEqual([["removePhase", 101, 33], ["addPhase", 101, 33]]);
    });

    it("says so instead of claiming an undo when the edit is already undone", async () => {
      state.assignment = { exposure_id: 101, state: "null", members: [] };
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "1" });
      await undo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Already that way, nothing to undo: ticked Pn3m on A1");
    });

    it("waits for the last edit to save before undoing it", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      server.saving = true;
      await undo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Saving your last edit first, then undoing: dropped image 2 of A1");
      server.saving = false;
      server.saved?.();
      await settle();
      expect(ops()).toEqual([["setStatus", 10, 101, null]]);
      expect(toasts.at(-1)).toBe("Undone: dropped image 2 of A1");
    });

    it("undoes a drop by restoring the frame", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      await press({ key: "z", metaKey: true });
      expect(ops()).toEqual([["setStatus", 10, 101, null]]);
      expect(toasts.at(-1)).toBe("Undone: dropped image 2 of A1");
    });

    it("leaves a frame alone when someone else restored it before the undo", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      server.status.set(101, null);
      await undo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Already that way, nothing to undo: dropped image 2 of A1");
    });

    it("leaves a label's colour alone when someone else recoloured it before the undo", async () => {
      state.labels = [{ id: 1, sample_id: 10, label: "37C", state: "confirmed", source: "user" }];
      renderAt("/experiments/1/workspace/10");
      fireEvent.contextMenu(screen.getByTestId("workspace-label-chips").querySelector("[data-label='37C']")!);
      fireEvent.click(within(screen.getByRole("dialog", { name: "Colour of label 37C" })).getByRole("radio", { name: "Teal" }));
      server.colors.set("37c", "green");
      await undo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Already that way, nothing to undo: recoloured the label 37C");
    });

    it("redoes a drop only while the frame is still how the undo left it", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      await undo();
      server.status.set(101, "rejected");
      calls.length = 0;
      await redo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Already that way, nothing to redo: dropped image 2 of A1");
    });

    it("redoes a label recolour only while the label still has the colour the undo restored", async () => {
      state.labels = [{ id: 1, sample_id: 10, label: "37C", state: "confirmed", source: "user" }];
      renderAt("/experiments/1/workspace/10");
      fireEvent.contextMenu(screen.getByTestId("workspace-label-chips").querySelector("[data-label='37C']")!);
      fireEvent.click(within(screen.getByRole("dialog", { name: "Colour of label 37C" })).getByRole("radio", { name: "Teal" }));
      await undo();
      server.colors.set("37c", "green");
      calls.length = 0;
      await redo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Already that way, nothing to redo: recoloured the label 37C");
    });

    it("redoes a rename only while the server still holds the name the undo restored", async () => {
      renderAt("/experiments/1/workspace/10");
      runFromPalette("settings", "name");
      const name = within(screen.getByRole("dialog", { name: "Name" })).getByRole("textbox", { name: "Name" });
      fireEvent.change(name, { target: { value: "BL-20" } });
      fireEvent.keyDown(name, { key: "Enter" });
      server.experiment = { name: "BL-20", description: null };
      await undo();
      server.experiment = { name: "BL-21", description: null };
      calls.length = 0;
      await redo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Already that way, nothing to redo: renamed the experiment to BL-20");
    });

    it("undoes recolouring a label while it still has the colour this edit set", async () => {
      state.labels = [{ id: 1, sample_id: 10, label: "37C", state: "confirmed", source: "user" }];
      renderAt("/experiments/1/workspace/10");
      fireEvent.contextMenu(screen.getByTestId("workspace-label-chips").querySelector("[data-label='37C']")!);
      fireEvent.click(within(screen.getByRole("dialog", { name: "Colour of label 37C" })).getByRole("radio", { name: "Teal" }));
      await undo();
      expect(ops()).toEqual([["setLabelColor", "37C", null]]);
      expect(toasts.at(-1)).toBe("Undone: recoloured the label 37C");
    });

    it("ignores a held-down Ctrl+Z", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      await press({ key: "z", ctrlKey: true, repeat: true });
      expect(ops()).toEqual([]);
      await undo();
      expect(ops()).toEqual([["setStatus", 10, 101, null]]);
    });

    it("goes back to the sample the undone edit was made on", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      fireEvent.keyDown(window, { key: "ArrowDown" });
      expect(screen.getByTestId("where")).toHaveTextContent("/experiments/1/workspace/11");
      await undo();
      expect(screen.getByTestId("where")).toHaveTextContent("/experiments/1/workspace/10");
      expect(ops()).toEqual([["setStatus", 10, 101, null]]);
    });

    it("undoes peak edits made on the plot, finding each peak by q on the server", async () => {
      state.trace = makeTrace(1);
      state.peaks = peakList(101, { auto: [0.1, 0.12], manual: [0.2] });
      renderAt("/experiments/1/workspace/10");
      clickPlotAt(0.1);
      clickPlotAt(0.2);
      state.peaks = peakList(101, { auto: [{ q: 0.1, excluded: true }, 0.12] });
      await undo();
      expect(ops()).toEqual([["addPeak", 101, 0.2]]);
      state.peaks = [...peakList(101, { auto: [{ q: 0.1, excluded: true }, 0.12] }), { ...peak(9, 0.2, "manual") }];
      await undo();
      expect(ops()).toEqual([["addPeak", 101, 0.2], ["setPeakExcluded", 101, 1, false]]);
    });

    it("undoing a peak edit that dropped the call ticks the phase again on the re-fitted candidate", async () => {
      state.trace = makeTrace(1);
      state.peaks = peakList(101, { auto: [0.1, 0.12] });
      state.indices = [index(2, "Pn3m", 0.9), index(3, "Im3m", 0.5)];
      state.assignment = { exposure_id: 101, state: "indexed", members: [2] };
      renderAt("/experiments/1/workspace/10");
      clickPlotAt(0.1);
      state.peaks = peakList(101, { auto: [{ q: 0.1, excluded: true }, 0.12] });
      server.indices.set(101, [{ ...index(33, "Pn3m", 0.9), basis: 0.102 }, { ...index(34, "Im3m", 0.5), basis: 0.3 }]);
      onServer(101, { state: "null", members: [] });
      await undo();
      expect(ops()).toEqual([["setPeakExcluded", 101, 1, false], ["addPhase", 101, 33]]);
      expect(toasts.at(-1)).toBe("Undone: turned off a peak on image 2 of A1");
    });

    it("does not record removing a peak that is still saving", async () => {
      state.trace = makeTrace(1);
      state.peaks = [peak(-5, 0.2, "manual")];
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(screen.getByRole("button", { name: /Manual peak at q = 0\.2000/ }), { key: "Enter" });
      expect(calls.filter((c) => c.hook === "removePeak")).toEqual([]);
      await undo();
      expect(toasts.at(-1)).toBe("Nothing to undo.");
    });

    it("leaves Ctrl+Z alone while the custom index is open", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "1" });
      fireEvent.keyDown(window, { key: "i" });
      expect(screen.getByTestId("workspace-custom-index")).toBeInTheDocument();
      await undo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Close the custom index first (Esc).");
    });

    it("leaves Ctrl+Z to a text field being typed in", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      const filter = screen.getByRole("searchbox", { name: "Filter samples" });
      filter.focus();
      fireEvent.keyDown(filter, { key: "z", ctrlKey: true });
      await settle();
      expect(ops()).toEqual([]);
    });

    it("says so when there is nothing to undo", async () => {
      renderAt("/experiments/1/workspace/10");
      await undo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Nothing to undo.");
    });

    it("does not undo before a user is set", async () => {
      const { rerender } = renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      state.username = undefined;
      rerender(page("/experiments/1/workspace/10"));
      await undo();
      expect(ops()).toEqual([]);
    });

    it("undoes removing a label from the label row", async () => {
      state.labels = [{ id: 1, sample_id: 10, label: "Buffer", state: "confirmed", source: "user" }];
      renderAt("/experiments/1/workspace/10");
      fireEvent.click(screen.getByRole("button", { name: "Remove label Buffer" }));
      expect(calls.some((c) => c.hook === "removeLabel")).toBe(true);
      state.labels = [];
      await undo();
      expect(ops()).toEqual([["addLabel", 10, "Buffer"]]);
    });

    it("undoes a label added from the palette", async () => {
      renderAt("/experiments/1/workspace/10");
      runFromPalette("label", "Calibrant");
      state.labels = [{ id: 1, sample_id: 10, label: "Calibrant", state: "confirmed", source: "user" }];
      await undo();
      expect(ops()).toEqual([["removeLabel", 10, "Calibrant"]]);
      expect(toasts.at(-1)).toBe("Undone: added the label Calibrant on A1");
    });

    it("undoes confirming a suggested label by putting it back to suggested, not deleting it", async () => {
      state.labels = [{ id: 1, sample_id: 10, label: "Buffer", state: "suggested", source: "ingest" }];
      renderAt("/experiments/1/workspace/10");
      fireEvent.click(screen.getByRole("button", { name: "Confirm label Buffer" }));
      state.labels = [{ id: 1, sample_id: 10, label: "Buffer", state: "confirmed", source: "user" }];
      await undo();
      expect(ops()).toEqual([["restoreLabelSuggestion", 10, "Buffer"]]);
      expect(toasts.at(-1)).toBe("Undone: confirmed the suggested label Buffer on A1");
    });

    it("undoes and redoes renaming the experiment and changing its description from the palette", async () => {
      renderAt("/experiments/1/workspace/10");
      runFromPalette("settings", "name");
      const name = within(screen.getByRole("dialog", { name: "Name" })).getByRole("textbox", { name: "Name" });
      fireEvent.change(name, { target: { value: "BL-20" } });
      fireEvent.keyDown(name, { key: "Enter" });
      runFromPalette("settings", "description");
      const description = within(screen.getByRole("dialog", { name: "Description" })).getByRole("textbox", { name: "Description" });
      fireEvent.change(description, { target: { value: "Lipid screen" } });
      fireEvent.keyDown(description, { key: "Enter" });
      server.experiment = { name: "BL-20", description: "Lipid screen" };
      calls.length = 0;
      await undo();
      expect(ops()).toEqual([["patchExperiment", { description: null }]]);
      expect(toasts.at(-1)).toBe("Undone: changed the experiment description");
      server.experiment = { name: "BL-20", description: null };
      calls.length = 0;
      await undo();
      expect(ops()).toEqual([["patchExperiment", { name: "BL-19" }]]);
      expect(toasts.at(-1)).toBe("Undone: renamed the experiment to BL-20");
      server.experiment = { name: "BL-19", description: null };
      calls.length = 0;
      await redo();
      expect(ops()).toEqual([["patchExperiment", { name: "BL-20" }]]);
      expect(toasts.at(-1)).toBe("Redone: renamed the experiment to BL-20");
    });

    const renameTo = (value: string) => {
      runFromPalette("settings", "name");
      const name = within(screen.getByRole("dialog", { name: "Name" })).getByRole("textbox", { name: "Name" });
      fireEvent.change(name, { target: { value } });
      fireEvent.keyDown(name, { key: "Enter" });
    };

    it("leaves a newer rename by someone else alone when undoing a rename", async () => {
      renderAt("/experiments/1/workspace/10");
      renameTo("BL-20");
      server.experiment = { name: "BL-21", description: null };
      calls.length = 0;
      await undo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Already that way, nothing to undo: renamed the experiment to BL-20");
    });

    it("sends nothing when the server already holds the name the undo would restore", async () => {
      renderAt("/experiments/1/workspace/10");
      renameTo("BL-20");
      server.experiment = { name: "BL-19", description: null };
      calls.length = 0;
      await undo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Already that way, nothing to undo: renamed the experiment to BL-20");
    });

    it("records nothing for a Settings save that leaves the name unchanged", async () => {
      renderAt("/experiments/1/workspace/10");
      renameTo("BL-19");
      expect(ops()).toEqual([]);
      await undo();
      expect(toasts.at(-1)).toBe("Nothing to undo.");
    });

    it("undoes creating a coloured label from the palette in one step", async () => {
      renderAt("/experiments/1/workspace/10");
      runFromPalette("label", "dehydrated");
      const form = screen.getByRole("dialog", { name: "New label" });
      fireEvent.keyDown(within(form).getByRole("textbox", { name: "Name" }), { key: "Enter" });
      state.labels = [{ id: 1, sample_id: 10, label: "dehydrated", state: "confirmed", source: "user" }];
      await undo();
      expect(ops()).toEqual([["removeLabel", 10, "dehydrated"], ["setLabelColor", "dehydrated", null]]);
      expect(toasts.at(-1)).toBe("Undone: created the label dehydrated on A1");
      calls.length = 0;
      await undo();
      expect(ops()).toEqual([]);
      expect(toasts.at(-1)).toBe("Nothing to undo.");
    });

    it("undoes saying a sample is not the same sample after all", async () => {
      state.bridges = [{ sample_id: 10, group_id: 1 }, { sample_id: 11, group_id: 1 }];
      renderAt("/experiments/1/workspace/11");
      runFromPalette("not the same sample as a1");
      expect(calls).toContainEqual({ hook: "unbridge", boundTo: 1, arg: { sampleId: 11 } });
      state.bridges = [];
      await undo();
      expect(ops()).toEqual([["bridge", 11, 10, "both"]]);
      expect(toasts.at(-1)).toBe("Undone: marked A2 as not the same sample");
    });

    it("undoes and redoes adding a comparison, and skips the redo once someone else re-added it", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "c" });
      fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
      calls.length = 0;
      await undo();
      expect(ops()).toEqual([["patchCompares", 10, { remove: [11] }]]);
      expect(toasts.at(-1)).toBe("Undone: compared A1 with A2");
      expect(screen.queryByRole("button", { name: "Remove A2" })).toBeNull();
      calls.length = 0;
      await redo();
      expect(ops()).toEqual([["patchCompares", 10, { add: [11] }]]);
      expect(screen.getByRole("button", { name: "Remove A2" })).toBeInTheDocument();
      await undo();
      act(() => setCompares(applyCompareChange(compareStore.data, 10, { add: [11] }, "ana")));
      calls.length = 0;
      await redo();
      expect(ops()).toEqual([]);
    });

    it("undoing a subtraction puts the raw trace back and leaves the comparison", async () => {
      state.trace = makeTrace(1);
      state.compareTraces = new Map([[11, makeTrace(2)]]);
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "c" });
      fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
      fireEvent.click(screen.getByRole("button", { name: "A2" }), { ctrlKey: true });
      expect(screen.getByTestId("workspace-subtraction-toggle")).toBeInTheDocument();
      calls.length = 0;
      await undo();
      expect(ops()).toEqual([["patchCompares", 10, { subtraction: null }]]);
      expect(toasts.at(-1)).toBe("Undone: subtracted A2 from A1");
      expect(screen.queryByTestId("workspace-subtraction-toggle")).toBeNull();
      expect(screen.getByRole("button", { name: "Remove A2" })).toBeInTheDocument();
    });

    it("undoing removing a subtracted comparison restores it as the subtraction", async () => {
      act(() => setCompares({ compares: [{ sample_id: 10, ref_sample_id: 11, added_by: "verify" }], subtractions: [{ sample_id: 10, ref_sample_id: 11 }] }));
      renderAt("/experiments/1/workspace/10");
      fireEvent.click(screen.getByRole("button", { name: "Remove A2" }));
      calls.length = 0;
      await undo();
      expect(ops()).toEqual([["patchCompares", 10, { add: [11], subtraction: 11 }]]);
      expect(toasts.at(-1)).toBe("Undone: stopped comparing A1 with A2");
    });

    it("Ctrl+K then Enter still goes to the next unreviewed sample when there is something to undo", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      fireEvent.keyDown(window, { key: "k", ctrlKey: true });
      fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Enter" });
      await settle();
      expect(ops()).toEqual([]);
      expect(screen.getByTestId("where")).not.toHaveTextContent("/experiments/1/workspace/10");
    });

    it("lists the next undo in the command palette", async () => {
      renderAt("/experiments/1/workspace/10");
      fireEvent.keyDown(window, { key: "x" });
      fireEvent.keyDown(window, { key: "k", ctrlKey: true });
      fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "undo" } });
      fireEvent.click(screen.getByRole("option", { name: /Undo: dropped image 2 of A1/ }));
      await settle();
      expect(ops()).toEqual([["setStatus", 10, 101, null]]);
    });
  });

  it("a dropped card offers Restore and the header counts the kept frames", () => {
    state.exposures = [exposure(100, { status: "rejected" }), exposure(101)];
    renderAt("/experiments/1/workspace/10");
    const strip = screen.getByTestId("workspace-exposures");
    expect(within(within(strip).getAllByRole("listitem")[0]!).getByRole("button", { name: "Restore" })).toBeInTheDocument();
    expect(screen.getByText("1 of 2 kept")).toBeInTheDocument();
  });

  it("refuses a phase call on a dropped frame and says why", () => {
    state.exposures = [exposure(100), exposure(101, { status: "rejected" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(currentFrame()).toBe("Image 2, dropped");
    fireEvent.keyDown(window, { key: "2" });
    expect(calls).toEqual([]);
    expect(screen.getByText("This image is dropped. Restore it to call a phase.")).toBeInTheDocument();
  });

  it("ticks a phase called on another frame and unticks it there, even when that frame is dropped", () => {
    state.exposures = [exposure(100), exposure(101, { status: "rejected" })];
    state.assignment = { exposure_id: 100, state: "indexed", members: [] };
    state.sampleCall = {
      sample_id: 10, state: "indexed", phase: "Pn3m", call_phases: ["Pn3m"],
      members: [{ index_id: 7, exposure_id: 101, phase: "Pn3m" }],
    };
    renderAt("/experiments/1/workspace/10");
    expect(currentFrame()).toBe("Image 1");
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows[0]).toHaveTextContent("Pn3m");
    expect(rows[0]).toHaveAttribute("aria-pressed", "true");
    expect(rows[1]).toHaveAttribute("aria-pressed", "false");
    fireEvent.keyDown(window, { key: "1" });
    expect(calls).toEqual([{ hook: "removeFrame", boundTo: 101, arg: 7 }]);
  });

  it("stops ticking a phase once its member on the other frame is gone from that frame's cache", () => {
    state.assignment = { exposure_id: 101, state: "indexed", members: [] };
    state.sampleCall = {
      sample_id: 10, state: "indexed", phase: "Pn3m", call_phases: ["Pn3m"],
      members: [{ index_id: 7, exposure_id: 100, phase: "Pn3m" }],
    };
    state.frameAssignments = new Map([[100, { exposure_id: 100, state: "indexed", members: [] }]]);
    renderAt("/experiments/1/workspace/10");
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    expect(rows[0]).toHaveAttribute("aria-pressed", "false");
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

  it("shows the sample's labels under its heading, with suggestions dashed and confirmed by a click", () => {
    state.labels = [
      { id: 1, sample_id: 10, label: "37C", state: "confirmed", source: "user" },
      { id: 2, sample_id: 10, label: "buffer", state: "suggested", source: "ingest" },
      { id: 3, sample_id: 11, label: "retake", state: "confirmed", source: "user" },
    ];
    renderAt("/experiments/1/workspace/10");
    const row = screen.getByTestId("workspace-labels");
    const chips = screen.getByTestId("workspace-label-chips");
    expect([...chips.querySelectorAll("[data-label]")].map((c) => [c.getAttribute("data-label"), c.getAttribute("data-label-state")])).toEqual([
      ["37C", "confirmed"],
      ["buffer", "suggested"],
    ]);
    const buffer = chips.querySelector<HTMLElement>("[data-label='buffer']")!;
    expect(buffer).toHaveAttribute("data-label-standard", "true");
    expect(within(buffer).getByText("Buffer")).toBeInTheDocument();
    expect(within(chips.querySelector<HTMLElement>("[data-label='37C']")!).queryByText(", suggested")).toBeNull();
    expect(chips.querySelector("[data-label='37C']")).not.toHaveAttribute("data-label-standard");
    fireEvent.click(within(row).getByRole("button", { name: "Confirm label Buffer" }));
    fireEvent.click(within(row).getByRole("button", { name: "Remove label 37C" }));
    expect(calls).toEqual([
      { hook: "addLabel", boundTo: 1, arg: { sampleId: 10, label: "buffer" } },
      { hook: "removeLabel", boundTo: 1, arg: { sampleId: 10, label: "37C" } },
    ]);
  });

  it("the palette's label page lists experiment labels and toggles one on Enter, and creates a new one with a colour", () => {
    state.labels = [
      { id: 3, sample_id: 11, label: "retake", state: "confirmed", source: "user" },
      { id: 4, sample_id: 10, label: "37C", state: "confirmed", source: "user" },
    ];
    renderAt("/experiments/1/workspace/10");
    runFromPalette("label");
    const page = screen.getByRole("dialog", { name: "Label" });
    expect(within(page).getAllByRole("option").map((o) => o.getAttribute("aria-label"))).toEqual([
      "37C, On this sample. Press Enter to remove it.",
      "retake",
      "Buffer",
      "Calibrant",
      "Empty capillary",
    ]);
    const input = within(page).getByRole("combobox", { name: "Name" });
    fireEvent.change(input, { target: { value: "retake" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(calls).toEqual([{ hook: "addLabel", boundTo: 1, arg: { sampleId: 10, label: "retake" } }]);
    calls.length = 0;
    runFromPalette("label", "37c");
    expect(calls).toEqual([{ hook: "removeLabel", boundTo: 1, arg: { sampleId: 10, label: "37C" } }]);
    calls.length = 0;
    runFromPalette("label", "dehydrated");
    const form = screen.getByRole("dialog", { name: "New label" });
    expect(within(form).getByRole("textbox", { name: "Name" })).toHaveValue("dehydrated");
    fireEvent.keyDown(within(form).getByRole("textbox", { name: "Name" }), { key: "Tab" });
    expect(within(form).getByRole("option", { name: "Rose" })).toHaveAttribute("aria-selected", "true");
    fireEvent.keyDown(within(form).getByRole("combobox", { name: "Colour" }), { key: "ArrowDown" });
    fireEvent.keyDown(within(form).getByRole("combobox", { name: "Colour" }), { key: "Enter" });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(calls).toEqual([
      { hook: "setLabelColor", boundTo: 0, arg: { label: "dehydrated", color: "plum" } },
      { hook: "addLabel", boundTo: 1, arg: { sampleId: 10, label: "dehydrated" } },
    ]);
  });

  it("creating a label that already has a colour elsewhere starts on that colour and leaves it alone", () => {
    state.labelColors = [{ label: "dried", color: "plum" }];
    renderAt("/experiments/1/workspace/10");
    runFromPalette("label", "Dried");
    const form = screen.getByRole("dialog", { name: "New label" });
    expect(within(form).getByRole("combobox", { name: "Colour" })).toHaveTextContent("Plum");
    fireEvent.keyDown(within(form).getByRole("textbox", { name: "Name" }), { key: "Enter" });
    expect(calls).toEqual([{ hook: "addLabel", boundTo: 1, arg: { sampleId: 10, label: "Dried" } }]);
  });

  it("L opens the label picker without typing the L, and Enter confirms a suggestion", () => {
    state.labels = [{ id: 2, sample_id: 10, label: "buffer", state: "suggested", source: "ingest" }];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "l" });
    const input = within(screen.getByRole("dialog", { name: "Add a label" })).getByLabelText("Label name");
    expect(input).toHaveFocus();
    expect(input).toHaveValue("");
    fireEvent.change(input, { target: { value: "buf" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(calls).toEqual([{ hook: "addLabel", boundTo: 1, arg: { sampleId: 10, label: "Buffer" } }]);
  });

  it("the label picker lists labels added recently first, and the choice stays on this browser only", () => {
    state.labels = [
      { id: 1, sample_id: 11, label: "37C", state: "confirmed", source: "user" },
      { id: 2, sample_id: 12, label: "37C", state: "confirmed", source: "user" },
      { id: 3, sample_id: 12, label: "dried", state: "confirmed", source: "user" },
    ];
    renderAt("/experiments/1/workspace/10");
    const options = () =>
      [...screen.getByTestId("workspace-label-picker").querySelectorAll("[data-label-option='existing']")].map((b) => b.textContent);
    fireEvent.keyDown(window, { key: "l" });
    expect(options()).toEqual(["37C", "dried", "Buffer", "Calibrant", "Empty capillary"]);
    const input = within(screen.getByRole("dialog", { name: "Add a label" })).getByLabelText("Label name");
    fireEvent.change(input, { target: { value: "cal" } });
    fireEvent.keyDown(input, { key: "Enter" });
    fireEvent.keyDown(window, { key: "l" });
    expect(options()).toEqual(["Calibrant", "37C", "dried", "Buffer", "Empty capillary"]);
    expect(JSON.parse(localStorage.getItem("himalaya-ui:recent:label")!)).toEqual(["calibrant"]);
  });

  it("a recent label the experiment does not use is never listed or attached, and the palette's Label page shares the recent order", () => {
    localStorage.setItem("himalaya-ui:recent:label", JSON.stringify(["ghost", "dried"]));
    state.labels = [
      { id: 1, sample_id: 11, label: "37C", state: "confirmed", source: "user" },
      { id: 2, sample_id: 12, label: "37C", state: "confirmed", source: "user" },
      { id: 3, sample_id: 12, label: "dried", state: "confirmed", source: "user" },
    ];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "l" });
    const picker = screen.getByTestId("workspace-label-picker");
    expect([...picker.querySelectorAll("[data-label-option='existing']")].map((b) => b.textContent)).toEqual(["dried", "37C", "Buffer", "Calibrant", "Empty capillary"]);
    fireEvent.keyDown(within(picker).getByLabelText("Label name"), { key: "Escape" });
    runFromPalette("label");
    const names = screen.getAllByRole("option").map((o) => o.textContent ?? "");
    expect(names.findIndex((n) => n.startsWith("dried"))).toBeLessThan(names.findIndex((n) => n.startsWith("37C")));
    expect(names.some((n) => n.toLowerCase().includes("ghost"))).toBe(false);
    expect(calls).toEqual([]);
  });

  it("L does nothing before a user is set", () => {
    state.username = undefined;
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "L" });
    expect(screen.queryByRole("dialog", { name: "Add a label" })).toBeNull();
  });

  it("a new label needs the Create row, and Backspace in the empty picker removes the last confirmed label, never a suggestion", () => {
    state.labels = [
      { id: 1, sample_id: 10, label: "37C", state: "confirmed", source: "user" },
      { id: 2, sample_id: 10, label: "buffer", state: "suggested", source: "ingest" },
    ];
    renderAt("/experiments/1/workspace/10");
    fireEvent.click(screen.getByRole("button", { name: "+ Label L" }));
    const input = within(screen.getByRole("dialog", { name: "Add a label" })).getByLabelText("Label name");
    fireEvent.keyDown(input, { key: "Backspace" });
    fireEvent.change(input, { target: { value: "heated" } });
    expect(screen.getByRole("button", { name: "Create “heated”" })).toBeInTheDocument();
    fireEvent.keyDown(input, { key: "Enter" });
    expect(calls).toEqual([
      { hook: "removeLabel", boundTo: 1, arg: { sampleId: 10, label: "37C" } },
      { hook: "addLabel", boundTo: 1, arg: { sampleId: 10, label: "heated" } },
    ]);
  });

  it("the sample filter matches labels as well as names", () => {
    state.labels = [{ id: 1, sample_id: 13, label: "empty capillary", state: "confirmed", source: "user" }];
    renderAt("/experiments/1/workspace/10");
    fireEvent.change(screen.getByLabelText("Filter samples"), { target: { value: "capill" } });
    const list = screen.getByTestId("workspace-samples");
    expect(within(list).getAllByRole("button")).toHaveLength(1);
    expect(within(list).getByText("A4")).toBeInTheDocument();
  });

  it("keeps the label row mounted with no labels and no user, so the plot below never moves", () => {
    state.username = undefined;
    renderAt("/experiments/1/workspace/11");
    const row = screen.getByTestId("workspace-labels");
    expect(row.querySelector("[data-testid='workspace-label-chips'] [data-label]")).toBeNull();
  });

  it("colours each label from the shared map, system labels their system colour and other labels violet by default", () => {
    state.labels = [
      { id: 1, sample_id: 10, label: "37C", state: "confirmed", source: "user" },
      { id: 2, sample_id: 10, label: "buffer", state: "confirmed", source: "user" },
      { id: 3, sample_id: 10, label: "retake", state: "suggested", source: "ingest" },
    ];
    state.labelColors = [{ label: "RETAKE", color: "rose" }];
    renderAt("/experiments/1/workspace/10");
    const chips = screen.getByTestId("workspace-label-chips");
    expect([...chips.querySelectorAll("[data-label]")].map((c) => [c.getAttribute("data-label"), c.getAttribute("data-label-color")])).toEqual([
      ["37C", "violet"],
      ["buffer", "blue"],
      ["retake", "rose"],
    ]);
  });

  it("right-clicking a label opens its colour menu, and a pick sets the colour for that name", () => {
    state.labels = [{ id: 1, sample_id: 10, label: "37C", state: "confirmed", source: "user" }];
    renderAt("/experiments/1/workspace/10");
    fireEvent.contextMenu(screen.getByTestId("workspace-label-chips").querySelector("[data-label='37C']")!);
    const menu = screen.getByRole("dialog", { name: "Colour of label 37C" });
    expect(within(menu).getByRole("radio", { name: "Violet" })).toHaveAttribute("aria-checked", "true");
    expect(within(menu).queryByRole("button", { name: "Default colour" })).toBeNull();
    fireEvent.click(within(menu).getByRole("radio", { name: "Teal" }));
    expect(screen.queryByRole("dialog", { name: "Colour of label 37C" })).toBeNull();
    expect(calls).toEqual([{ hook: "setLabelColor", boundTo: 0, arg: { label: "37C", color: "teal" } }]);
  });

  it("the palette's label colour page recolours any label of the experiment through a colour form", () => {
    state.labels = [
      { id: 2, sample_id: 10, label: "buffer", state: "confirmed", source: "user" },
      { id: 3, sample_id: 11, label: "retake", state: "confirmed", source: "user" },
    ];
    state.labelColors = [{ label: "buffer", color: "stone" }];
    renderAt("/experiments/1/workspace/10");
    runFromPalette("label colour", "retake");
    const form = screen.getByRole("dialog", { name: "retake" });
    expect(within(form).getByRole("option", { name: "Violet" })).toHaveAttribute("aria-selected", "true");
    const colour = within(form).getByRole("combobox", { name: "Colour" });
    expect(colour).toHaveFocus();
    fireEvent.click(within(form).getByRole("option", { name: "Green" }));
    fireEvent.keyDown(colour, { key: "Enter" });
    expect(calls).toEqual([{ hook: "setLabelColor", boundTo: 0, arg: { label: "retake", color: "green" } }]);
  });

  it("system labels keep their own colours: no colour menu, and no row on the palette's colour page", () => {
    state.labels = [
      { id: 1, sample_id: 10, label: "buffer", state: "confirmed", source: "user" },
      { id: 2, sample_id: 10, label: "Calibrant", state: "confirmed", source: "user" },
      { id: 3, sample_id: 10, label: "retake", state: "confirmed", source: "user" },
    ];
    state.labelColors = [{ label: "buffer", color: "rose" }];
    renderAt("/experiments/1/workspace/10");
    const chips = screen.getByTestId("workspace-label-chips");
    expect(chips.querySelector("[data-label='buffer']")).toHaveAttribute("data-label-color", "blue");
    expect(chips.querySelector("[data-label='Calibrant']")).toHaveAttribute("data-label-color", "ochre");
    fireEvent.contextMenu(chips.querySelector("[data-label='buffer']")!);
    expect(screen.queryByRole("dialog", { name: /^Colour of label/ })).toBeNull();
    runFromPalette("label colour");
    const page = screen.getByRole("dialog", { name: "Label colour" });
    expect(within(page).getAllByRole("option").map((o) => o.getAttribute("aria-label"))).toEqual(["retake, 1 sample"]);
  });

  it("the palette's Settings page renames the experiment and sets its description", () => {
    renderAt("/experiments/1/workspace/10");
    runFromPalette("settings", "name");
    const name = within(screen.getByRole("dialog", { name: "Name" })).getByRole("textbox", { name: "Name" });
    expect(name).toHaveValue("BL-19");
    fireEvent.change(name, { target: { value: "BL-20" } });
    fireEvent.keyDown(name, { key: "Enter" });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(calls).toEqual([{ hook: "op:patchExperiment", boundTo: 0, arg: [{ name: "BL-20" }] }]);
    calls.length = 0;
    runFromPalette("settings", "description");
    const description = within(screen.getByRole("dialog", { name: "Description" })).getByRole("textbox", { name: "Description" });
    fireEvent.change(description, { target: { value: "Lipid screen" } });
    fireEvent.keyDown(description, { key: "Enter" });
    expect(calls).toEqual([{ hook: "op:patchExperiment", boundTo: 0, arg: [{ description: "Lipid screen" }] }]);
  });

  it("offers no Settings before a user is set", () => {
    state.username = undefined;
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    expect(document.querySelector("[data-palette-item='settings']")).toBeNull();
  });

  it("offers no colour menu before a user is set", () => {
    state.username = undefined;
    state.labels = [{ id: 1, sample_id: 10, label: "37C", state: "confirmed", source: "user" }];
    renderAt("/experiments/1/workspace/10");
    fireEvent.contextMenu(screen.getByTestId("workspace-label-chips").querySelector("[data-label='37C']")!);
    expect(screen.queryByRole("dialog", { name: "Colour of label 37C" })).toBeNull();
  });

  it("offers no label editing before a user is set", () => {
    state.username = undefined;
    state.labels = [{ id: 1, sample_id: 10, label: "37C", state: "confirmed", source: "user" }];
    renderAt("/experiments/1/workspace/10");
    const row = screen.getByTestId("workspace-labels");
    expect(within(row).queryAllByRole("button")).toHaveLength(0);
    expect(row.querySelector("[data-label='37C']")).not.toBeNull();
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
    state.exposures = [exposure(100), exposure(101, { status: "rejected" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(currentFrame()).toBe("Image 2, dropped");
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
    runFromPalette("plot linear");
    expect(currentView()).toBe("Linear");
    expect(ticks()).not.toBe(logTicks);
    runFromPalette("plot log-log");
    expect(currentView()).toBe("Log-log");
    runFromPalette("plot q²I");
    expect(currentView()).toBe("q²I (Kratky)");
    expect(screen.getByTestId("workspace-trace").textContent).toContain("q²·I");
    runFromPalette("plot view", "log");
    expect(currentView()).toBe("Log");
    expect(ticks()).toBe(logTicks);
  });

  it("V steps forward through the plot views and Shift+V steps back, wrapping", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    const pressed = (): string | null =>
      within(screen.getByRole("group", { name: "Plot view" })).getByRole("button", { pressed: true }).textContent;
    expect(pressed()).toBe("Log");
    fireEvent.keyDown(window, { key: "v" });
    expect(pressed()).toBe("Linear");
    fireEvent.keyDown(window, { key: "v" });
    fireEvent.keyDown(window, { key: "v" });
    expect(pressed()).toBe("q²I");
    expect(screen.getByTestId("workspace-trace").textContent).toContain("q²·I");
    fireEvent.keyDown(window, { key: "v" });
    expect(pressed()).toBe("Log");
    fireEvent.keyDown(window, { key: "V", shiftKey: true });
    expect(pressed()).toBe("q²I");
    expect(currentView()).toBe("q²I (Kratky)");
  });

  it("V does nothing while typing in the sample filter", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    const filter = document.getElementById("workspace-filter")!;
    filter.focus();
    fireEvent.keyDown(filter, { key: "v" });
    expect(currentView()).toBe("Log");
  });

  it("the plot view buttons switch the view directly and mark the current one", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    const group = screen.getByRole("group", { name: "Plot view" });
    fireEvent.click(within(group).getByRole("button", { name: "Log-log plot" }));
    expect(currentView()).toBe("Log-log");
    fireEvent.click(within(group).getByRole("button", { name: "q²I (Kratky) plot" }));
    expect(currentView()).toBe("q²I (Kratky)");
    expect(within(group).getByRole("button", { name: "q²I (Kratky) plot" })).toHaveAttribute("aria-pressed", "true");
    expect(within(group).getByRole("button", { name: "Log plot" })).toHaveAttribute("aria-pressed", "false");
    expect(fireEvent.mouseDown(within(group).getByRole("button", { name: "Linear plot" }))).toBe(false);
  });

  it("the difference replaces the trace and follows the plot view", () => {
    state.trace = makeTrace(1);
    const flat = makeTrace(1);
    state.compareTraces = new Map([[11, { ...flat, I: flat.q.map(() => 1) }]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    fireEvent.click(screen.getByRole("button", { name: "A2" }), { ctrlKey: true });
    const strip = (): string => screen.getByTestId("workspace-trace").querySelector("[data-role='difference-line'] path")!.getAttribute("d")!;
    const log = strip();
    runFromPalette("plot log-log");
    const loglog = strip();
    expect(loglog).not.toBe(log);
    runFromPalette("plot q²I");
    expect(strip()).not.toBe(log);
    expect(strip()).not.toBe(loglog);
  });

  it("C opens a picker, and Enter adds the sample scaled to the current trace without subtracting it", () => {
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
    expect(screen.getByRole("button", { name: "A2" })).not.toHaveAttribute("data-subtracted");
    expect(screen.queryByTestId("workspace-subtraction-toggle")).toBeNull();
    expect(screen.getByTestId("workspace-compares").querySelector("[data-compare='11']")?.getAttribute("d")).toMatch(/^M/);
  });

  it("the compare picker lists buffers first, marks which are only suggested, and picks nothing by itself", () => {
    state.labels = [
      { id: 1, sample_id: 13, label: "Buffer", state: "suggested", source: "ingest" },
      { id: 2, sample_id: 12, label: "buffer", state: "confirmed", source: "user" },
    ];
    renderAt("/experiments/1/workspace/10");
    expect(screen.queryByRole("button", { name: /^Minus / })).toBeNull();
    fireEvent.keyDown(window, { key: "c" });
    const rows = within(screen.getByRole("dialog", { name: "Compare with a sample" })).getAllByRole("button");
    expect(rows.map((r) => r.textContent)).toEqual(["A3Buffer", "A4Suggested buffer", "A2"]);
    expect(rows.map((r) => r.getAttribute("aria-current"))).toEqual(["true", null, null]);
    expect(screen.queryByRole("button", { name: /^Remove/ })).toBeNull();
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
    const monitors = (transmitted_counts: number) => ({ transmitted_counts, exposure_time: 15, counting_time: 15.8 });
    state.trace = { ...makeTrace(1), monitors: monitors(6000) };
    state.compareTraces = new Map([[11, { ...makeTrace(2), monitors: monitors(8000) }]]);
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

  it("Ctrl-clicking the subtracted chip stops subtracting, the palette hides comparisons, and × removes one", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    fireEvent.click(screen.getByRole("button", { name: "A2" }), { ctrlKey: true });
    expect(screen.getByTestId("workspace-subtraction-toggle")).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: "A2" }), { metaKey: true });
    expect(screen.queryByTestId("workspace-subtraction-toggle")).toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "A2" }), { metaKey: true });
    fireEvent.click(screen.getByRole("button", { name: "Raw" }));
    runFromPalette("Hide comparisons");
    expect(screen.queryByTestId("workspace-subtraction-toggle")).toBeNull();
    expect(screen.getByTestId("workspace-compares").children).toHaveLength(0);
    expect(screen.getByRole("button", { name: "A2, hidden" })).toBeInTheDocument();
    runFromPalette("Show comparisons");
    fireEvent.click(screen.getByRole("button", { name: "Remove A2" }));
    expect(screen.queryByTestId("workspace-subtraction-toggle")).toBeNull();
    expect(screen.getByTestId("workspace-compares").children).toHaveLength(0);
  });

  it("a plain click on a chip hides or shows just that trace and keeps the subtraction", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)], [12, makeTrace(3)]]);
    renderAt("/experiments/1/workspace/10");
    for (const id of [11, 12]) {
      fireEvent.keyDown(window, { key: "c" });
      fireEvent.click(screen.getByRole("dialog").querySelector(`[data-sample-id='${id}']`)!);
    }
    const lines = (): string[] => Array.from(screen.getByTestId("workspace-compares").querySelectorAll("[data-compare]")).map((n) => n.getAttribute("data-compare")!);
    expect(lines()).toEqual(["11", "12"]);
    fireEvent.click(screen.getByRole("button", { name: "A2" }), { ctrlKey: true });
    fireEvent.click(screen.getByRole("button", { name: "Raw" }));
    fireEvent.click(screen.getByRole("button", { name: "A2" }));
    expect(lines()).toEqual(["12"]);
    expect(screen.getByRole("button", { name: "A2, hidden" })).toHaveAttribute("aria-pressed", "false");
    expect(screen.getByRole("button", { name: "A2, hidden" })).toHaveAttribute("data-subtracted");
    expect(screen.getByRole("button", { name: "Raw" })).toHaveAttribute("aria-pressed", "true");
    fireEvent.click(screen.getByRole("button", { name: "A2, hidden" }));
    expect(lines()).toEqual(["11", "12"]);
  });

  it("Ctrl/Cmd-clicking a sidebar row toggles it as a comparison without leaving the sample", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    const row = (): HTMLElement => screen.getByTestId("workspace-samples").querySelector<HTMLElement>("[data-sample-id='11']")!;
    fireEvent.click(row(), { metaKey: true });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/10");
    expect(screen.getByRole("button", { name: "Remove A2" })).toBeInTheDocument();
    expect(row()).toHaveAttribute("data-compared");
    fireEvent.click(row(), { ctrlKey: true });
    expect(screen.queryByRole("button", { name: "Remove A2" })).toBeNull();
    expect(row()).not.toHaveAttribute("data-compared");
  });

  it("Ctrl/Cmd-clicking the open sample's row never makes it a comparison, even after moving on", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[10, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.click(screen.getByTestId("workspace-samples").querySelector<HTMLElement>("[data-sample-id='10']")!, { ctrlKey: true });
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
    expect(screen.queryByRole("button", { name: "Remove A1" })).toBeNull();
  });

  it("Ctrl/Cmd-clicking a row of the same sample compares it like any other sample", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[12, makeTrace(3)]]);
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    expect(screen.queryByRole("button", { name: "Remove A3" })).toBeNull();
    fireEvent.click(screen.getByTestId("workspace-samples").querySelector<HTMLElement>("[data-sample-id='12']")!, { metaKey: true });
    expect(screen.getByRole("button", { name: "Remove A3" })).toBeInTheDocument();
  });

  it("a hidden comparison shows again when it is removed and added back", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    const lines = (): number => screen.getByTestId("workspace-compares").querySelectorAll("[data-compare='11']").length;
    runFromPalette("compare with a sample", "a2");
    fireEvent.click(screen.getByRole("button", { name: "A2" }));
    expect(lines()).toBe(0);
    fireEvent.click(screen.getByRole("button", { name: "Remove A2" }));
    runFromPalette("compare with a sample", "a2");
    expect(lines()).toBe(1);
    fireEvent.click(screen.getByRole("button", { name: "A2" }));
    const row = screen.getByTestId("workspace-samples").querySelector<HTMLElement>("[data-sample-id='11']")!;
    fireEvent.click(row, { ctrlKey: true });
    fireEvent.click(row, { ctrlKey: true });
    expect(lines()).toBe(1);
  });

  it("clicking a chip while comparisons are hidden shows them all again", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)], [12, makeTrace(3)]]);
    renderAt("/experiments/1/workspace/10");
    runFromPalette("compare with a sample", "a2");
    runFromPalette("compare with a sample", "a3");
    runFromPalette("Hide comparisons");
    expect(screen.getByTestId("workspace-compares").children).toHaveLength(0);
    fireEvent.click(screen.getByRole("button", { name: "A2, hidden" }));
    const ids = Array.from(screen.getByTestId("workspace-compares").querySelectorAll("[data-compare]")).map((n) => n.getAttribute("data-compare"));
    expect(ids).toEqual(["11", "12"]);
  });

  it("the palette subtracts a comparison and stops, so keyboard users can pick the background", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)], [12, makeTrace(3)]]);
    renderAt("/experiments/1/workspace/10");
    runFromPalette("compare with a sample", "a2");
    runFromPalette("compare with a sample", "a3");
    expect(screen.queryByTestId("workspace-subtraction-toggle")).toBeNull();
    runFromPalette("subtract a comparison", "a2");
    expect(screen.getByRole("button", { name: "Minus A2" })).toHaveAttribute("aria-pressed", "true");
    runFromPalette("subtract a comparison", "a3");
    expect(screen.getByRole("button", { name: "Minus A3" })).toBeInTheDocument();
    runFromPalette("subtract a comparison", "a3");
    expect(screen.queryByTestId("workspace-subtraction-toggle")).toBeNull();
  });

  it("Ctrl-click subtracts in the main plot, Raw brings the trace and comparisons back, and a dimmed chip goes back to Raw", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)], [12, makeTrace(3)]]);
    renderAt("/experiments/1/workspace/10");
    runFromPalette("compare with a sample", "a2");
    runFromPalette("compare with a sample", "a3");
    const plot = (): HTMLElement => screen.getByTestId("workspace-trace");
    const lines = (): string[] => Array.from(screen.getByTestId("workspace-compares").querySelectorAll("[data-compare]")).map((n) => n.getAttribute("data-compare")!);
    expect(screen.queryByTestId("workspace-subtraction-toggle")).toBeNull();
    expect(plot().querySelector("[data-role='trace-line']")).not.toBeNull();
    fireEvent.click(screen.getByRole("button", { name: "A3" }), { ctrlKey: true });
    expect(screen.getByRole("button", { name: "Minus A3" })).toHaveAttribute("aria-pressed", "true");
    expect(plot().querySelector("[data-role='trace-line']")).toBeNull();
    expect(plot().querySelector("[data-role='difference-line']")).not.toBeNull();
    expect(lines()).toEqual([]);
    expect(screen.getByRole("button", { name: "A2, hidden" })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: "A3" })).toHaveAttribute("data-subtracted");
    const mark = (id: number): string | null | undefined => document.querySelector(`[data-compared][data-sample-id='${id}'] [role='img']`)?.getAttribute("aria-label");
    expect(mark(11)).toBe("Compared, hidden");
    expect(mark(12)).toBe("Compared");
    fireEvent.click(screen.getByRole("button", { name: "Raw" }));
    expect(plot().querySelector("[data-role='trace-line']")).not.toBeNull();
    expect(lines()).toEqual(["11", "12"]);
    expect(mark(11)).toBe("Compared");
    fireEvent.click(screen.getByRole("button", { name: "Minus A3" }));
    fireEvent.click(screen.getByRole("button", { name: "A2, hidden" }));
    expect(screen.getByRole("button", { name: "Raw" })).toHaveAttribute("aria-pressed", "true");
    expect(lines()).toEqual(["11", "12"]);
    expect(screen.getByRole("button", { name: "A3" })).toHaveAttribute("data-subtracted");
  });

  it("the palette switches between the raw trace and the difference, and marks the subtracted sample", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    runFromPalette("compare with a sample", "a2");
    runFromPalette("subtract a comparison", "a2");
    runFromPalette("show the raw trace");
    runFromPalette("show minus a2");
    expect(screen.getByRole("button", { name: "Minus A2" })).toHaveAttribute("aria-pressed", "true");
    runFromPalette("show the raw trace");
    expect(screen.getByRole("button", { name: "Raw" })).toHaveAttribute("aria-pressed", "true");
    runFromPalette("subtract a comparison");
    expect(screen.getByRole("option", { name: /A2/ })).toHaveTextContent("Subtracted, Enter stops");
  });

  it("Ctrl-clicking a chip while comparisons are hidden shows them and subtracts it", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    runFromPalette("compare with a sample", "a2");
    runFromPalette("Hide comparisons");
    fireEvent.click(screen.getByRole("button", { name: "A2, hidden" }), { ctrlKey: true });
    expect(screen.getByRole("button", { name: "Minus A2" })).toHaveAttribute("aria-pressed", "true");
    expect(screen.getByRole("button", { name: "A2" })).toHaveAttribute("data-subtracted");
  });

  it("the palette stops comparing with a sample, and offers that only while there is one", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    const offered = (): boolean => {
      fireEvent.keyDown(window, { key: "k", ctrlKey: true });
      fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "stop comparing" } });
      const found = screen.queryAllByRole("option").some((o) => o.getAttribute("data-palette-item") === "uncompare");
      fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Escape" });
      return found;
    };
    expect(offered()).toBe(false);
    runFromPalette("compare with a sample", "a2");
    expect(screen.getByRole("button", { name: "Remove A2" })).toBeInTheDocument();
    expect(offered()).toBe(true);
    runFromPalette("stop comparing", "a2");
    expect(screen.queryByRole("button", { name: "Remove A2" })).toBeNull();
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
    fireEvent.click(screen.getByRole("button", { name: "A3" }), { ctrlKey: true });
    fireEvent.click(screen.getByRole("button", { name: "Remove A2" }));
    pick(11);
    expect(screen.getByRole("button", { name: "A2, hidden" })).not.toHaveAttribute("data-subtracted");
    expect(screen.getByRole("button", { name: "A3" })).toHaveAttribute("data-subtracted");
    expect(screen.getByRole("button", { name: "Minus A3" })).toBeInTheDocument();
  });

  it("comparisons and the subtraction belong to the sample: moving on carries nothing, and coming back restores them", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "c" });
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A2"));
    fireEvent.click(screen.getByRole("button", { name: "A2" }), { ctrlKey: true });
    expect(calls.filter((c) => c.hook === "patchCompares").map((c) => [c.boundTo, c.arg])).toEqual([
      [10, { add: [11] }],
      [10, { subtraction: 11 }],
    ]);
    const overlay = (): Element | null => screen.getByTestId("workspace-compares").querySelector("[data-compare='11']");
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
    expect(screen.queryByRole("button", { name: "Remove A2" })).toBeNull();
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/12");
    expect(screen.queryByRole("button", { name: /^Remove/ })).toBeNull();
    expect(overlay()).toBeNull();
    fireEvent.keyDown(window, { key: "ArrowUp" });
    fireEvent.keyDown(window, { key: "ArrowUp" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/10");
    expect(screen.getByRole("button", { name: "A2" })).toHaveAttribute("data-subtracted");
    expect(screen.getByRole("button", { name: "Minus A2" })).toHaveAttribute("aria-pressed", "true");
  });

  it("a comparison someone else added arrives hidden, and yours start shown", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)], [12, makeTrace(3)]]);
    renderAt("/experiments/1/workspace/10");
    const lines = (): string[] => Array.from(screen.getByTestId("workspace-compares").querySelectorAll("[data-compare]")).map((n) => n.getAttribute("data-compare")!);
    act(() => setCompares({ compares: [{ sample_id: 10, ref_sample_id: 11, added_by: "ana" }], subtractions: [{ sample_id: 10, ref_sample_id: 11 }] }));
    expect(screen.getByRole("button", { name: "A2, hidden" })).toHaveAttribute("data-subtracted");
    expect(lines()).toEqual([]);
    expect(screen.queryByTestId("workspace-subtraction-toggle")).not.toBeNull();
    expect(screen.getByRole("button", { name: "Raw" })).toHaveAttribute("aria-pressed", "true");
    fireEvent.click(screen.getByRole("button", { name: "A2, hidden" }));
    expect(lines()).toEqual(["11"]);
    runFromPalette("compare with a sample", "a3");
    expect(lines()).toEqual(["11", "12"]);
  });

  it("a sample holds five comparisons, one per colour, and the picker, palette and sample list say so at the cap", () => {
    state.corpus = [...state.corpus, sample(14, "A5"), sample(15, "A6"), sample(16, "A7")];
    expect(COMPARE_COLORS).toHaveLength(MAX_COMPARES);
    act(() => setCompares({ compares: [11, 12, 13, 14, 15].map((r) => ({ sample_id: 10, ref_sample_id: r, added_by: "jc" })), subtractions: [] }));
    renderAt("/experiments/1/workspace/10");
    const colours = Array.from(screen.getByTestId("workspace-compare-bar").querySelectorAll<HTMLElement>("[data-compare-chip]")).map((n) => n.style.color);
    expect(new Set(colours).size).toBe(MAX_COMPARES);
    fireEvent.keyDown(window, { key: "c" });
    const picker = screen.getByRole("dialog", { name: "Compare with a sample" });
    expect(within(picker).queryAllByRole("button")).toHaveLength(0);
    expect(picker).toHaveTextContent(COMPARES_FULL);
    fireEvent.keyDown(within(picker).getByLabelText("Find a sample to compare"), { key: "Escape" });
    runFromPalette("compare with a sample");
    expect(screen.getByRole("dialog")).toHaveTextContent(COMPARES_FULL);
    fireEvent.keyDown(screen.getByRole("combobox"), { key: "Escape" });
    calls.length = 0;
    fireEvent.click(within(screen.getByTestId("workspace-samples")).getByRole("button", { name: /A7/ }), { ctrlKey: true });
    expect(calls.filter((c) => c.hook === "patchCompares")).toEqual([]);
    expect(toasts.at(-1)).toBe(COMPARES_FULL);
  });

  it("showing or hiding a comparison is kept per sample, even when two samples compare the same one", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[12, makeTrace(2)]]);
    act(() => setCompares({ compares: [{ sample_id: 10, ref_sample_id: 12, added_by: "jc" }, { sample_id: 11, ref_sample_id: 12, added_by: "jc" }], subtractions: [] }));
    renderAt("/experiments/1/workspace/10");
    fireEvent.click(screen.getByRole("button", { name: "A3" }));
    expect(screen.getByRole("button", { name: "A3, hidden" })).toBeInTheDocument();
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
    expect(screen.getByRole("button", { name: "A3" })).toHaveAttribute("aria-pressed", "true");
    fireEvent.keyDown(window, { key: "ArrowUp" });
    expect(screen.getByRole("button", { name: "A3, hidden" })).toBeInTheDocument();
  });

  it("Minus is kept per sample: another sample subtracting the same comparison opens on Raw", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[12, makeTrace(2)]]);
    act(() => setCompares({
      compares: [{ sample_id: 10, ref_sample_id: 12, added_by: "jc" }, { sample_id: 11, ref_sample_id: 12, added_by: "jc" }],
      subtractions: [{ sample_id: 11, ref_sample_id: 12 }],
    }));
    renderAt("/experiments/1/workspace/10");
    fireEvent.click(screen.getByRole("button", { name: "A3" }), { ctrlKey: true });
    expect(screen.getByRole("button", { name: "Minus A3" })).toHaveAttribute("aria-pressed", "true");
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/11");
    expect(screen.getByRole("button", { name: "Raw" })).toHaveAttribute("aria-pressed", "true");
  });

  it("when someone else changes the subtraction, Minus stays on and uses the new background", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)], [12, makeTrace(3)]]);
    act(() => setCompares({ compares: [11, 12].map((r) => ({ sample_id: 10, ref_sample_id: r, added_by: "jc" })), subtractions: [] }));
    renderAt("/experiments/1/workspace/10");
    fireEvent.click(screen.getByRole("button", { name: "A2" }), { ctrlKey: true });
    expect(screen.getByRole("button", { name: "Minus A2" })).toHaveAttribute("aria-pressed", "true");
    act(() => setCompares(applyCompareChange(compareStore.data, 10, { subtraction: 12 }, "ana")));
    expect(screen.getByRole("button", { name: "Minus A3" })).toHaveAttribute("aria-pressed", "true");
  });

  it("a hidden comparison someone else added fetches no trace until it is shown", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[11, makeTrace(2)]]);
    renderAt("/experiments/1/workspace/10");
    act(() => setCompares({ compares: [{ sample_id: 10, ref_sample_id: 11, added_by: "ana" }], subtractions: [] }));
    expect(compareFetches.at(-1)).toEqual([]);
    fireEvent.click(screen.getByRole("button", { name: "A2, hidden" }));
    expect(compareFetches.at(-1)).toEqual([11]);
  });

  it("the compare picker lists samples you compared recently first, without adding them anywhere", () => {
    state.trace = makeTrace(1);
    state.labels = [{ id: 1, sample_id: 12, label: "Buffer", state: "confirmed", source: "user" }];
    renderAt("/experiments/1/workspace/10");
    const rows = (): string[] => within(screen.getByRole("dialog", { name: "Compare with a sample" })).getAllByRole("button").map((r) => r.textContent ?? "");
    fireEvent.keyDown(window, { key: "c" });
    expect(rows()[0]).toBe("A3Buffer");
    fireEvent.click(within(screen.getByRole("dialog")).getByText("A4"));
    fireEvent.click(screen.getByRole("button", { name: "Remove A4" }));
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(screen.queryByRole("button", { name: /^Remove/ })).toBeNull();
    fireEvent.keyDown(window, { key: "c" });
    expect(rows()[0]).toBe("A4");
  });

  it("drops a kept zoom that falls outside the next sample's q range, and ⌘/Ctrl+0 resets the zoom", () => {
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

  it("typing a sample's exact name and Enter goes there rather than comparing with it", () => {
    renderAt("/experiments/1/workspace/10");
    runFromPalette("A3");
    expect(screen.getByTestId("where").textContent).toBe("/experiments/1/workspace/12");
    expect(calls).toEqual([]);
  });

  it("Ctrl+K closes the palette from inside a page", () => {
    renderAt("/experiments/1/workspace/10");
    runFromPalette("compare with a sample");
    const page = screen.getByRole("combobox", { name: "Sample" });
    fireEvent.keyDown(page, { key: "k", ctrlKey: true });
    expect(screen.queryByRole("dialog")).toBeNull();
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
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "plot view" } });
    fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Enter" });
    expect(screen.getByRole("option", { name: "q²I (Kratky), Current" })).toBeInTheDocument();
    fireEvent.keyDown(screen.getByRole("combobox", { name: "View" }), { key: "Escape" });
    expect(screen.getByRole("combobox", { name: "Command or sample" })).toHaveValue("plot view");
    fireEvent.keyDown(screen.getByRole("combobox"), { key: "Escape" });
    expect(screen.queryByRole("dialog")).toBeNull();
  });

  it("the palette's theme page switches the workspace to dark and back, remembers it, and can follow the system again", () => {
    const { unmount } = renderAt("/experiments/1/workspace/10");
    const root = screen.getByTestId("workspace");
    expect(root).toHaveAttribute("data-theme", "system");
    runFromPalette("theme", "dark");
    expect(root).toHaveAttribute("data-theme", "dark");
    unmount();
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace")).toHaveAttribute("data-theme", "dark");
    runFromPalette("theme");
    expect(screen.getByRole("option", { name: "Dark, Current" })).toBeInTheDocument();
    const input = screen.getByRole("combobox", { name: "Theme" });
    fireEvent.change(input, { target: { value: "light" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(screen.getByTestId("workspace")).toHaveAttribute("data-theme", "light");
    runFromPalette("theme", "system");
    expect(screen.getByTestId("workspace")).toHaveAttribute("data-theme", "system");
    expect(window.localStorage.getItem("himalaya-ui:theme")).toBeNull();
    runFromPalette("theme");
    expect(screen.getByRole("option", { name: "Follow the system, Current" })).toBeInTheDocument();
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
    expect(currentFrame()).toBe("Image 2");
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

  it("Tab and Shift+Tab step the custom index's symmetry with wrap, from the page or from inside the panel", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    const panel = screen.getByTestId("workspace-custom-index");
    const pressed = (): string | null => within(panel).getAllByRole("button", { pressed: true })[0]!.textContent;
    const symmetries = within(within(panel).getByRole("group", { name: "Symmetry" })).getAllByRole("button").map((b) => b.textContent);
    expect(pressed()).toBe("Pn3m");
    const at = symmetries.indexOf("Pn3m");
    fireEvent.keyDown(window, { key: "Tab" });
    expect(pressed()).toBe(symmetries[(at + 1) % symmetries.length]);
    fireEvent.keyDown(window, { key: "Tab", shiftKey: true });
    fireEvent.keyDown(window, { key: "Tab", shiftKey: true });
    expect(pressed()).toBe(symmetries[(at - 1 + symmetries.length) % symmetries.length]);
    const field = within(panel).getByRole("spinbutton");
    field.focus();
    fireEvent.keyDown(field, { key: "Tab" });
    expect(pressed()).toBe("Pn3m");
    expect(document.activeElement).toBe(field);
    expect(calls).toEqual([]);
  });

  it("Escape cancels a custom index without writing", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.queryByTestId("workspace-custom-index")).toBeNull();
    expect(calls).toEqual([]);
  });

  it("a custom index on an unselected frame lands on that frame without selecting it", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    fireEvent.keyDown(window, { key: "i" });
    fireEvent.click(within(screen.getByTestId("workspace-custom-index")).getByRole("button", { name: /^Add/ }));
    expect(calls).toHaveLength(1);
    expect(calls[0]).toMatchObject({ hook: "custom", boundTo: 100 });
  });

  it("does not open a custom index on a dropped frame or before a user is set", () => {
    state.username = undefined;
    const { unmount } = renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    expect(screen.queryByTestId("workspace-custom-index")).toBeNull();
    unmount();
    state.username = "jc";
    state.exposures = [exposure(100), exposure(101, { status: "rejected" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(currentFrame()).toBe("Image 2, dropped");
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

  it("the indexing plot waits for the trace before counting reflections", () => {
    state.trace = undefined;
    state.peaks = [peak(7, 0.1, "auto"), peak(8, 0.14, "auto")];
    state.indices = [
      {
        ...index(2, "Pn3m", 0.9),
        peaks: [
          { peak_id: 7, ratio_position: 1, residual: 0, q_observed: 0.1 },
          { peak_id: 8, ratio_position: 2, residual: 0, q_observed: 0.14 },
        ],
      },
    ];
    renderAt("/experiments/1/workspace/10");
    expect(screen.queryByTestId("workspace-indexing-fit")).toBeNull();
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
    expect(screen.getByRole("heading", { name: "Assignment" })).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: "Image" })).toBeInTheDocument();
  });

  it("folds Assignment and Image as an accordion on a short window, and opening a custom index from the folded Assignment header shows it in place of the list", () => {
    vi.stubGlobal("matchMedia", (query: string) => ({ matches: query === "(max-height: 860px)", media: query, addEventListener: () => {}, removeEventListener: () => {} }));
    try {
      renderAt("/experiments/1/workspace/10");
      const tools = screen.getByRole("complementary", { name: "Tools" });
      const call = within(tools).getByRole("button", { name: /Assignment/ });
      const exposure = within(tools).getByRole("button", { name: /^Image/ });
      expect(call).toHaveAttribute("aria-expanded", "true");
      expect(exposure).toHaveAttribute("aria-expanded", "false");
      expect(within(tools).queryByText("File")).toBeNull();
      fireEvent.click(exposure);
      expect(exposure).toHaveAttribute("aria-expanded", "true");
      expect(within(tools).getByRole("button", { name: /Assignment/ })).toHaveAttribute("aria-expanded", "false");
      expect(within(tools).getByText("File")).toBeInTheDocument();
      fireEvent.click(screen.getByRole("button", { name: /Custom index/ }));
      expect(screen.getByTestId("workspace-custom-index")).toBeInTheDocument();
      expect(within(tools).queryByRole("button", { name: /^Image/ })).toBeNull();
      expect(within(tools).getByRole("heading", { name: /^Image/ })).toBeInTheDocument();
      expect(within(tools).queryByText("File")).toBeNull();
    } finally {
      vi.unstubAllGlobals();
    }
  });

  it("on a short window the folded Frame header names the shown frame's sample for a bridged sample, and only then", () => {
    vi.stubGlobal("matchMedia", (query: string) => ({ matches: query === "(max-height: 860px)", media: query, addEventListener: () => {}, removeEventListener: () => {} }));
    try {
      state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
      state.corpus = state.corpus.map((s) => (s.id === 12 ? { ...s, assignment_state: "indexed" } : s));
      state.mateExposures = new Map([[12, [exposure(200, { sample_id: 12 })]]]);
      const { unmount } = renderAt("/experiments/1/workspace/10");
      const tools = () => screen.getByRole("complementary", { name: "Tools" });
      const frame = () => within(tools()).getByRole("button", { name: /^Image/ }).closest("section")!;
      expect(frame().querySelector("[data-section-folded]")).toHaveTextContent("A1");
      fireEvent.click(screen.getByRole("button", { name: "Image 1 of A3" }));
      expect(frame().querySelector("[data-section-folded]")).toHaveTextContent("A3");
      fireEvent.click(within(tools()).getByRole("button", { name: /^Image/ }));
      expect(frame().querySelector("[data-section-folded]")).toBeNull();
      unmount();
      state.bridges = [];
      renderAt("/experiments/1/workspace/10");
      expect(frame().querySelector("[data-section-folded]")).toBeNull();
    } finally {
      vi.unstubAllGlobals();
    }
  });

  it("puts the Assignment list first in the tools column, then the indexing plot, then the exposure", () => {
    renderAt("/experiments/1/workspace/10");
    const tools = screen.getByRole("complementary", { name: "Tools" });
    const headings = within(tools).getAllByRole("heading").map((h) => h.textContent);
    expect(headings).toEqual(["Assignment", "Indexing", "Image"]);
  });

  it("the Custom index button sits in the Assignment header and stays there while the custom index it opens replaces the list", () => {
    renderAt("/experiments/1/workspace/10");
    const header = () => screen.getByRole("heading", { name: "Assignment" }).closest("header")!;
    const button = within(header()).getByRole("button", { name: /Custom index/ });
    expect(within(screen.getByRole("heading", { name: "Indexing" }).closest("header")!).queryByRole("button", { name: /Custom index/ })).toBeNull();
    expect(button).toHaveTextContent("Custom index I");
    fireEvent.click(button);
    expect(screen.getByTestId("workspace-custom-index")).toBeInTheDocument();
    expect(within(header()).getByRole("button", { name: /Custom index/ })).toBe(button);
    expect(button).toHaveAttribute("aria-pressed", "true");
    expect(button).toHaveTextContent("Custom index Esc");
  });

  it("disables the Custom index button before a user is set", () => {
    state.username = undefined;
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByRole("button", { name: /Custom index/ })).toBeDisabled();
  });

  it("lists every exposure field, with a dash for the ones the frame lacks, so the section keeps its height", () => {
    state.exposures = [exposure(100), exposure(101)];
    renderAt("/experiments/1/workspace/10");
    const section = screen.getByRole("heading", { name: "Image" }).closest("section")!;
    const text = section.textContent ?? "";
    expect(text).toContain("KindSingle");
    expect(text).toContain("Exposure—");
    expect(text).toContain("I₀—");
    expect(text).toContain("I signal—");
  });

  it("lists the shown exposure's time, I₀ and I signal under the .prp names", () => {
    state.exposures = [exposure(100), exposure(101, { exposure_time: 10, incident_monitor: 1.5007, transmitted_counts: 7798.4 })];
    renderAt("/experiments/1/workspace/10");
    const section = screen.getByRole("heading", { name: "Image" }).closest("section")!;
    const text = section.textContent ?? "";
    expect(text).toContain("Exposure10.0 s");
    expect(text).toContain("I₀1.501");
    expect(text).toContain("I signal7,798");
  });

  it("lists the shown exposure's acquisition time, flight tube, energy and detector from the .prp", () => {
    state.exposures = [exposure(100), exposure(101, {
      timestamp: "2026-04-26T18:20:27", pipe_length_m: 0.7, beam_energy_ev: 11000.2, detector: "Area detector",
    })];
    renderAt("/experiments/1/workspace/10");
    const section = screen.getByRole("heading", { name: "Image" }).closest("section")!;
    const text = section.textContent ?? "";
    expect(text).toContain("Acquired2026-04-26 18:20:27");
    expect(text).toContain("Flight tube700 mm");
    expect(text).toContain("Energy11.000 keV");
    expect(text).toContain("DetectorArea detector");
  });

  it("keeps the new .prp rows with a dash when the frame has no .prp", () => {
    renderAt("/experiments/1/workspace/10");
    const text = screen.getByRole("heading", { name: "Image" }).closest("section")!.textContent ?? "";
    for (const k of ["Acquired", "Flight tube", "Energy", "Detector"]) expect(text).toContain(`${k}—`);
  });

  it("the indexing plot follows the custom index while it is open", () => {
    state.trace = makeTrace(1);
    const a = (2 * Math.PI * Math.SQRT2) / 0.1;
    state.peaks = [peak(7, 0.1, "auto"), peak(8, 0.1 * Math.sqrt(1.5), "auto")];
    state.indices = [{ ...index(2, "Pn3m", 0.9), lattice_d: a }];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    expect(screen.getByTestId("workspace-custom-landed")).toHaveTextContent("2 of 9");
    expect(screen.getByTestId("workspace-indexing-fit")).toHaveTextContent("2/9 on peaks");
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
    fireEvent.click(screen.getByRole("button", { name: "Image 1" }));
    expect(currentFrame()).toBe("Image 1");
    expect(screen.queryByTestId("workspace-custom-index")).toBeNull();
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
  it("Space while hovering a candidate draws its half-rings on the enlarged image only, and keeps them after the mouse leaves the list", () => {
    state.setups = [{ id: 1, name: "setup1", folder: "setup1", included: true,
      beam_center_x: 421.409, beam_center_y: 836.946, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8095, q_units: "A-1" }];
    state.corpus = state.corpus.map((s) => (s.id === 10 ? { ...s, setup_id: 1 } : s));
    state.exposures = [exposure(100, { image_path: "/d/100.tif", image_version: "v" }), exposure(101, { image_path: "/d/101.tif", image_version: "v" })];
    state.detectorSize = { width: 981, height: 1043 };
    renderAt("/experiments/1/workspace/10");
    const ringSets = () => [...document.querySelectorAll('[data-testid="detector-rings"]')];
    expect(new Set(state.sizeRequests)).toEqual(new Set(["/api/exposures/100/image?thumb=1&v=v", "/api/exposures/101/image?thumb=1&v=v"]));
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    fireEvent.mouseEnter(rows[1]!);
    expect(ringSets()).toHaveLength(0);
    fireEvent.keyDown(window, { key: " " });
    fireEvent.mouseLeave(screen.getByTestId("workspace-candidates"));
    const big = screen.getByTestId("workspace-detector-big");
    expect(ringSets()).toHaveLength(1);
    const arcs = [...big.querySelectorAll('[data-role="det-ring"]')];
    expect(arcs).toHaveLength(2);
    const [, x0, y0, , r, , , , sweep, x1, y1] = arcs[0]!.getAttribute("d")!.split(" ");
    expect(Number(r)).toBeCloseTo(230.7, 0);
    expect(sweep).toBe("0");
    expect([x1, Number(y1) - Number(y0)]).toEqual([x0, 2 * Number(r)]);
    expect(arcs[0]!.getAttribute("stroke")).toBe(phaseColor("Im3m"));
    expect(within(screen.getByTestId("workspace-exposures")).queryByTestId("detector-rings")).toBeNull();
  });

  it("Space with no candidate hovered draws the ticked phases' rings on the enlarged image", () => {
    state.setups = [{ id: 1, name: "setup1", folder: "setup1", included: true,
      beam_center_x: 421.409, beam_center_y: 836.946, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8095, q_units: "A-1" }];
    state.corpus = state.corpus.map((s) => (s.id === 10 ? { ...s, setup_id: 1 } : s));
    state.exposures = [exposure(100, { image_path: "/d/100.tif", image_version: "v" }), exposure(101, { image_path: "/d/101.tif", image_version: "v" })];
    state.detectorSize = { width: 981, height: 1043 };
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: " " });
    const strokes = [...screen.getByTestId("workspace-detector-big").querySelectorAll('[data-role="det-ring"]')].map((a) => a.getAttribute("stroke"));
    expect(strokes.length).toBeGreaterThan(0);
    expect(new Set(strokes)).toEqual(new Set([phaseColor("Pn3m")]));
  });

  it("draws a same sample's rings from another setup in that setup's geometry but in this sample's q units", () => {
    const geometry = { beam_center_x: 421.409, beam_center_y: 836.946, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8095 };
    state.setups = [{ id: 1, name: "Short pipe", folder: "s1", included: true, ...geometry, q_units: "A-1" },
      { id: 2, name: "Long pipe", folder: "s2", included: true, ...geometry, q_units: "nm-1" }];
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    state.corpus = state.corpus.map((s) => (s.id === 10 ? { ...s, setup_id: 1 } : s.id === 12 ? { ...s, setup_id: 2, assignment_state: "indexed" } : s));
    state.exposures = [exposure(100, { image_path: "/d/100.tif", image_version: "v" })];
    state.mateExposures = new Map([[12, [exposure(200, { sample_id: 12, image_path: "/d/200.tif", image_version: "v" })]]]);
    state.detectorSize = { width: 981, height: 1043 };
    renderAt("/experiments/1/workspace/10");
    fireEvent.mouseEnter(within(screen.getByTestId("workspace-candidates")).getAllByRole("button")[1]!);
    fireEvent.keyDown(window, { key: " " });
    const radius = () => Number(screen.getByTestId("workspace-detector-big").querySelector('[data-role="det-ring"]')!.getAttribute("d")!.split(" ")[4]);
    const first = radius();
    expect(first).toBeCloseTo(230.7, 0);
    fireEvent.keyDown(screen.getByTestId("workspace-detector-big"), { key: "ArrowRight" });
    expect(screen.getByTestId("workspace-detector-caption").textContent).toMatch(/200|Image 2/);
    expect(radius()).toBeCloseTo(first, 6);
  });

  it("Space opens the enlarged image while a custom index is open, with the draft's rings, and Esc there keeps the draft", () => {
    state.setups = [{ id: 1, name: "setup1", folder: "setup1", included: true,
      beam_center_x: 421.409, beam_center_y: 836.946, pixel_size_um: 172, energy_kev: 9, flight_path_m: 1.8095, q_units: "A-1" }];
    state.corpus = state.corpus.map((s) => (s.id === 10 ? { ...s, setup_id: 1 } : s));
    state.exposures = [exposure(100, { image_path: "/d/100.tif", image_version: "v" }), exposure(101, { image_path: "/d/101.tif", image_version: "v" })];
    state.detectorSize = { width: 981, height: 1043 };
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "i" });
    fireEvent.keyDown(window, { key: " " });
    const big = screen.getByTestId("workspace-detector-big");
    expect(big.querySelectorAll('[data-testid="detector-rings"] [data-role="det-ring"]').length).toBeGreaterThan(0);
    fireEvent.keyDown(big, { key: "Escape" });
    expect(screen.queryByTestId("workspace-detector-big")).toBeNull();
    expect(screen.getByTestId("workspace-custom-index")).toBeInTheDocument();
  });

  it("draws no rings for a setup that lacks its geometry", () => {
    state.setups = [{ id: 1, name: "setup1", folder: "setup1", included: true, beam_center_x: null }];
    state.corpus = state.corpus.map((s) => (s.id === 10 ? { ...s, setup_id: 1 } : s));
    state.exposures = [exposure(100, { image_path: "/d/100.tif", image_version: "v" })];
    state.detectorSize = { width: 981, height: 1043 };
    renderAt("/experiments/1/workspace/10");
    fireEvent.mouseEnter(within(screen.getByTestId("workspace-candidates")).getAllByRole("button")[1]!);
    expect(screen.getByTestId("workspace-exposures").querySelectorAll('[data-testid="detector-rings"]')).toHaveLength(0);
    expect(state.sizeRequests).toEqual([]);
  });

  it("the plot keeps the called phase's reflections faded, and draws a hovered candidate's at full strength instead", () => {
    state.trace = makeTrace(1);
    state.peaks = [peak(7, 0.1, "auto")];
    renderAt("/experiments/1/workspace/10");
    const combs = (): Element[] => Array.from(screen.getByTestId("workspace-combs").querySelectorAll("[data-comb]"));
    const ids = (): Set<string> => new Set(combs().map((l) => l.getAttribute("data-comb") ?? ""));
    const faded = (): boolean => combs().every((l) => l.getAttribute("stroke-opacity") !== null);
    expect(combs().length).toBeGreaterThan(0);
    expect(ids()).toEqual(new Set(["2"]));
    expect(faded()).toBe(true);
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    fireEvent.mouseEnter(rows[1]!);
    expect(ids()).toEqual(new Set(["1"]));
    expect(combs().some((l) => l.getAttribute("stroke-opacity") !== null)).toBe(false);
    fireEvent.mouseLeave(screen.getByTestId("workspace-candidates"));
    expect(ids()).toEqual(new Set(["2"]));
    expect(faded()).toBe(true);
  });

  it("the plot draws no reflections when nothing is called or hovered", () => {
    state.trace = makeTrace(1);
    state.peaks = [peak(7, 0.1, "auto")];
    state.assignment = { exposure_id: 101, state: "null", members: [] };
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-combs").querySelectorAll("[data-comb]")).toHaveLength(0);
  });
  it("hovering a candidate, a custom row or a Bonnet group header colours the peaks it claims on the trace", () => {
    state.trace = makeTrace(1);
    state.peaks = [peak(7, 0.05, "auto"), peak(8, 0.06, "auto"), peak(9, 0.07, "manual"), peak(10, 0.2, "auto"), peak(11, 0.08, "auto")];
    const claims = (pairs: [number, number][]) => pairs.map(([id, q], i) => ({ peak_id: id, ratio_position: i + 1, residual: 0, q_observed: q }));
    state.indices = [
      { ...index(1, "Im3m", 0.9), lattice_d: 127.9, peaks: claims([[7, 0.05], [9, 0.07]]) },
      { ...index(2, "Pn3m", 0.6), lattice_d: 100, peaks: claims([[8, 0.06], [11, 0.08]]) },
      { ...index(3, "Ia3d", 0), kind: "speculative", score: null, lattice_d: 300, peaks: claims([[10, 0.2]]) },
    ];
    state.assignment = { exposure_id: 101, state: "null", members: [] };
    renderAt("/experiments/1/workspace/10");
    const glyph = (key: string): Element => screen.getByTestId("workspace-trace").querySelector(`[data-peak-key="${key}"]`)!;
    const dimmed = (key: string): boolean => glyph(key).getAttribute("opacity") !== null;
    const list = screen.getByTestId("workspace-candidates");
    const row = (name: RegExp): HTMLElement => within(list).getAllByRole("button").find((b) => name.test(b.textContent ?? ""))!;
    expect(["auto:7", "auto:8", "manual:9", "auto:10"].some(dimmed)).toBe(false);

    fireEvent.mouseEnter(row(/Im3m/));
    expect(glyph("auto:7").getAttribute("fill")).toBe(phaseColor("Im3m"));
    expect(glyph("manual:9").getAttribute("fill")).toBe(phaseColor("Im3m"));
    expect(dimmed("auto:8")).toBe(true);
    expect(dimmed("auto:10")).toBe(true);

    fireEvent.mouseEnter(row(/Custom/));
    expect(glyph("auto:10").getAttribute("fill")).toBe(phaseColor("Ia3d"));
    expect(dimmed("auto:7")).toBe(true);

    fireEvent.mouseEnter(screen.getByTestId("workspace-bonnet-group"));
    expect(glyph("auto:7").getAttribute("fill")).toBe(phaseColor("Im3m"));
    expect(glyph("auto:8").getAttribute("fill")).toBe(phaseColor("Pn3m"));
    expect(dimmed("auto:10")).toBe(true);

    fireEvent.mouseLeave(list);
    expect(["auto:7", "auto:8", "manual:9", "auto:10"].some(dimmed)).toBe(false);
    expect(glyph("auto:7").getAttribute("fill")).toBe("var(--color-ink)");
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
    expect(currentFrame()).toBe("Image 1");
    unmount();
    renderAt("/experiments/1/workspace/10?exposure=999");
    expect(currentFrame()).toBe("Image 2");
  });

  it("Tab and Shift+Tab cycle only the phase choices, wrapping", () => {
    renderAt("/experiments/1/workspace/10");
    const rows = within(screen.getByTestId("workspace-candidates")).getAllByRole("button");
    const ff = screen.getByTestId("workspace-form-factor");
    const choices = [...rows, ff];
    expect(fireEvent.keyDown(document.body, { key: "Tab" })).toBe(false);
    expect(document.activeElement).toBe(choices[0]);
    for (let i = 1; i < choices.length; i++) {
      fireEvent.keyDown(document.activeElement!, { key: "Tab" });
      expect(document.activeElement).toBe(choices[i]);
    }
    fireEvent.keyDown(document.activeElement!, { key: "Tab" });
    expect(document.activeElement).toBe(choices[0]);
    fireEvent.keyDown(document.activeElement!, { key: "Tab", shiftKey: true });
    expect(document.activeElement).toBe(ff);
    fireEvent.keyDown(screen.getByRole("button", { name: "Image 1" }), { key: "Tab", shiftKey: true });
    expect(document.activeElement).toBe(ff);
  });

  it("leaves Enter on a focused choice to the button, and Space keeps opening the detector view", () => {
    state.exposures = [exposure(100), exposure(101, { image_path: "/x.tif", image_version: "7" })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(document.body, { key: "Tab" });
    expect(fireEvent.keyDown(document.activeElement!, { key: "Enter" })).toBe(true);
    fireEvent.keyDown(document.activeElement!, { key: " " });
    expect(screen.getByRole("dialog", { name: "Detector image" })).toBeInTheDocument();
  });

  it("a stray Esc outside the phase choices leaves Tab on them", () => {
    renderAt("/experiments/1/workspace/10");
    expect(fireEvent.keyDown(document.body, { key: "Escape" })).toBe(true);
    expect(fireEvent.keyDown(document.body, { key: "Tab" })).toBe(false);
    expect(document.activeElement?.closest("[data-call-choices]")).not.toBeNull();
  });

  it("Esc leaves the phase choices and lets Tab reach every control until another sample opens", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(document.body, { key: "Tab" });
    const first = document.activeElement;
    expect(fireEvent.keyDown(first!, { key: "Escape" })).toBe(false);
    expect(document.activeElement).toBe(document.body);
    expect(fireEvent.keyDown(document.body, { key: "Tab" })).toBe(true);
    (first as HTMLElement).focus();
    expect(fireEvent.keyDown(first!, { key: "Tab" })).toBe(true);
    fireEvent.keyDown(window, { key: "ArrowDown" });
    expect(fireEvent.keyDown(document.body, { key: "Tab" })).toBe(false);
    expect(document.activeElement?.closest("[data-call-choices]")).not.toBeNull();
  });

  it("gives Tab to the filter's matches and to a custom index's symmetry, never to the phase choices", () => {
    state.trace = makeTrace(1);
    renderAt("/experiments/1/workspace/10");
    const inChoices = (): boolean => document.activeElement?.closest("[data-call-choices]") != null;
    const filter = document.getElementById("workspace-filter")!;
    act(() => filter.focus());
    expect(fireEvent.keyDown(filter, { key: "Tab" })).toBe(false);
    expect(document.activeElement).toBe(filter);
    expect(inChoices()).toBe(false);
    expect(document.querySelector("[data-testid=workspace-samples] [data-highlighted]")).not.toBeNull();
    act(() => filter.blur());
    fireEvent.keyDown(window, { key: "i" });
    const panel = screen.getByTestId("workspace-custom-index");
    const pressed = (): string | null => within(panel).getAllByRole("button", { pressed: true })[0]!.textContent;
    const before = pressed();
    expect(fireEvent.keyDown(document.body, { key: "Tab" })).toBe(false);
    expect(pressed()).not.toBe(before);
    expect(inChoices()).toBe(false);
  });

  it("leaves Tab in a text field to the field: a plain input keeps native Tab, and the label picker keeps it for its own list", () => {
    renderAt("/experiments/1/workspace/10");
    const inChoices = (): boolean => document.activeElement?.closest("[data-call-choices]") != null;
    const plain = document.createElement("input");
    plain.type = "text";
    document.body.appendChild(plain);
    try {
      act(() => plain.focus());
      expect(fireEvent.keyDown(plain, { key: "Tab" })).toBe(true);
      expect(inChoices()).toBe(false);
    } finally {
      plain.remove();
    }
    fireEvent.keyDown(window, { key: "l" });
    const input = within(screen.getByRole("dialog", { name: "Add a label" })).getByLabelText("Label name");
    expect(fireEvent.keyDown(input, { key: "Tab" })).toBe(false);
    expect(input).toHaveFocus();
    expect(inChoices()).toBe(false);
  });

  it("does not trap Tab when no phase choice can be picked", () => {
    state.username = undefined;
    renderAt("/experiments/1/workspace/10");
    expect(fireEvent.keyDown(document.body, { key: "Tab" })).toBe(true);
  });

  it("on a short window with Call folded, Tab opens Call and lands on its first choice", () => {
    vi.stubGlobal("matchMedia", (query: string) => ({ matches: query === "(max-height: 860px)", media: query, addEventListener: () => {}, removeEventListener: () => {} }));
    try {
      renderAt("/experiments/1/workspace/10");
      const tools = screen.getByRole("complementary", { name: "Tools" });
      fireEvent.click(within(tools).getByRole("button", { name: /^Image/ }));
      expect(screen.queryByTestId("workspace-candidates")).toBeNull();
      fireEvent.keyDown(document.body, { key: "Tab", shiftKey: true });
      expect(document.activeElement).toBe(screen.getByTestId("workspace-form-factor"));
    } finally {
      vi.unstubAllGlobals();
    }
  });

  it("0 marks the frame as no phase (form factor), and 0 again clears it", () => {
    state.indices = [];
    renderAt("/experiments/1/workspace/10");
    const row = screen.getByTestId("workspace-form-factor");
    expect(row).toHaveAttribute("aria-pressed", "false");
    fireEvent.keyDown(window, { key: "0" });
    expect(calls).toEqual([{ hook: "state", boundTo: 10, arg: "form_factor" }]);
    calls.length = 0;
    state.assignment = { exposure_id: 101, state: "form_factor", members: [] };
    fireEvent.click(screen.getByRole("button", { name: "Image 1" }));
    fireEvent.click(screen.getByRole("button", { name: "Image 2" }));
    expect(screen.getByTestId("workspace-form-factor")).toHaveAttribute("aria-pressed", "true");
    fireEvent.click(screen.getByTestId("workspace-form-factor"));
    expect(calls).toEqual([{ hook: "state", boundTo: 10, arg: "null" }]);
  });

  it("form factor is the sample's call, set the same from any frame", () => {
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    fireEvent.keyDown(window, { key: "0" });
    expect(calls).toEqual([{ hook: "state", boundTo: 10, arg: "form_factor" }]);
  });

  it("form factor is not offered on a dropped frame or before a user is set", () => {
    state.username = undefined;
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByTestId("workspace-form-factor")).toBeDisabled();
    fireEvent.keyDown(window, { key: "0" });
    expect(calls).toEqual([]);
  });
  it("shows a same sample's frames under their own heading after this sample's, and picking one shows it", () => {
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    state.setups = [{ id: 1, name: "Short pipe", folder: "s1", included: true }, { id: 2, name: "Long pipe", folder: "s2", included: true }];
    state.corpus = state.corpus.map((s) => (s.id === 10 ? { ...s, setup_id: 1 } : s.id === 12 ? { ...s, setup_id: 2, assignment_state: "indexed" } : s));
    state.mateExposures = new Map([[12, [exposure(200, { sample_id: 12 }), exposure(201, { sample_id: 12, image_path: "/m.tif" }), exposure(202, { sample_id: 12 })]]]);
    renderAt("/experiments/1/workspace/10");
    const frames = screen.getByTestId("workspace-exposures");
    expect(within(frames).getAllByTestId("frame-group-header").map((h) => h.textContent)).toEqual(["A1 in Short pipe", "A3 in Long pipe"]);
    expect(within(frames).getAllByRole("button", { name: /^Image \d/ }).map((b) => b.getAttribute("aria-label"))).toEqual([
      "Image 1", "Image 2", "Image 1 of A3 in Long pipe", "Image 2 of A3 in Long pipe", "Image 3 of A3 in Long pipe",
    ]);
    expect(screen.getByText("Image 2 of 2")).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: "Image" }).closest("section")!.textContent).toContain("SampleA1 in Short pipe");
    fireEvent.click(screen.getByRole("button", { name: "Image 2 of A3 in Long pipe" }));
    expect(screen.getByRole("button", { name: "Image 2 of A3 in Long pipe" })).toHaveAttribute("aria-current", "true");
    expect(screen.getByText("Image 2 of 3")).toBeInTheDocument();
    expect(screen.getByRole("heading", { name: "Image" }).closest("section")!.textContent).toContain("SampleA3 in Long pipe");
    expect(screen.getByText("Same as A3 in Long pipe")).toBeInTheDocument();
    fireEvent.keyDown(window, { key: " " });
    expect(screen.getByTestId("workspace-detector-caption")).toHaveTextContent("Image 2 of 3, A3 in Long pipe");
    fireEvent.keyDown(window, { key: " " });
  });

  it("draws no bridge to a sample whose setup is hidden, and draws it again once the setup is shown", () => {
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    state.setups = [{ id: 1, name: "Short pipe", folder: "s1", included: true }, { id: 2, name: "Long pipe", folder: "s2", included: false }];
    state.corpus = state.corpus.map((s) => (s.id === 12 ? { ...s, setup_id: 2 } : { ...s, setup_id: 1 }));
    state.mateExposures = new Map([[12, [exposure(200, { sample_id: 12 })]]]);
    const { rerender } = renderAt("/experiments/1/workspace/10");
    const frames = screen.getByTestId("workspace-exposures");
    expect(within(frames).queryByTestId("frame-group-header")).toBeNull();
    expect(within(frames).queryByRole("button", { name: /^Image 1 of A3/ })).toBeNull();
    expect(screen.queryByText(/^Same as A3/)).toBeNull();
    expect(screen.getByTestId("workspace-samples").querySelector("[data-bridged]")).toBeNull();

    state.setups = state.setups.map((st) => ({ ...st, included: true }));
    rerender(page("/experiments/1/workspace/10"));
    expect(within(screen.getByTestId("workspace-exposures")).getAllByTestId("frame-group-header")).toHaveLength(2);
    expect(screen.getByText("Same as A3 in Long pipe")).toBeInTheDocument();
    expect(screen.getByTestId("workspace-samples").querySelector("[data-sample-id='10'] [data-bridged]")).not.toBeNull();
  });

  it("the kept count under Frames counts only this sample's own frames", () => {
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    state.corpus = state.corpus.map((s) => (s.id === 12 ? { ...s, assignment_state: "indexed" } : s));
    state.exposures = [exposure(100), exposure(101, { status: "rejected" })];
    state.mateExposures = new Map([[12, [exposure(200, { sample_id: 12 }), exposure(201, { sample_id: 12, status: "rejected" })]]]);
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByText("1 of 2 kept")).toBeInTheDocument();
  });

  it("arrow keys step from this sample's last frame into the same sample's frames", () => {
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    state.corpus = state.corpus.map((s) => (s.id === 12 ? { ...s, assignment_state: "indexed" } : s));
    state.mateExposures = new Map([[12, [exposure(200, { sample_id: 12 })]]]);
    renderAt("/experiments/1/workspace/10?exposure=101");
    fireEvent.keyDown(window, { key: "ArrowRight" });
    expect(screen.getByRole("button", { name: "Image 1 of A3" })).toHaveAttribute("aria-current", "true");
    fireEvent.keyDown(window, { key: "ArrowLeft" });
    expect(screen.getByRole("button", { name: "Image 2" })).toHaveAttribute("aria-current", "true");
  });

  it("dropping a same sample's frame drops it on that sample, and undo restores it there", async () => {
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    state.corpus = state.corpus.map((s) => (s.id === 12 ? { ...s, assignment_state: "indexed" } : s));
    state.mateExposures = new Map([[12, [exposure(200, { sample_id: 12 }), exposure(201, { sample_id: 12 })]]]);
    state.exposures = [...state.exposures, exposure(200, { sample_id: 12 }), exposure(201, { sample_id: 12 })];
    renderAt("/experiments/1/workspace/10");
    fireEvent.click(screen.getByRole("button", { name: "Image 1 of A3" }));
    const tile = screen.getByRole("button", { name: "Image 1 of A3" }).closest("li")!;
    fireEvent.click(within(tile).getByRole("button", { name: "Drop" }));
    expect(calls.filter((c) => c.hook === "op:setStatus").map((c) => c.arg)).toEqual([[12, 200, "rejected"]]);
    expect(calls.filter((c) => c.hook === "status")).toHaveLength(0);
    await act(async () => { fireEvent.keyDown(window, { key: "z", ctrlKey: true }); for (let i = 0; i < 10; i++) await Promise.resolve(); });
    expect(calls.filter((c) => c.hook === "op:setStatus").map((c) => c.arg)).toEqual([[12, 200, "rejected"], [12, 200, null]]);
  });

  it("an unbridged sample's frames have no group heading", () => {
    renderAt("/experiments/1/workspace/10");
    expect(within(screen.getByTestId("workspace-exposures")).queryByTestId("frame-group-header")).toBeNull();
  });

  it("names the same sample in the heading and the list, and keeps it out of the comparisons", () => {
    state.trace = makeTrace(1);
    state.compareTraces = new Map([[12, makeTrace(3)]]);
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByText("Same as A3")).toBeInTheDocument();
    expect(screen.getByTestId("workspace-compares").querySelector("[data-compare='12']")).toBeNull();
    expect(screen.getByTestId("workspace-compare-bar").querySelector("[data-compare-chip]")).toBeNull();
    const list = screen.getByTestId("workspace-samples");
    expect(list.querySelector("[data-sample-id='10'] [data-bridged]")).toHaveTextContent("1 bridge");
    expect(list.querySelector("[data-sample-id='12'] [data-bridged]")).not.toBeNull();
    expect(list.querySelector("[data-sample-id='11'] [data-bridged]")).toBeNull();
    expect(screen.queryByRole("button", { name: /^Minus / })).toBeNull();
  });

  it("the palette's same-sample page leaves out this sample and its group, and Enter links the pick", () => {
    state.corpus = state.corpus.map((s) => (s.id === 12 ? { ...s, assignment_state: "indexed" } : s));
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    runFromPalette("same sample as…");
    const page = screen.getByRole("dialog", { name: "Same sample as" });
    expect(within(page).getAllByRole("option").map((o) => o.textContent)).toEqual(["A2", "A4"]);
    const input = within(page).getByRole("combobox", { name: "Sample" });
    fireEvent.change(input, { target: { value: "a4" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(calls).toContainEqual({ hook: "bridge", boundTo: 1, arg: { sampleId: 10, toSampleId: 13 } });
  });

  it("bridging to a sample with a different call asks which call to keep and sends that choice", async () => {
    state.corpus = state.corpus.map((s) => (s.id === 13 ? { ...s, assignment_state: "indexed", phase: "Im3m", call_phases: ["Im3m"] } : s));
    renderAt("/experiments/1/workspace/10");
    runFromPalette("same sample as…");
    const page = screen.getByRole("dialog", { name: "Same sample as" });
    const input = within(page).getByRole("combobox", { name: "Sample" });
    fireEvent.change(input, { target: { value: "a4" } });
    fireEvent.keyDown(input, { key: "Enter" });
    const form = screen.getByRole("dialog", { name: "Same sample as A4" });
    const keep = within(form).getByRole("combobox", { name: "Keep the call of" });
    expect(keep).toHaveTextContent("A1: Pn3m");
    expect(calls.filter((c) => c.hook === "bridge")).toHaveLength(0);
    fireEvent.keyDown(keep, { key: "ArrowDown" });
    fireEvent.keyDown(keep, { key: "Enter" });
    expect(screen.queryByRole("dialog")).toBeNull();
    await waitFor(() => expect(calls).toContainEqual({ hook: "bridge", boundTo: 1, arg: { sampleId: 10, toSampleId: 13, keepCall: 13 } }));
  });

  it("opening a bridged sample whose call state differs from its mate's asks once which call to keep", async () => {
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    state.setups = [{ id: 1, name: "Short pipe", folder: "s1", included: true }, { id: 2, name: "Long pipe", folder: "s2", included: true }];
    state.corpus = state.corpus.map((s) => (s.id === 10 ? { ...s, setup_id: 1 } : s.id === 12 ? { ...s, setup_id: 2 } : s));
    renderAt("/experiments/1/workspace/10");
    const form = screen.getByRole("dialog", { name: "Calls differ: same sample as A3 in Long pipe" });
    const keep = within(form).getByRole("combobox", { name: "Keep the call of" });
    expect(keep).toHaveTextContent("A1 in Short pipe: Pn3m");
    fireEvent.keyDown(keep, { key: "ArrowDown" });
    fireEvent.keyDown(keep, { key: "Enter" });
    expect(screen.queryByRole("dialog")).toBeNull();
    await waitFor(() => expect(calls.filter((c) => c.hook === "op:bridge").map((c) => c.arg)).toEqual([[10, 12, 12]]));
    fireEvent.click(screen.getByRole("button", { name: "Image 1" }));
    expect(screen.queryByRole("dialog")).toBeNull();
  });

  it("keeping the form-factor call when settling an old bridge can be undone, which re-ticks the phases it cleared", async () => {
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    const keep = within(screen.getByRole("dialog", { name: "Calls differ: same sample as A3" })).getByRole("combobox", { name: "Keep the call of" });
    fireEvent.keyDown(keep, { key: "ArrowDown" });
    fireEvent.keyDown(keep, { key: "Enter" });
    await waitFor(() => expect(calls.filter((c) => c.hook === "op:bridge")).toHaveLength(1));
    server.callMembers.set(10, []);
    await act(async () => { fireEvent.keyDown(window, { key: "z", ctrlKey: true }); for (let i = 0; i < 10; i++) await Promise.resolve(); });
    await waitFor(() => expect(calls.filter((c) => c.hook === "op:addPhase").map((c) => c.arg)).toEqual([[101, 2]]));
  });

  it("Esc on the settle prompt closes the palette instead of opening the command list", () => {
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    const keep = within(screen.getByRole("dialog", { name: "Calls differ: same sample as A3" })).getByRole("combobox", { name: "Keep the call of" });
    fireEvent.keyDown(keep, { key: "Escape" });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(calls.filter((c) => c.hook === "op:bridge")).toHaveLength(0);
  });

  it("the palette says the open sample is not the same sample as every other sample in its group", () => {
    state.corpus = state.corpus.map((s) => (s.id === 12 ? { ...s, assignment_state: "indexed" } : s));
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }, { sample_id: 13, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "not the same" } });
    expect(screen.getByRole("option", { name: "Not the same sample as A3, A4" })).toBeInTheDocument();
    fireEvent.keyDown(screen.getByLabelText("Command or sample"), { key: "Enter" });
    expect(calls).toContainEqual({ hook: "unbridge", boundTo: 1, arg: { sampleId: 10 } });
  });

  it("the same sample takes no comparison slot, so a fifth comparison goes in past a bridge", () => {
    state.corpus = [...state.corpus, sample(14, "A5"), sample(15, "A6"), sample(16, "A7")];
    act(() => setCompares({ compares: [11, 13, 14, 15].map((r) => ({ sample_id: 10, ref_sample_id: r, added_by: "jc" })), subtractions: [] }));
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    calls.length = 0;
    fireEvent.click(within(screen.getByTestId("workspace-samples")).getByRole("button", { name: /A7/ }), { ctrlKey: true });
    expect(calls.filter((c) => c.hook === "patchCompares")).toHaveLength(1);
    expect(toasts.at(-1)).not.toBe(COMPARES_FULL);
    expect(screen.getByTestId("workspace-compare-bar").querySelectorAll("[data-compare-chip]")).toHaveLength(5);
  });

  it("tells same-named samples in two shown setups apart in the picker, the chips and the same-sample note", () => {
    state.setups = [
      { id: 1, name: "setup1", folder: "setup1", included: true },
      { id: 2, name: "setup2", folder: "setup2", included: true },
    ];
    state.corpus = [{ ...sample(10, "MO_H2O"), setup_id: 1 }, { ...sample(11, "MO_H2O"), setup_id: 2 }, { ...sample(12, "MO_PEG5"), setup_id: 2 }];
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 11, group_id: 4 }];
    act(() => setCompares({ compares: [{ sample_id: 10, ref_sample_id: 12, added_by: "jc" }], subtractions: [] }));
    renderAt("/experiments/1/workspace/10");
    expect(screen.getByText("Same as MO_H2O in setup2")).toBeInTheDocument();
    expect(within(screen.getByTestId("workspace-compare-bar")).getByRole("button", { name: "MO_PEG5, setup2" })).toBeInTheDocument();
    fireEvent.keyDown(window, { key: "c" });
    const picker = screen.getByRole("dialog", { name: "Compare with a sample" });
    expect(within(picker).getByRole("button", { name: "MO_H2O, setup2" })).toBeInTheDocument();
    fireEvent.change(within(picker).getByLabelText("Find a sample to compare"), { target: { value: "setup2" } });
    expect(within(picker).getAllByRole("button").map((b) => b.getAttribute("data-sample-id"))).toEqual(["11"]);
  });

  it("leaves a same sample whose setup is hidden out of the note and the leave command, though the bridge stays stored", () => {
    state.setups = [
      { id: 2, name: "setup2", folder: "setup2", included: true },
      { id: 1, name: "setup1", folder: "setup1", included: false },
    ];
    state.corpus = [{ ...sample(10, "B1"), setup_id: 2 }, { ...sample(12, "A1"), setup_id: 1 }];
    state.bridges = [{ sample_id: 10, group_id: 4 }, { sample_id: 12, group_id: 4 }];
    renderAt("/experiments/1/workspace/10");
    expect(screen.queryByText(/^Same as A1/)).toBeNull();
    fireEvent.keyDown(window, { key: "k", ctrlKey: true });
    fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "not the same" } });
    expect(screen.queryByRole("option", { name: /^Not the same sample as/ })).toBeNull();
    expect(calls.filter((c) => c.hook.includes("bridge"))).toEqual([]);
  });

  it("B is not a shortcut, the palette offers no unbridging without bridges, and no bridging before a user is set", () => {
    const options = (): string[] => {
      fireEvent.keyDown(window, { key: "k", ctrlKey: true });
      fireEvent.change(screen.getByLabelText("Command or sample"), { target: { value: "same sample" } });
      const ids = screen.queryAllByRole("option").map((o) => o.getAttribute("data-palette-item") ?? "");
      fireEvent.keyDown(screen.getByRole("dialog"), { key: "Escape" });
      return ids.filter((id) => !id.startsWith("bridge-"));
    };
    const first = renderAt("/experiments/1/workspace/11");
    fireEvent.keyDown(window, { key: "b" });
    fireEvent.keyDown(window, { key: "B", shiftKey: true });
    expect(screen.queryByRole("dialog")).toBeNull();
    expect(calls.filter((c) => c.hook === "unbridge")).toHaveLength(0);
    expect(options()).toEqual(["bridge"]);
    first.unmount();
    state.username = undefined;
    renderAt("/experiments/1/workspace/11");
    expect(options()).toEqual([]);
  });
});
