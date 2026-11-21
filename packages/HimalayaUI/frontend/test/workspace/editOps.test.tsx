import { describe, it, expect, vi, beforeEach } from "vitest";
import { act, renderHook, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Assignment, IndexEntry, Peak, SampleBridge, SampleLabel } from "../../src/api";
import { queryKeys, useAddPeak, useEditOps } from "../../src/queries";
import { useUndoHistory } from "../../src/workspace/undo";
import { bridgeEntry, CANNOT_UNDO, customIndexEntry, labelEntry, peakEntry, tickEntry, unbridgeEntry } from "../../src/workspace/undoEntries";
import { pendingDeferreds } from "../../src/lib/queue/deferred";

const server = {
  peaks: [] as Peak[],
  indices: new Map<number, IndexEntry[]>(),
  assignments: new Map<number, Assignment>(),
  labels: [] as SampleLabel[],
  bridges: [] as SampleBridge[],
};

vi.mock("../../src/state", () => ({
  useAppState: (sel: (s: { username: string }) => unknown) => sel({ username: "jc" }),
}));

vi.mock("../../src/lib/toast", () => ({ showToast: () => {} }));

vi.mock("../../src/api", async (importOriginal) => {
  const actual = await importOriginal<typeof import("../../src/api")>();
  const ok = (extra: object = {}) => vi.fn(async () => ({ event_id: 1, view_row_id: null, ...extra }));
  return {
    ...actual,
    addAssignmentPhase: vi.fn(async (eid: number) => ({ event_id: 1, view_row_id: null, exposure_id: eid, state: "indexed", members: [] })),
    removeAssignmentPhase: vi.fn(async (eid: number) => ({ event_id: 1, view_row_id: null, exposure_id: eid, state: "indexed", members: [] })),
    setSampleCallState: vi.fn(async (sid: number, state: string) => ({
      event_id: 1, view_row_id: null, sample_id: sid, state, phase: null, call_phases: [], members: [],
    })),
    setExposureStatus: ok(),
    addPeak: vi.fn(),
    removeManualPeak: ok(),
    setAutoPeakExcluded: ok(),
    createCustomIndex: ok(),
    deleteIndex: vi.fn(async () => ({ deleted: 1 })),
    addLabel: vi.fn(async () => []),
    removeLabel: vi.fn(async () => []),
    setLabelColor: vi.fn(async () => []),
    bridgeSamples: vi.fn(async () => []),
    unbridgeSample: vi.fn(async () => []),
    listPeaks: vi.fn(async () => server.peaks),
    listIndices: vi.fn(async (eid: number) => server.indices.get(eid) ?? []),
    getAssignment: vi.fn(async (eid: number) => server.assignments.get(eid) ?? { exposure_id: eid, state: "null", members: [] }),
    listLabels: vi.fn(async () => server.labels),
    listBridges: vi.fn(async () => server.bridges),
  };
});

const api = await import("../../src/api");
const mocked = <T,>(fn: T) => fn as unknown as ReturnType<typeof vi.fn>;

function index(id: number, phase: string, basis: number, kind: IndexEntry["kind"] = "auto"): IndexEntry {
  return {
    id, exposure_id: 101, phase, basis, score: 0.9, r_squared: 0.99, lattice_d: 100, ngc: null,
    status: "candidate", kind, inputs_hash: null, peaks: [], predicted_q: [],
  } as IndexEntry;
}

function manual(id: number, q: number): Peak {
  return { id, exposure_id: 101, q, intensity: null, prominence: null, sharpness: null, source: "manual", excluded: false };
}

let qc: QueryClient;
function wrap({ children }: { children: ReactNode }) {
  return <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
}

beforeEach(() => {
  qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  pendingDeferreds.clear();
  vi.clearAllMocks();
  server.peaks = [];
  server.indices = new Map();
  server.assignments = new Map();
  server.labels = [];
  server.bridges = [];
});

const firstArgs = (fn: unknown) => mocked(fn).mock.calls.map((c) => c.slice(0, 2));

describe("useEditOps", () => {
  it("sends each edit to the exposure or sample it names, not the open one", async () => {
    const { result } = renderHook(() => useEditOps(1), { wrapper: wrap });
    const ops = result.current;
    act(() => {
      ops.addPhase(101, 5);
      ops.removePhase(102, 6);
      ops.setCallState(10, "form_factor");
      ops.setStatus(10, 103, "rejected");
      ops.removePeak(104, 7);
      ops.setPeakExcluded(104, 8, true);
      ops.commitCustom(105, "Ia3d", 0.05, [1, 2]);
      ops.deleteIndex(105, 9);
      ops.addLabel(11, "Buffer");
      ops.removeLabel(11, "Buffer");
      ops.bridge(11, 12);
      ops.unbridge(12);
    });
    await act(() => ops.whenSaved());
    expect(firstArgs(api.addAssignmentPhase)).toEqual([[101, 5]]);
    expect(firstArgs(api.removeAssignmentPhase)).toEqual([[102, 6]]);
    expect(firstArgs(api.setSampleCallState)).toEqual([[10, "form_factor"]]);
    expect(firstArgs(api.setExposureStatus)).toEqual([[103, "rejected"]]);
    expect(firstArgs(api.removeManualPeak)).toEqual([[7, expect.anything()]]);
    expect(firstArgs(api.setAutoPeakExcluded)).toEqual([[8, true]]);
    expect(firstArgs(api.createCustomIndex)).toEqual([[105, "Ia3d"]]);
    expect(firstArgs(api.deleteIndex)).toEqual([[9, expect.anything()]]);
    expect(firstArgs(api.addLabel)).toEqual([[11, "Buffer"]]);
    expect(firstArgs(api.removeLabel)).toEqual([[11, "Buffer"]]);
    expect(firstArgs(api.bridgeSamples)).toEqual([[11, 12]]);
    expect(firstArgs(api.unbridgeSample)).toEqual([[12, expect.anything()]]);
  });

  it("reads from the server, not a cache that a peak edit has left stale", async () => {
    qc.setQueryData(queryKeys.indices(101), [index(2, "Pn3m", 0.1)]);
    server.indices.set(101, [index(33, "Pn3m", 0.101)]);
    const { result } = renderHook(() => useEditOps(1), { wrapper: wrap });
    expect((await result.current.read.indices(101)).map((i) => i.id)).toEqual([33]);
  });
});

describe("undo against the server", () => {
  function setup() {
    return renderHook(() => {
      const ops = useEditOps(1);
      return { ops, add: useAddPeak(101), history: useUndoHistory(1, ops.whenSaved) };
    }, { wrapper: wrap });
  }

  it("unticks a phase whose index was renumbered by a peak edit, and ticks it again on redo", async () => {
    const { result } = setup();
    server.indices.set(101, [index(33, "Pn3m", 0.101), index(34, "Im3m", 0.2)]);
    server.assignments.set(101, { exposure_id: 101, state: "indexed", members: [33] });
    act(() => result.current.history.record(tickEntry(result.current.ops, {
      label: "ticked Pn3m", sampleId: 10, exposureId: 101, id: { exposureId: 101, phase: "Pn3m", basis: 0.1, kind: "auto" }, wasFormFactor: false,
    })));
    let sent = 0;
    await act(async () => { sent = await result.current.history.undo()!.done; });
    expect(sent).toBe(1);
    expect(firstArgs(api.removeAssignmentPhase)).toEqual([[101, 33]]);

    server.assignments.set(101, { exposure_id: 101, state: "null", members: [] });
    await act(async () => { await result.current.history.redo()!.done; });
    expect(firstArgs(api.addAssignmentPhase)).toEqual([[101, 33]]);
  });

  it("unticks the auto index, not a custom index of the same phase at a nearby lattice", async () => {
    const { result } = setup();
    server.indices.set(101, [index(50, "Pn3m", 0.1005, "speculative"), index(51, "Pn3m", 0.102)]);
    server.assignments.set(101, { exposure_id: 101, state: "indexed", members: [50, 51] });
    act(() => result.current.history.record(tickEntry(result.current.ops, {
      label: "ticked Pn3m", sampleId: 10, exposureId: 101, id: { exposureId: 101, phase: "Pn3m", basis: 0.1, kind: "auto" }, wasFormFactor: false,
    })));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.removeAssignmentPhase)).toEqual([[101, 51]]);
  });

  it("waits for a slow peak add to land, then removes the peak the server made", async () => {
    let land: (p: unknown) => void = () => {};
    mocked(api.addPeak).mockImplementation(() => new Promise((resolve) => { land = resolve; }));
    const { result } = setup();
    act(() => {
      result.current.add.mutate(0.2165);
      result.current.history.record(peakEntry(result.current.ops, {
        label: "added a peak", sampleId: 10, exposureId: 101, q: 0.2165, added: true,
      }));
    });
    let done: Promise<number> | undefined;
    act(() => { done = result.current.history.undo()!.done; });
    await act(async () => { await Promise.resolve(); });
    expect(api.listPeaks).not.toHaveBeenCalled();
    expect(api.removeManualPeak).not.toHaveBeenCalled();

    server.peaks = [manual(6, 0.2165)];
    await waitFor(() => expect(api.addPeak).toHaveBeenCalled());
    await act(async () => { land({ event_id: 1, view_row_id: null, ...manual(6, 0.2165) }); });
    let sent = 0;
    await act(async () => { sent = await done!; });
    expect(sent).toBe(1);
    expect(firstArgs(api.removeManualPeak)).toEqual([[6, expect.anything()]]);
  });

  it("reports nothing sent when the edit is already undone on the server", async () => {
    const { result } = setup();
    act(() => result.current.history.record(peakEntry(result.current.ops, {
      label: "added a peak", sampleId: 10, exposureId: 101, q: 0.2, added: true,
    })));
    let sent = -1;
    await act(async () => { sent = await result.current.history.undo()!.done; });
    expect(sent).toBe(0);
    expect(api.removeManualPeak).not.toHaveBeenCalled();
  });

  it("deletes a custom index on undo and adds it again on redo", async () => {
    const { result } = setup();
    server.indices.set(101, [index(40, "Ia3d", 0.05, "speculative")]);
    act(() => result.current.history.record(customIndexEntry(result.current.ops, {
      label: "added a custom Ia3d index", sampleId: 10, exposureId: 101, phase: "Ia3d", basis: 0.05, ratios: [1, 2], wasFormFactor: true,
    })));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.deleteIndex)).toEqual([[40, expect.anything()]]);
    expect(firstArgs(api.setSampleCallState)).toEqual([[10, "form_factor"]]);

    server.indices.set(101, []);
    await act(async () => { await result.current.history.redo()!.done; });
    expect(firstArgs(api.createCustomIndex)).toEqual([[101, "Ia3d"]]);
  });

  it("puts back a removed label and removes an added one", async () => {
    const { result } = setup();
    act(() => result.current.history.record(labelEntry(result.current.ops, {
      label: "removed Buffer", sampleId: 11, exposureId: undefined, tag: "Buffer", added: false, prev: "confirmed",
    })));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.addLabel)).toEqual([[11, "Buffer"]]);
    server.labels = [{ id: 1, sample_id: 11, label: "buffer", state: "confirmed", source: "user" }];
    await act(async () => { await result.current.history.redo()!.done; });
    expect(firstArgs(api.removeLabel)).toEqual([[11, "Buffer"]]);
  });

  it("never deletes a label to undo confirming or dismissing a suggestion", async () => {
    const { result } = setup();
    server.labels = [{ id: 1, sample_id: 11, label: "Buffer", state: "confirmed", source: "suggested" }];
    for (const added of [true, false]) {
      act(() => result.current.history.record(labelEntry(result.current.ops, {
        label: "suggestion", sampleId: 11, exposureId: undefined, tag: "Buffer", added, prev: "suggested",
      })));
      let sent = 0;
      await act(async () => { sent = await result.current.history.undo()!.done; });
      expect(sent).toBe(CANNOT_UNDO);
    }
    expect(api.removeLabel).not.toHaveBeenCalled();
    expect(api.addLabel).not.toHaveBeenCalled();
  });

  it("splits a bridge back into the two groups it joined", async () => {
    const { result } = setup();
    const before: SampleBridge[] = [
      { sample_id: 1, group_id: 1 }, { sample_id: 2, group_id: 1 }, { sample_id: 3, group_id: 1 },
      { sample_id: 4, group_id: 2 }, { sample_id: 5, group_id: 2 },
    ];
    act(() => result.current.history.record(bridgeEntry(result.current.ops, {
      label: "bridged", sampleId: 1, exposureId: undefined, from: 1, to: 4, bridges: before,
    })));
    server.bridges = before.map((b) => ({ ...b, group_id: 1 }));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.unbridgeSample).map((c) => c[0])).toEqual([4, 5]);
    expect(firstArgs(api.bridgeSamples)).toEqual([[4, 5]]);
  });

  it("rejoins an unbridged sample to the group it left", async () => {
    const { result } = setup();
    const before: SampleBridge[] = [{ sample_id: 1, group_id: 1 }, { sample_id: 2, group_id: 1 }];
    act(() => result.current.history.record(unbridgeEntry(result.current.ops, {
      label: "unbridged", sampleId: 1, exposureId: undefined, target: 2, bridges: before,
    })));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.bridgeSamples)).toEqual([[2, 1]]);
  });
});
