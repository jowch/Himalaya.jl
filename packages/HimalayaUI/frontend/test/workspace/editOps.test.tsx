import { describe, it, expect, vi, beforeEach } from "vitest";
import { act, renderHook, waitFor } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Assignment, CompareState, IndexEntry, Peak, SampleBridge, SampleCall, SampleLabel } from "../../src/api";
import { queryKeys, useAddPeak, useEditOps } from "../../src/queries";
import { useUndoHistory, type UndoEntry } from "../../src/workspace/undo";
import { bridgeEntry, compareEntry, combinedEntry, customIndexEntry, droppedByEdit, excludeEntry, labelEntry, peakEntry, tickEntry, unbridgeEntry } from "../../src/workspace/undoEntries";
import { pendingDeferreds } from "../../src/lib/queue/deferred";

const server = {
  peaks: [] as Peak[],
  indices: new Map<number, IndexEntry[]>(),
  assignments: new Map<number, Assignment>(),
  labels: [] as SampleLabel[],
  bridges: [] as SampleBridge[],
  calls: new Map<number, SampleCall>(),
  compares: { compares: [], subtractions: [] } as CompareState,
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
    restoreLabelSuggestion: vi.fn(async () => []),
    setLabelColor: vi.fn(async () => []),
    bridgeSamples: vi.fn(async () => []),
    unbridgeSample: vi.fn(async () => []),
    listPeaks: vi.fn(async () => server.peaks),
    listIndices: vi.fn(async (eid: number) => server.indices.get(eid) ?? []),
    getAssignment: vi.fn(async (eid: number) => server.assignments.get(eid) ?? { exposure_id: eid, state: "null", members: [] }),
    listLabels: vi.fn(async () => server.labels),
    listBridges: vi.fn(async () => server.bridges),
    getSampleCall: vi.fn(async (sid: number) => server.calls.get(sid) ?? { sample_id: sid, state: "indexed", phase: null, call_phases: [], members: [] }),
    patchCompares: vi.fn(async () => server.compares),
    listCompares: vi.fn(async () => server.compares),
  };
});

const api = await import("../../src/api");
const mocked = <T,>(fn: T) => fn as unknown as ReturnType<typeof vi.fn>;

function index(id: number, phase: string, basis: number, kind: IndexEntry["kind"] = "auto"): IndexEntry {
  return {
    id, exposure_id: 101, phase, basis, score: 0.9, r_squared: 0.99, lattice_d: 100, 
    status: "candidate", kind, inputs_hash: null, peaks: [], predicted_q: [],
  } as IndexEntry;
}

function auto(id: number, q: number, excluded: boolean): Peak {
  return { id, exposure_id: 101, q, intensity: null, prominence: null, sharpness: null, source: "auto", excluded };
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
  server.calls = new Map();
  server.compares = { compares: [], subtractions: [] };
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

  it("undoing a peak edit that dropped the call ticks the phase again once the re-fit has landed", async () => {
    let land: (p: unknown) => void = () => {};
    mocked(api.setAutoPeakExcluded).mockImplementationOnce(() => new Promise((resolve) => { land = resolve; }));
    const { result } = setup();
    server.peaks = [auto(3, 0.0644, true)];
    server.assignments.set(101, { exposure_id: 101, state: "null", members: [] });
    act(() => result.current.history.record(excludeEntry(result.current.ops, {
      label: "turned off a peak", sampleId: 10, exposureId: 101, q: 0.0644, excluded: true,
      dropped: Promise.resolve([{ exposureId: 101, phase: "Im3m", basis: 0.0455, kind: "auto" }]),
    })));
    let done: Promise<number> | undefined;
    act(() => { done = result.current.history.undo()!.done; });
    await waitFor(() => expect(api.setAutoPeakExcluded).toHaveBeenCalled());
    expect(firstArgs(api.setAutoPeakExcluded)).toEqual([[3, false]]);
    await act(async () => { await Promise.resolve(); });
    expect(api.addAssignmentPhase).not.toHaveBeenCalled();

    server.indices.set(101, [index(70, "Im3m", 0.0456), index(71, "Pn3m", 0.07)]);
    await act(async () => { land({ event_id: 1, view_row_id: null }); });
    let sent = 0;
    await act(async () => { sent = await done!; });
    expect(sent).toBe(2);
    expect(firstArgs(api.addAssignmentPhase)).toEqual([[101, 70]]);
  });

  it("undoing a peak edit leaves the call alone when the phase is still ticked or no longer fits", async () => {
    const { result } = setup();
    server.peaks = [auto(3, 0.0644, true), manual(6, 0.2)];
    server.indices.set(101, [index(70, "Im3m", 0.0456), index(71, "Pn3m", 0.09)]);
    server.assignments.set(101, { exposure_id: 101, state: "indexed", members: [70] });
    act(() => result.current.history.record(excludeEntry(result.current.ops, {
      label: "turned off a peak", sampleId: 10, exposureId: 101, q: 0.0644, excluded: true,
      dropped: Promise.resolve([{ exposureId: 101, phase: "Im3m", basis: 0.0455, kind: "auto" }, { exposureId: 101, phase: "Pn3m", basis: 0.07, kind: "auto" }]),
    })));
    let sent = 0;
    await act(async () => { sent = await result.current.history.undo()!.done; });
    expect(sent).toBe(1);
    expect(api.addAssignmentPhase).not.toHaveBeenCalled();
  });

  it("undoing a removed peak brings back the call it dropped, and an undo the server already made touches no call", async () => {
    const { result } = setup();
    mocked(api.addPeak).mockImplementationOnce(async () => ({ event_id: 1, view_row_id: null, ...manual(6, 0.2) }));
    server.indices.set(101, [index(80, "Ia3d", 0.05)]);
    act(() => result.current.history.record(peakEntry(result.current.ops, {
      label: "removed a peak", sampleId: 10, exposureId: 101, q: 0.2, added: false,
      dropped: Promise.resolve([{ exposureId: 101, phase: "Ia3d", basis: 0.05, kind: "auto" }]),
    })));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.addAssignmentPhase)).toEqual([[101, 80]]);

    mocked(api.addAssignmentPhase).mockClear();
    server.peaks = [manual(6, 0.2)];
    act(() => result.current.history.record(peakEntry(result.current.ops, {
      label: "removed a peak", sampleId: 10, exposureId: 101, q: 0.2, added: false,
      dropped: Promise.resolve([{ exposureId: 101, phase: "Ia3d", basis: 0.05, kind: "auto" }]),
    })));
    let sent = -1;
    await act(async () => { sent = await result.current.history.undo()!.done; });
    expect(sent).toBe(0);
    expect(api.addAssignmentPhase).not.toHaveBeenCalled();
  });

  it("after the edit saves, names only the phases its re-fit dropped", async () => {
    let land: (p: unknown) => void = () => {};
    mocked(api.setAutoPeakExcluded).mockImplementationOnce(() => new Promise((resolve) => { land = resolve; }));
    const { result } = setup();
    act(() => result.current.ops.setPeakExcluded(101, 3, true));
    let dropped: Promise<readonly unknown[]> | undefined;
    act(() => {
      dropped = droppedByEdit(result.current.ops, [
        { exposureId: 101, phase: "Im3m", basis: 0.0455, kind: "auto" },
        { exposureId: 101, phase: "Pn3m", basis: 0.07, kind: "auto" },
      ]);
    });
    await waitFor(() => expect(api.setAutoPeakExcluded).toHaveBeenCalled());
    server.indices.set(101, [index(70, "Im3m", 0.0456), index(71, "Pn3m", 0.09)]);
    server.assignments.set(101, { exposure_id: 101, state: "indexed", members: [70] });
    await act(async () => { land({ event_id: 1, view_row_id: null }); });
    expect(await dropped).toEqual([{ exposureId: 101, phase: "Pn3m", basis: 0.07, kind: "auto" }]);
  });

  it("names no dropped phase when the server cannot be read after the edit, so undo re-ticks nothing", async () => {
    mocked(api.listIndices).mockRejectedValueOnce(new Error("offline"));
    const { result } = setup();
    expect(await droppedByEdit(result.current.ops, [{ exposureId: 101, phase: "Im3m", basis: 0.0455, kind: "auto" }])).toEqual([]);
  });

  it("undoing a peak edit leaves alone a phase the edit kept, even when another tab unticked it since", async () => {
    const { result } = setup();
    server.peaks = [auto(3, 0.0644, false)];
    server.indices.set(101, [index(70, "Im3m", 0.0456)]);
    server.assignments.set(101, { exposure_id: 101, state: "indexed", members: [70] });
    const dropped = droppedByEdit(result.current.ops, [{ exposureId: 101, phase: "Im3m", basis: 0.0455, kind: "auto" }]);
    expect(await dropped).toEqual([]);
    act(() => result.current.history.record(excludeEntry(result.current.ops, {
      label: "turned off a peak", sampleId: 10, exposureId: 101, q: 0.0644, excluded: true, dropped,
    })));
    server.peaks = [auto(3, 0.0644, true)];
    server.assignments.set(101, { exposure_id: 101, state: "null", members: [] });
    let sent = 0;
    await act(async () => { sent = await result.current.history.undo()!.done; });
    expect(sent).toBe(1);
    expect(api.addAssignmentPhase).not.toHaveBeenCalled();
  });

  it("redoing a peak edit repeats only the peak edit and leaves the restored call to the server's re-fit", async () => {
    const { result } = setup();
    server.peaks = [auto(3, 0.0644, true)];
    server.indices.set(101, [index(70, "Im3m", 0.0456)]);
    act(() => result.current.history.record(excludeEntry(result.current.ops, {
      label: "turned off a peak", sampleId: 10, exposureId: 101, q: 0.0644, excluded: true,
      dropped: Promise.resolve([{ exposureId: 101, phase: "Im3m", basis: 0.0455, kind: "auto" as const }]),
    })));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.addAssignmentPhase)).toEqual([[101, 70]]);
    server.peaks = [auto(3, 0.0644, false)];
    server.assignments.set(101, { exposure_id: 101, state: "indexed", members: [70] });
    let sent = 0;
    await act(async () => { sent = await result.current.history.redo()!.done; });
    expect(sent).toBe(1);
    expect(firstArgs(api.setAutoPeakExcluded)).toEqual([[3, false], [3, true]]);
    expect(api.removeAssignmentPhase).not.toHaveBeenCalled();
    expect(firstArgs(api.addAssignmentPhase)).toEqual([[101, 70]]);
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

  it("puts a confirmed or dismissed suggestion back to suggested instead of deleting it, and redoes the call", async () => {
    const { result } = setup();
    server.labels = [{ id: 1, sample_id: 11, label: "Buffer", state: "confirmed", source: "user" }];
    act(() => result.current.history.record(labelEntry(result.current.ops, {
      label: "confirmed", sampleId: 11, exposureId: undefined, tag: "Buffer", added: true, prev: "suggested",
    })));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.restoreLabelSuggestion)).toEqual([[11, "Buffer"]]);
    server.labels = [{ id: 1, sample_id: 11, label: "Buffer", state: "suggested", source: "ingest" }];
    await act(async () => { await result.current.history.redo()!.done; });
    expect(firstArgs(api.addLabel)).toEqual([[11, "Buffer"]]);

    mocked(api.restoreLabelSuggestion).mockClear();
    server.labels = [];
    act(() => result.current.history.record(labelEntry(result.current.ops, {
      label: "dismissed", sampleId: 12, exposureId: undefined, tag: "Calibrant", added: false, prev: "suggested",
    })));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.restoreLabelSuggestion)).toEqual([[12, "Calibrant"]]);
    expect(api.removeLabel).not.toHaveBeenCalled();
    server.labels = [{ id: 2, sample_id: 12, label: "Calibrant", state: "suggested", source: "ingest" }];
    await act(async () => { await result.current.history.redo()!.done; });
    expect(firstArgs(api.removeLabel)).toEqual([[12, "Calibrant"]]);
  });

  it("says a suggestion is already back when the server shows it suggested", async () => {
    const { result } = setup();
    server.labels = [{ id: 1, sample_id: 11, label: "buffer", state: "suggested", source: "ingest" }];
    act(() => result.current.history.record(labelEntry(result.current.ops, {
      label: "confirmed", sampleId: 11, exposureId: undefined, tag: "Buffer", added: true, prev: "suggested",
    })));
    let sent = -1;
    await act(async () => { sent = await result.current.history.undo()!.done; });
    expect(sent).toBe(0);
    expect(api.restoreLabelSuggestion).not.toHaveBeenCalled();
  });

  it("undoes the parts of one gesture together, last first, and redoes them in order", async () => {
    const { result } = setup();
    const order: string[] = [];
    const part = (name: string, sent: number): UndoEntry => ({
      label: name, sampleId: 11, exposureId: undefined,
      undo: async () => { order.push(`undo ${name}`); return sent; },
      redo: async () => { order.push(`redo ${name}`); return sent; },
    });
    act(() => result.current.history.record(combinedEntry({ label: "created", sampleId: 11, exposureId: undefined }, [part("colour", 1), part("add", 0)])));
    let sent = -1;
    await act(async () => { sent = await result.current.history.undo()!.done; });
    expect(sent).toBe(1);
    await act(async () => { await result.current.history.redo()!.done; });
    expect(order).toEqual(["undo add", "undo colour", "redo colour", "redo add"]);
  });

  it("splits a bridge back into the two groups it joined", async () => {
    const { result } = setup();
    const before: SampleBridge[] = [
      { sample_id: 1, group_id: 1 }, { sample_id: 2, group_id: 1 }, { sample_id: 3, group_id: 1 },
      { sample_id: 4, group_id: 2 }, { sample_id: 5, group_id: 2 },
    ];
    act(() => result.current.history.record(bridgeEntry(result.current.ops, {
      label: "bridged", sampleId: 1, exposureId: undefined, from: 1, to: 4, bridges: before, labels: [],
      keep: undefined, states: { from: "indexed", to: "indexed" }, lost: [],
    })));
    server.bridges = before.map((b) => ({ ...b, group_id: 1 }));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.unbridgeSample).map((c) => c[0])).toEqual([4, 5]);
    expect(firstArgs(api.bridgeSamples)).toEqual([[4, 5]]);
    expect(mocked(api.bridgeSamples).mock.calls[0]![2]).toBe("both");
  });

  it("undoing a bridge that kept one call brings back the other side's call state and frames' phases", async () => {
    const { result } = setup();
    const lost = [{ index_id: 21, exposure_id: 201, phase: "Im3m" }];
    act(() => result.current.history.record(bridgeEntry(result.current.ops, {
      label: "bridged", sampleId: 1, exposureId: undefined, from: 1, to: 2, bridges: [], labels: [],
      keep: 1, states: { from: "form_factor", to: "indexed" }, lost,
    })));
    server.bridges = [{ sample_id: 1, group_id: 1 }, { sample_id: 2, group_id: 1 }];
    server.calls.set(1, { sample_id: 1, state: "form_factor", phase: null, call_phases: [], members: [] });
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.unbridgeSample).map((c) => c[0])).toEqual([1]);
    expect(firstArgs(api.setSampleCallState)).toEqual([[2, "indexed"]]);
    expect(firstArgs(api.addAssignmentPhase)).toEqual([[201, 21]]);
    server.bridges = [];
    await act(async () => { await result.current.history.redo()!.done; });
    expect(mocked(api.bridgeSamples).mock.calls.at(-1)!.slice(0, 3)).toEqual([1, 2, 1]);
  });

  it("undoing a bridge takes back only the labels the bridge shared, and leaves flight paths and later labels", async () => {
    const { result } = setup();
    const l = (id: number, sample_id: number, label: string, state: SampleLabel["state"] = "confirmed"): SampleLabel =>
      ({ id, sample_id, label, state, source: state === "suggested" ? "ingest" : "user" });
    const before = [l(1, 1, "Buffer"), l(2, 2, "37C"), l(3, 2, "buffer", "suggested"), l(4, 2, "1.7 m")];
    act(() => result.current.history.record(bridgeEntry(result.current.ops, {
      label: "bridged", sampleId: 1, exposureId: undefined, from: 1, to: 2, bridges: [], labels: before,
      keep: undefined, states: { from: "indexed", to: "indexed" }, lost: [],
    })));
    server.bridges = [{ sample_id: 1, group_id: 1 }, { sample_id: 2, group_id: 1 }];
    server.labels = [l(1, 1, "Buffer"), l(6, 1, "37C"), l(7, 1, "later"), l(2, 2, "37C"), l(3, 2, "buffer"), l(4, 2, "1.7 m"), l(8, 2, "later")];
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.unbridgeSample).map((c) => c[0])).toEqual([1]);
    expect(firstArgs(api.removeLabel)).toEqual([[1, "37C"]]);
    expect(firstArgs(api.restoreLabelSuggestion)).toEqual([[2, "buffer"]]);
  });

  it("undoing a bridge brings back a cleared call that the other side's form factor replaced", async () => {
    const { result } = setup();
    act(() => result.current.history.record(bridgeEntry(result.current.ops, {
      label: "bridged", sampleId: 1, exposureId: undefined, from: 1, to: 2, bridges: [], labels: [],
      keep: undefined, states: { from: "null", to: "form_factor" }, lost: [],
    })));
    server.bridges = [{ sample_id: 1, group_id: 1 }, { sample_id: 2, group_id: 1 }];
    server.calls.set(1, { sample_id: 1, state: "form_factor", phase: null, call_phases: [], members: [] });
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.setSampleCallState)).toEqual([[1, "null"]]);
  });

  it("undoes and redoes a comparison only while the server still matches what the edit left", async () => {
    const { result } = setup();
    const at = { label: "compared", sampleId: 1, exposureId: undefined };
    act(() => result.current.history.record(compareEntry(result.current.ops, { ...at, edit: { kind: "add", ref: 2 } })));
    server.compares = { compares: [{ sample_id: 1, ref_sample_id: 2, added_by: "jc" }], subtractions: [] };
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.patchCompares)).toEqual([[1, { remove: [2] }]]);
    mocked(api.patchCompares).mockClear();
    server.compares = { compares: [{ sample_id: 1, ref_sample_id: 2, added_by: "ana" }], subtractions: [] };
    await act(async () => { await result.current.history.redo()!.done; });
    expect(api.patchCompares).not.toHaveBeenCalled();
  });

  it("undoing an added comparison sends nothing once someone else already removed it", async () => {
    const { result } = setup();
    const at = { label: "compared", sampleId: 1, exposureId: undefined };
    act(() => result.current.history.record(compareEntry(result.current.ops, { ...at, edit: { kind: "add", ref: 2 } })));
    server.compares = { compares: [{ sample_id: 1, ref_sample_id: 3, added_by: "ana" }], subtractions: [] };
    await act(async () => { await result.current.history.undo()!.done; });
    expect(api.patchCompares).not.toHaveBeenCalled();
  });

  it("undoing a removed comparison brings back its subtraction unless another one was set since", async () => {
    const { result } = setup();
    const at = { label: "stopped comparing", sampleId: 1, exposureId: undefined };
    act(() => result.current.history.record(compareEntry(result.current.ops, { ...at, edit: { kind: "remove", ref: 2, subtracted: true } })));
    server.compares = { compares: [], subtractions: [] };
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.patchCompares)).toEqual([[1, { add: [2], subtraction: 2 }]]);
    mocked(api.patchCompares).mockClear();
    act(() => result.current.history.record(compareEntry(result.current.ops, { ...at, edit: { kind: "remove", ref: 2, subtracted: true } })));
    server.compares = { compares: [{ sample_id: 1, ref_sample_id: 3, added_by: "ana" }], subtractions: [{ sample_id: 1, ref_sample_id: 3 }] };
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.patchCompares)).toEqual([[1, { add: [2] }]]);
  });

  it("undoes a subtraction change only while that subtraction still stands", async () => {
    const { result } = setup();
    const at = { label: "subtracted", sampleId: 1, exposureId: undefined };
    act(() => result.current.history.record(compareEntry(result.current.ops, { ...at, edit: { kind: "subtract", from: 3, to: 2 } })));
    server.compares = { compares: [], subtractions: [{ sample_id: 1, ref_sample_id: 2 }] };
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.patchCompares)).toEqual([[1, { subtraction: 3 }]]);
    mocked(api.patchCompares).mockClear();
    act(() => result.current.history.record(compareEntry(result.current.ops, { ...at, edit: { kind: "subtract", from: undefined, to: 2 } })));
    server.compares = { compares: [], subtractions: [{ sample_id: 1, ref_sample_id: 4 }] };
    await act(async () => { await result.current.history.undo()!.done; });
    expect(api.patchCompares).not.toHaveBeenCalled();
  });

  it("rejoins an unbridged sample to the group it left", async () => {
    const { result } = setup();
    const before: SampleBridge[] = [{ sample_id: 1, group_id: 1 }, { sample_id: 2, group_id: 1 }];
    act(() => result.current.history.record(unbridgeEntry(result.current.ops, {
      label: "unbridged", sampleId: 1, exposureId: undefined, target: 2, bridges: before,
    })));
    await act(async () => { await result.current.history.undo()!.done; });
    expect(firstArgs(api.bridgeSamples)).toEqual([[2, 1]]);
    expect(mocked(api.bridgeSamples).mock.calls[0]![2]).toBe("both");
  });
});
