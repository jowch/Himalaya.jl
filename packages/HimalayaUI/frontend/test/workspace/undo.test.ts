import { describe, it, expect, vi } from "vitest";
import { act, renderHook } from "@testing-library/react";
import { emptyHistory, recordEntry, stepBack, stepForward, UNDO_LIMIT, useUndoHistory, type UndoEntry, type UndoStep } from "../../src/workspace/undo";

function entry(label: string): UndoEntry {
  return { label, sampleId: 1, exposureId: 10, undo: vi.fn(async () => 1), redo: vi.fn(async () => 1) };
}

describe("undo history", () => {
  it("steps back through recorded edits newest first and forward again", () => {
    const a = entry("a");
    const b = entry("b");
    let h = recordEntry(recordEntry(emptyHistory, a), b);
    const back1 = stepBack(h);
    expect(back1.entry).toBe(b);
    const back2 = stepBack(back1.history);
    expect(back2.entry).toBe(a);
    expect(stepBack(back2.history).entry).toBeUndefined();
    h = back2.history;
    expect(stepForward(h).entry).toBe(a);
    expect(stepForward(stepForward(h).history).entry).toBe(b);
  });

  it("drops the redo branch once a new edit is recorded after undoing", () => {
    const back = stepBack(recordEntry(emptyHistory, entry("a")));
    const h = recordEntry(back.history, entry("c"));
    expect(stepForward(h).entry).toBeUndefined();
  });

  it("keeps only the most recent edits", () => {
    let h = emptyHistory;
    for (let i = 0; i < UNDO_LIMIT + 5; i++) h = recordEntry(h, entry(String(i)));
    expect(h.past).toHaveLength(UNDO_LIMIT);
    expect(h.past[0]!.label).toBe("5");
  });
});

const saved = () => Promise.resolve();

describe("useUndoHistory", () => {
  it("records an edit without running it, runs its inverse on undo and the edit again on redo", async () => {
    const { result } = renderHook(() => useUndoHistory(1, saved));
    const e = entry("drop");
    act(() => result.current.record(e));
    expect(e.redo).not.toHaveBeenCalled();
    expect(result.current.nextUndo).toBe(e);
    let step: UndoStep | undefined;
    act(() => { step = result.current.undo(); });
    expect(step?.entry).toBe(e);
    expect(result.current.nextUndo).toBeUndefined();
    expect(result.current.nextRedo).toBe(e);
    await act(async () => { await step!.done; });
    expect(e.undo).toHaveBeenCalledTimes(1);
    await act(async () => { await result.current.redo()!.done; });
    expect(e.redo).toHaveBeenCalledTimes(1);
    expect(result.current.nextUndo).toBe(e);
  });

  it("waits for saving before each step, and runs steps one after another", async () => {
    let release: () => void = () => {};
    const gate = vi.fn(() => new Promise<void>((r) => { release = r; }));
    const { result } = renderHook(() => useUndoHistory(1, gate));
    const a = entry("a");
    const b = entry("b");
    act(() => { result.current.record(a); result.current.record(b); });
    let first: UndoStep | undefined;
    let second: UndoStep | undefined;
    act(() => { first = result.current.undo(); second = result.current.undo(); });
    expect(first?.entry).toBe(b);
    expect(second?.entry).toBe(a);
    await act(async () => { await Promise.resolve(); });
    expect(b.undo).not.toHaveBeenCalled();
    await act(async () => { release(); await first!.done; });
    expect(b.undo).toHaveBeenCalledTimes(1);
    expect(a.undo).not.toHaveBeenCalled();
    await act(async () => { release(); await second!.done; });
    expect(a.undo).toHaveBeenCalledTimes(1);
  });

  it("does nothing when there is nothing to undo or redo", () => {
    const { result } = renderHook(() => useUndoHistory(1, saved));
    let r: UndoStep | undefined;
    act(() => { r = result.current.undo(); });
    expect(r).toBeUndefined();
    act(() => { r = result.current.redo(); });
    expect(r).toBeUndefined();
  });

  it("forgets the history when the experiment or user changes", () => {
    const { result, rerender } = renderHook(({ k }) => useUndoHistory(k, saved), { initialProps: { k: "1:jc" } });
    act(() => result.current.record(entry("a")));
    rerender({ k: "1:mk" });
    expect(result.current.nextUndo).toBeUndefined();
  });
});
