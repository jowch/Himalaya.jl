import { describe, it, expect, vi } from "vitest";
import { act, renderHook } from "@testing-library/react";
import { emptyHistory, recordEntry, stepBack, stepForward, UNDO_LIMIT, useUndoHistory, type UndoEntry } from "../../src/workspace/undo";

function entry(label: string): UndoEntry {
  return { label, sampleId: 1, exposureId: 10, undo: vi.fn(), redo: vi.fn() };
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

describe("useUndoHistory", () => {
  it("records an edit without running it, runs its inverse on undo and the edit again on redo", () => {
    const { result } = renderHook(() => useUndoHistory(1));
    const e = entry("drop");
    act(() => result.current.record(e));
    expect(e.redo).not.toHaveBeenCalled();
    expect(result.current.nextUndo).toBe(e);
    let undone: UndoEntry | undefined;
    act(() => { undone = result.current.undo(); });
    expect(undone).toBe(e);
    expect(e.undo).toHaveBeenCalledTimes(1);
    expect(result.current.nextUndo).toBeUndefined();
    expect(result.current.nextRedo).toBe(e);
    act(() => { result.current.redo(); });
    expect(e.redo).toHaveBeenCalledTimes(1);
    expect(result.current.nextUndo).toBe(e);
  });

  it("does nothing when there is nothing to undo or redo", () => {
    const { result } = renderHook(() => useUndoHistory(1));
    let r: UndoEntry | undefined = entry("x");
    act(() => { r = result.current.undo(); });
    expect(r).toBeUndefined();
    act(() => { r = result.current.redo(); });
    expect(r).toBeUndefined();
  });

  it("forgets the history when the experiment changes", () => {
    const { result, rerender } = renderHook(({ k }) => useUndoHistory(k), { initialProps: { k: 1 } });
    act(() => result.current.record(entry("a")));
    rerender({ k: 2 });
    expect(result.current.nextUndo).toBeUndefined();
  });
});
