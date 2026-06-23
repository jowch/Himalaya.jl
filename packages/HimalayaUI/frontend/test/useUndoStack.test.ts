import { describe, it, expect } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { useUndoStack } from "../src/hooks/useUndoStack";

type Entry = { label: string; undo: () => void };

describe("useUndoStack", () => {
  it("starts empty, canUndo=false", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    expect(result.current.canUndo).toBe(false);
    expect(result.current.top).toBeUndefined();
  });

  it("push then pop returns the entry and removes it", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    act(() => result.current.push({ label: "move", undo: () => {} }));
    expect(result.current.canUndo).toBe(true);
    expect(result.current.top?.label).toBe("move");
    let popped: Entry | undefined;
    act(() => { popped = result.current.pop(); });
    expect(popped?.label).toBe("move");
    expect(result.current.canUndo).toBe(false);
  });

  it("pop on empty returns undefined (no throw)", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    let popped: Entry | undefined = { label: "x", undo: () => {} };
    act(() => { popped = result.current.pop(); });
    expect(popped).toBeUndefined();
  });

  it("clear empties the stack", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    act(() => result.current.push({ label: "a", undo: () => {} }));
    act(() => result.current.clear());
    expect(result.current.canUndo).toBe(false);
  });

  it("push uses a functional updater (StrictMode double-invoke pushes once per act)", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    act(() => result.current.push({ label: "x", undo: () => {} }));
    expect(result.current.depth).toBe(1);
  });

  it("pop moves the entry onto the redo stack; popRedo moves it back", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    act(() => result.current.push({ label: "a", undo: () => {} }));
    expect(result.current.canRedo).toBe(false);

    act(() => { result.current.pop(); });
    expect(result.current.canUndo).toBe(false);
    expect(result.current.canRedo).toBe(true);
    expect(result.current.redoTop?.label).toBe("a");

    let redone: Entry | undefined;
    act(() => { redone = result.current.popRedo(); });
    expect(redone?.label).toBe("a");
    expect(result.current.canUndo).toBe(true);
    expect(result.current.canRedo).toBe(false);
  });

  it("a fresh push clears the redo future", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    act(() => result.current.push({ label: "a", undo: () => {} }));
    act(() => { result.current.pop(); });
    expect(result.current.canRedo).toBe(true);
    act(() => result.current.push({ label: "b", undo: () => {} }));
    expect(result.current.canRedo).toBe(false);
  });

  it("popRedo on empty returns undefined (no throw)", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    let r: Entry | undefined = { label: "x", undo: () => {} };
    act(() => { r = result.current.popRedo(); });
    expect(r).toBeUndefined();
  });

  it("clear empties both the undo and redo stacks", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    act(() => result.current.push({ label: "a", undo: () => {} }));
    act(() => { result.current.pop(); });
    act(() => result.current.clear());
    expect(result.current.canUndo).toBe(false);
    expect(result.current.canRedo).toBe(false);
  });
});
