import { useCallback, useRef, useState } from "react";

/**
 * A generic LIFO undo/redo stack. Entries are opaque to the hook — the caller
 * decides each entry's shape and applies the effect itself: an undo entry must
 * carry BOTH the prior state (to undo) and the post-edit state (to redo), so
 * the caller can replay either direction from a returned entry.
 *
 * Every state transition uses a FUNCTIONAL updater so React StrictMode's
 * double-invoke in dev does NOT double-push or double-pop (the SeriesScoping
 * lesson: an impure updater that reads a captured `stack` pushed twice).
 */
export interface UndoStack<T> {
  /** Push one entry. Functional-updater safe (StrictMode double-invoke pushes
   *  once per `act` because the updater is pure). A new action invalidates the
   *  redo future, so this also clears the redo stack (standard semantics). */
  push: (entry: T) => void;
  /**
   * Remove and RETURN the most recent entry (or undefined when empty), moving
   * it onto the redo stack. The caller runs the inverse effect OUTSIDE any
   * setState — `pop` performs NO side effect of its own, so a StrictMode
   * double-invoke of the internal updater cannot double-apply the caller's
   * effect (the SeriesScoping lesson: never run a side effect inside a setState
   * updater).
   */
  pop: () => T | undefined;
  /**
   * Remove and RETURN the most recent REDONE entry (or undefined when the redo
   * stack is empty), moving it back onto the undo stack. Same no-side-effect
   * contract as `pop`: the caller replays the FORWARD effect outside any
   * setState.
   */
  popRedo: () => T | undefined;
  /** Empty BOTH stacks (e.g. on route leave or successful build). */
  clear: () => void;
  canUndo: boolean;
  canRedo: boolean;
  /** The entry that `pop` would return next, or undefined when empty. */
  top: T | undefined;
  /** The entry that `popRedo` would return next, or undefined when empty. */
  redoTop: T | undefined;
  depth: number;
}

export function useUndoStack<T>(): UndoStack<T> {
  const [stack, setStack] = useState<T[]>([]);
  const [redo, setRedo] = useState<T[]>([]);
  // Ref mirrors let `pop`/`popRedo` return the moved entry synchronously
  // without running a side effect inside a setState updater.
  const ref = useRef<T[]>(stack);
  ref.current = stack;
  const redoRef = useRef<T[]>(redo);
  redoRef.current = redo;

  const push = useCallback((entry: T) => {
    setStack((prev) => [...prev, entry]);
    setRedo([]); // a fresh action invalidates the redo future
  }, []);

  const pop = useCallback((): T | undefined => {
    const cur = ref.current;
    if (cur.length === 0) return undefined;
    const entry = cur[cur.length - 1] as T;
    setStack((prev) => prev.slice(0, -1)); // pure updaters — no side effect
    setRedo((prev) => [...prev, entry]);
    return entry;
  }, []);

  const popRedo = useCallback((): T | undefined => {
    const cur = redoRef.current;
    if (cur.length === 0) return undefined;
    const entry = cur[cur.length - 1] as T;
    setRedo((prev) => prev.slice(0, -1)); // pure updaters — no side effect
    setStack((prev) => [...prev, entry]);
    return entry;
  }, []);

  const clear = useCallback(() => {
    setStack([]);
    setRedo([]);
  }, []);

  return {
    push,
    pop,
    popRedo,
    clear,
    canUndo: stack.length > 0,
    canRedo: redo.length > 0,
    top: stack.length > 0 ? (stack[stack.length - 1] as T) : undefined,
    redoTop: redo.length > 0 ? (redo[redo.length - 1] as T) : undefined,
    depth: stack.length,
  };
}
