import { useCallback, useRef, useState } from "react";

/**
 * A generic LIFO undo stack. Entries are opaque to the hook — the caller
 * decides each entry's shape and provides an `apply` callback at undo time.
 *
 * Every state transition uses a FUNCTIONAL updater so React StrictMode's
 * double-invoke in dev does NOT double-push or double-pop (the SeriesScoping
 * lesson: an impure updater that reads a captured `stack` pushed twice).
 */
export interface UndoStack<T> {
  /** Push one entry. Functional-updater safe (StrictMode double-invoke pushes
   *  once per `act` because the updater is pure). */
  push: (entry: T) => void;
  /**
   * Remove and RETURN the most recent entry (or undefined when empty). The
   * caller runs the inverse effect OUTSIDE any setState — `pop` performs NO
   * side effect of its own, so a StrictMode double-invoke of the internal
   * updater cannot double-apply the caller's effect (the SeriesScoping lesson:
   * never run a side effect inside a setState updater).
   */
  pop: () => T | undefined;
  /** Empty the stack (e.g. on route leave or successful build). */
  clear: () => void;
  canUndo: boolean;
  /** The entry that `pop` would return next, or undefined when empty. */
  top: T | undefined;
  depth: number;
}

export function useUndoStack<T>(): UndoStack<T> {
  const [stack, setStack] = useState<T[]>([]);
  // A ref mirror lets `pop` return the popped entry synchronously without
  // running a side effect inside the setState updater.
  const ref = useRef<T[]>(stack);
  ref.current = stack;

  const push = useCallback((entry: T) => {
    setStack((prev) => [...prev, entry]);
  }, []);

  const pop = useCallback((): T | undefined => {
    const cur = ref.current;
    if (cur.length === 0) return undefined;
    const entry = cur[cur.length - 1] as T;
    setStack((prev) => prev.slice(0, -1)); // pure updater — no side effect
    return entry;
  }, []);

  const clear = useCallback(() => setStack([]), []);

  return {
    push,
    pop,
    clear,
    canUndo: stack.length > 0,
    top: stack.length > 0 ? (stack[stack.length - 1] as T) : undefined,
    depth: stack.length,
  };
}
