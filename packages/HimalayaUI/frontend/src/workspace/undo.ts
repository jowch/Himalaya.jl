import { useCallback, useEffect, useRef, useState } from "react";

export interface UndoEntry {
  label: string;
  sampleId: number;
  exposureId: number | undefined;
  undo: () => Promise<number>;
  redo: () => Promise<number>;
}

export interface UndoHistory {
  past: readonly UndoEntry[];
  future: readonly UndoEntry[];
}

export const UNDO_LIMIT = 200;

export const emptyHistory: UndoHistory = { past: [], future: [] };

export function recordEntry(history: UndoHistory, entry: UndoEntry): UndoHistory {
  return { past: [...history.past, entry].slice(-UNDO_LIMIT), future: [] };
}

export function stepBack(history: UndoHistory): { entry: UndoEntry | undefined; history: UndoHistory } {
  const entry = history.past.at(-1);
  if (!entry) return { entry: undefined, history };
  return { entry, history: { past: history.past.slice(0, -1), future: [...history.future, entry] } };
}

export function stepForward(history: UndoHistory): { entry: UndoEntry | undefined; history: UndoHistory } {
  const entry = history.future.at(-1);
  if (!entry) return { entry: undefined, history };
  return { entry, history: { past: [...history.past, entry], future: history.future.slice(0, -1) } };
}

export interface UndoStep {
  entry: UndoEntry;
  done: Promise<number>;
}

export interface UndoControls {
  record: (entry: UndoEntry) => void;
  undo: () => UndoStep | undefined;
  redo: () => UndoStep | undefined;
  nextUndo: UndoEntry | undefined;
  nextRedo: UndoEntry | undefined;
}

export function useUndoHistory(resetKey: unknown, whenSaved: () => Promise<void>): UndoControls {
  const ref = useRef<UndoHistory>(emptyHistory);
  const chain = useRef<Promise<unknown>>(Promise.resolve());
  const saved = useRef(whenSaved);
  saved.current = whenSaved;
  const [history, setHistory] = useState<UndoHistory>(emptyHistory);
  const commit = useCallback((next: UndoHistory) => {
    ref.current = next;
    setHistory(next);
  }, []);
  useEffect(() => commit(emptyHistory), [resetKey, commit]);
  const record = useCallback((entry: UndoEntry) => commit(recordEntry(ref.current, entry)), [commit]);
  const step = useCallback(
    (dir: "undo" | "redo"): UndoStep | undefined => {
      const moved = dir === "undo" ? stepBack(ref.current) : stepForward(ref.current);
      const entry = moved.entry;
      if (!entry) return undefined;
      commit(moved.history);
      const done = chain.current
        .then(() => saved.current())
        .then(() => (dir === "undo" ? entry.undo() : entry.redo()));
      chain.current = done.catch(() => undefined);
      return { entry, done };
    },
    [commit],
  );
  const undo = useCallback(() => step("undo"), [step]);
  const redo = useCallback(() => step("redo"), [step]);
  return { record, undo, redo, nextUndo: history.past.at(-1), nextRedo: history.future.at(-1) };
}
