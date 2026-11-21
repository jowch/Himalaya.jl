import { useCallback, useEffect, useRef, useState } from "react";

export interface UndoEntry {
  label: string;
  sampleId: number;
  exposureId: number | undefined;
  undo: () => void;
  redo: () => void;
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

export interface UndoControls {
  record: (entry: UndoEntry) => void;
  undo: () => UndoEntry | undefined;
  redo: () => UndoEntry | undefined;
  nextUndo: UndoEntry | undefined;
  nextRedo: UndoEntry | undefined;
}

export function useUndoHistory(resetKey: unknown): UndoControls {
  const ref = useRef<UndoHistory>(emptyHistory);
  const [history, setHistory] = useState<UndoHistory>(emptyHistory);
  const commit = useCallback((next: UndoHistory) => {
    ref.current = next;
    setHistory(next);
  }, []);
  useEffect(() => commit(emptyHistory), [resetKey, commit]);
  const record = useCallback((entry: UndoEntry) => commit(recordEntry(ref.current, entry)), [commit]);
  const undo = useCallback(() => {
    const step = stepBack(ref.current);
    step.entry?.undo();
    commit(step.history);
    return step.entry;
  }, [commit]);
  const redo = useCallback(() => {
    const step = stepForward(ref.current);
    step.entry?.redo();
    commit(step.history);
    return step.entry;
  }, [commit]);
  return { record, undo, redo, nextUndo: history.past.at(-1), nextRedo: history.future.at(-1) };
}
