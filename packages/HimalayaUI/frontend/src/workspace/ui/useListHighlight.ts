import { useCallback, useEffect, useState, type KeyboardEvent } from "react";

interface Options {
  count: number;
  onCommit: (index: number) => void;
  onEscape: () => void;
  keys?: readonly (string | number)[];
  resetKey?: unknown;
  enabled?: boolean;
}

export interface ListHighlight {
  at: number;
  setAt: (index: number) => void;
  onKeyDown: (e: KeyboardEvent) => boolean;
}

export function stepHighlight(at: number, count: number, dir: 1 | -1): number {
  return count === 0 ? 0 : (at + dir + count) % count;
}

export function useListHighlight({ count, onCommit, onEscape, keys, resetKey, enabled = true }: Options): ListHighlight {
  const [cursor, setCursor] = useState<{ index: number; key: string | number | undefined }>({ index: 0, key: undefined });
  const kept = cursor.key === undefined || !keys ? -1 : keys.indexOf(cursor.key);
  const at = kept >= 0 ? kept : Math.min(cursor.index, Math.max(0, count - 1));
  const setAt = useCallback((index: number) => setCursor({ index, key: keys?.[index] }), [keys]);
  useEffect(() => setCursor({ index: 0, key: undefined }), [resetKey]);
  const onKeyDown = useCallback(
    (e: KeyboardEvent): boolean => {
      if (!enabled || e.metaKey || e.ctrlKey || e.altKey || e.nativeEvent.isComposing) return false;
      const move = e.key === "ArrowDown" || (e.key === "Tab" && !e.shiftKey) ? 1 : e.key === "ArrowUp" || (e.key === "Tab" && e.shiftKey) ? -1 : 0;
      if (move !== 0) {
        e.preventDefault();
        setAt(stepHighlight(at, count, move));
        return true;
      }
      if (e.key !== "Enter" && e.key !== "Escape") return false;
      e.preventDefault();
      if (e.repeat) return true;
      if (e.key === "Escape") onEscape();
      else if (count > 0) onCommit(at);
      return true;
    },
    [enabled, at, count, onCommit, onEscape, setAt],
  );
  return { at, setAt, onKeyDown };
}
