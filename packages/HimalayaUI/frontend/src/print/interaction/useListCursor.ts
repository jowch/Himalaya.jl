import { useCallback, useEffect, useMemo, useRef, useState } from "react";
import type { ListCursor, RowProps, CursorStepperProps } from "./types";
import { safeScrollIntoView } from "../../lib/safeScrollIntoView";

interface Opts {
  ids: number[];
  onActivate?: (id: number) => void;
  stepperLabel?: string;
  stepperTestIdBase?: string;
  axis?: "vertical" | "horizontal";
}

export function useListCursor(opts: Opts): ListCursor {
  const { ids, onActivate, stepperLabel = "Item", stepperTestIdBase = "item", axis = "vertical" } = opts;

  const [cursorId, setCursorId] = useState<number | null>(ids[0] ?? null);
  const [selected, setSelected] = useState<Set<number>>(() => new Set());
  const elements = useRef<Map<number, HTMLElement>>(new Map());

  // Keep the cursor on a live id. If the cursored item vanished (SSE remove),
  // fall back to the nearest surviving index; never silently jump to index 0.
  const prevIds = useRef<number[]>(ids);
  useEffect(() => {
    if (cursorId !== null && ids.includes(cursorId)) {
      prevIds.current = ids;
      return;
    }
    if (cursorId === null) {
      if (ids.length > 0) setCursorId(ids[0]!);
      prevIds.current = ids;
      return;
    }
    const oldIndex = prevIds.current.indexOf(cursorId);
    const fallback = ids[Math.min(oldIndex, ids.length - 1)] ?? ids[ids.length - 1] ?? null;
    setCursorId(fallback);
    prevIds.current = ids;
  }, [ids, cursorId]);

  // Selection is ID-based too: prune any selected id that left the list (an SSE
  // remove, or a re-scope that swaps the whole `ids` set). Without this a sample
  // selected under one scope would linger after switching scope with no row left
  // to clear it (the SA-STALESELECT hazard). Orthogonality is preserved — this
  // fires only on `ids` membership change, never on cursor movement.
  useEffect(() => {
    setSelected((prev) => {
      let changed = false;
      const next = new Set<number>();
      for (const id of prev) {
        if (ids.includes(id)) next.add(id);
        else changed = true;
      }
      return changed ? next : prev;
    });
  }, [ids]);

  // Cursor === DOM focus: move focus to the cursored row.
  // Guard: if a focused interactive child (e.g. an edit input or a grip button)
  // already lives inside the row, defer to it — stealing focus would blur the
  // child and fire its onBlur (commit / cancel), which is wrong.
  useEffect(() => {
    if (cursorId === null) return;
    const el = elements.current.get(cursorId);
    if (!el) return;
    if (el !== document.activeElement && el.contains(document.activeElement)) {
      safeScrollIntoView(el, { block: "nearest" });
      return;
    }
    el.focus();
    safeScrollIntoView(el, { block: "nearest" });
  }, [cursorId]);

  const indexOf = useCallback((id: number | null) => (id === null ? -1 : ids.indexOf(id)), [ids]);

  const setCursor = useCallback((id: number) => setCursorId(id), []);

  const moveBy = useCallback(
    (delta: number) => {
      if (ids.length === 0) return;
      const from = Math.max(0, indexOf(cursorId));
      const next = Math.min(Math.max(from + delta, 0), ids.length - 1);
      setCursorId(ids[next]!);
    },
    [ids, cursorId, indexOf],
  );

  const activate = useCallback(() => {
    if (cursorId !== null) onActivate?.(cursorId);
  }, [cursorId, onActivate]);

  const toggleSelect = useCallback(
    (id?: number) => {
      const target = id ?? cursorId;
      if (target === null) return;
      setSelected((prev) => {
        const next = new Set(prev);
        if (next.has(target)) next.delete(target);
        else next.add(target);
        return next;
      });
    },
    [cursorId],
  );

  const rowProps = useCallback(
    (id: number): RowProps => ({
      ref: (el) => {
        if (el) elements.current.set(id, el);
        else elements.current.delete(id);
      },
      tabIndex: id === cursorId ? 0 : -1,
      onClick: () => setCursorId(id),
      onDoubleClick: () => {
        setCursorId(id);
        onActivate?.(id);
      },
      role: "row",
      "aria-current": id === cursorId ? "true" : undefined,
      "data-cursored": id === cursorId ? "true" : "false",
    }),
    [cursorId, onActivate],
  );

  const stepperProps = useCallback((): CursorStepperProps => {
    const i = indexOf(cursorId);
    return {
      label: stepperLabel,
      axis,
      testIdBase: stepperTestIdBase,
      count: ids.length > 0 ? `${Math.max(0, i) + 1} / ${ids.length}` : "0 / 0",
      onPrev: () => moveBy(-1),
      onNext: () => moveBy(1),
      prevDisabled: i <= 0,
      nextDisabled: i < 0 || i >= ids.length - 1,
    };
  }, [ids, cursorId, indexOf, moveBy, stepperLabel, stepperTestIdBase, axis]);

  return useMemo<ListCursor>(
    () => ({ cursorId, selected, setCursor, moveBy, activate, toggleSelect, rowProps, stepperProps }),
    [cursorId, selected, setCursor, moveBy, activate, toggleSelect, rowProps, stepperProps],
  );
}
