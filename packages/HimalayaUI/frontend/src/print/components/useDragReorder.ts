import { useState } from "react";
import type { DragEvent } from "react";

/** Pure: return a new array with the item at `from` moved to `to`. Out-of-range
 *  or no-op indices return the input unchanged. */
export function reorder<T>(list: readonly T[], from: number, to: number): T[] {
  if (from === to || from < 0 || to < 0 || from >= list.length || to >= list.length) {
    return list.slice();
  }
  const next = list.slice();
  const [moved] = next.splice(from, 1);
  next.splice(to, 0, moved as T);
  return next;
}

export interface DragItemProps {
  draggable: true;
  onDragStart: (e: DragEvent) => void;
  onDragOver: (e: DragEvent) => void;
  onDrop: (e: DragEvent) => void;
  onDragEnd: () => void;
  "data-dragging": boolean;
}

export interface UseDragReorder {
  draggingIndex: number | null;
  /** Spread onto each row wrapper: `<div {...dragItemProps(i)}>`. */
  dragItemProps: (index: number) => DragItemProps;
}

/** Native HTML5 drag-reorder wiring for a list whose order the caller owns.
 *  `onReorder(from, to)` is called on drop; the caller updates its order state
 *  (use the exported `reorder` helper). Pointer-only; keyboard reorder is a
 *  documented follow-up (the row grips are aria-hidden visual handles). */
export function useDragReorder(onReorder: (from: number, to: number) => void): UseDragReorder {
  const [draggingIndex, setDraggingIndex] = useState<number | null>(null);
  const dragItemProps = (index: number): DragItemProps => ({
    draggable: true,
    onDragStart: (e) => {
      setDraggingIndex(index);
      e.dataTransfer.effectAllowed = "move";
      // Some browsers require data to be set for a drag to start.
      try {
        e.dataTransfer.setData("text/plain", String(index));
      } catch {
        /* jsdom */
      }
    },
    onDragOver: (e) => {
      e.preventDefault();
      e.dataTransfer.dropEffect = "move";
    },
    onDrop: (e) => {
      e.preventDefault();
      if (draggingIndex !== null && draggingIndex !== index) onReorder(draggingIndex, index);
      setDraggingIndex(null);
    },
    onDragEnd: () => setDraggingIndex(null),
    "data-dragging": draggingIndex === index,
  });
  return { draggingIndex, dragItemProps };
}
