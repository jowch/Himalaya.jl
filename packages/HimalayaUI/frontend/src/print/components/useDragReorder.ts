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

/** Which edge of row `index` should show the drop indicator while dragging,
 *  or null. Matches reorder() placement: dragging from BELOW the target (a
 *  smaller drop index) inserts ABOVE it ("top"); from ABOVE inserts BELOW
 *  ("bottom"). No indicator on the dragged row itself or when not over it. */
export function dropEdgeFor(
  draggingIndex: number | null,
  overIndex: number | null,
  index: number,
): "top" | "bottom" | null {
  if (draggingIndex === null || overIndex === null) return null;
  if (overIndex !== index || draggingIndex === index) return null;
  return draggingIndex > index ? "top" : "bottom";
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
  overIndex: number | null;
  /** Spread onto each row wrapper: `<div {...dragItemProps(i)}>`. */
  dragItemProps: (index: number) => DragItemProps;
  /** Edge of row `index` to draw the insertion line on, or null. */
  dropEdge: (index: number) => "top" | "bottom" | null;
}

/** Native HTML5 drag-reorder wiring for a list whose order the caller owns.
 *  `onReorder(from, to)` is called on drop; the caller updates its order state
 *  (use the exported `reorder` helper). This hook is the POINTER path only:
 *  the keyboard path is the row's `onMoveBy` contract (see ScopeSampleRow),
 *  which the page must route into the SAME `onReorder` handler so both
 *  converge on one order mutation (done on SeriesScopingPage; the builder's
 *  adoption is tracked as BU-ORDERINERT). */
export function useDragReorder(onReorder: (from: number, to: number) => void): UseDragReorder {
  const [draggingIndex, setDraggingIndex] = useState<number | null>(null);
  const [overIndex, setOverIndex] = useState<number | null>(null);
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
      setOverIndex(index);
    },
    onDrop: (e) => {
      e.preventDefault();
      if (draggingIndex !== null && draggingIndex !== index) onReorder(draggingIndex, index);
      setDraggingIndex(null);
      setOverIndex(null);
    },
    onDragEnd: () => {
      setDraggingIndex(null);
      setOverIndex(null);
    },
    "data-dragging": draggingIndex === index,
  });
  return {
    draggingIndex,
    overIndex,
    dragItemProps,
    dropEdge: (index: number) => dropEdgeFor(draggingIndex, overIndex, index),
  };
}
