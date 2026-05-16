/**
 * Drag-vs-click disambiguation helper (spec §7.3 / Compare UX refinement).
 *
 * Usage:
 *   const state = makeDragThresholdState({ thresholdPx: 4 });
 *   onPointerDown:  state.onPointerDown(x, y)
 *   onPointerMove:  if (state.onPointerMove(x, y) === "drag-start") startDrag();
 *   onPointerUp:    const outcome = state.onPointerUp(x, y)
 *                    // outcome = "click" or "drag"
 *
 * Returns "click" when the pointer-up lands within `thresholdPx` Manhattan
 * distance of the down point AND no intervening move crossed the threshold.
 * Returns "drag" otherwise.
 */
export type DragOutcome = "click" | "drag";

export interface DragThresholdState {
  onPointerDown(x: number, y: number): void;
  onPointerMove(x: number, y: number): "drag-start" | "below-threshold";
  onPointerUp(x: number, y: number): DragOutcome;
  isDragging(): boolean;
  /**
   * Abandon the in-progress gesture without resolving it to a click or
   * drag. Called on `pointercancel` (native drag takeover, OS gesture,
   * touch interruption) — when no `pointerup` will follow — so the machine
   * is left in a neutral state ready for the next `onPointerDown`.
   */
  reset(): void;
}

export function makeDragThresholdState(
  opts: { thresholdPx: number },
): DragThresholdState {
  let downX: number | null = null;
  let downY: number | null = null;
  let dragging = false;

  return {
    onPointerDown(x, y) { downX = x; downY = y; dragging = false; },
    onPointerMove(x, y) {
      if (downX === null || downY === null) return "below-threshold";
      const dx = Math.abs(x - downX);
      const dy = Math.abs(y - downY);
      if (!dragging && (dx + dy) > opts.thresholdPx) {
        dragging = true;
        return "drag-start";
      }
      return "below-threshold";
    },
    onPointerUp(x, y) {
      const wasDragging = dragging;
      if (!wasDragging && downX !== null && downY !== null) {
        const dx = Math.abs(x - downX);
        const dy = Math.abs(y - downY);
        if ((dx + dy) > opts.thresholdPx) {
          dragging = true;
        }
      }
      const outcome: DragOutcome = dragging ? "drag" : "click";
      downX = null; downY = null; dragging = false;
      return outcome;
    },
    isDragging() { return dragging; },
    reset() { downX = null; downY = null; dragging = false; },
  };
}
