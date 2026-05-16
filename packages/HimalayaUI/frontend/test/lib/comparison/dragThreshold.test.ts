import { describe, it, expect } from "vitest";
import { makeDragThresholdState, type DragOutcome } from "../../../src/lib/comparison/dragThreshold";

describe("dragThreshold — Compare UX C-1", () => {
  it("treats a move equal to the 4px threshold as click", () => {
    const s = makeDragThresholdState({ thresholdPx: 4 });
    s.onPointerDown(10, 10);
    s.onPointerMove(12, 12);
    const outcome: DragOutcome = s.onPointerUp(12, 12);
    expect(outcome).toBe("click");
  });

  it("treats > 4px move as drag", () => {
    const s = makeDragThresholdState({ thresholdPx: 4 });
    s.onPointerDown(10, 10);
    s.onPointerMove(15, 10);
    expect(s.onPointerUp(15, 10)).toBe("drag");
  });

  it("uses Manhattan distance for the threshold check (axis-agnostic)", () => {
    const s = makeDragThresholdState({ thresholdPx: 4 });
    s.onPointerDown(0, 0);
    // 3px in x AND 3px in y → euclidean ~4.24, manhattan 6 → drag either way
    s.onPointerMove(3, 3);
    expect(s.onPointerUp(3, 3)).toBe("drag");
  });

  it("reset() abandons the gesture — a later move without a fresh pointerDown is neutral", () => {
    const s = makeDragThresholdState({ thresholdPx: 4 });
    s.onPointerDown(10, 10);
    s.reset();
    expect(s.isDragging()).toBe(false);
    // No fresh onPointerDown — a move past the threshold must NOT register
    // as a drag because the machine is back to neutral (no down point).
    expect(s.onPointerMove(100, 100)).toBe("below-threshold");
    expect(s.isDragging()).toBe(false);
  });
});
