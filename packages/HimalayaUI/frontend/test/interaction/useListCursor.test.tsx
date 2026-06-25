import { describe, it, expect, vi } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { useListCursor } from "../../src/print/interaction/useListCursor";

describe("useListCursor", () => {
  it("defaults the cursor to the first id", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30] }));
    expect(result.current.cursorId).toBe(10);
  });

  it("moveBy steps by id and clamps at the ends", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30] }));
    act(() => result.current.moveBy(1));
    expect(result.current.cursorId).toBe(20);
    act(() => result.current.moveBy(-5));
    expect(result.current.cursorId).toBe(10);
    act(() => result.current.moveBy(99));
    expect(result.current.cursorId).toBe(30);
  });

  it("setCursor parks the cursor on a clicked id", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30] }));
    act(() => result.current.setCursor(30));
    expect(result.current.cursorId).toBe(30);
  });

  it("activate calls onActivate with the cursor id", () => {
    const onActivate = vi.fn();
    const { result } = renderHook(() => useListCursor({ ids: [10, 20], onActivate }));
    act(() => result.current.moveBy(1));
    act(() => result.current.activate());
    expect(onActivate).toHaveBeenCalledWith(20);
  });

  it("SSE-safety: removing an id before the cursor keeps the SAME item cursored", () => {
    const { result, rerender } = renderHook(({ ids }) => useListCursor({ ids }), {
      initialProps: { ids: [10, 20, 30] },
    });
    act(() => result.current.setCursor(30));
    rerender({ ids: [20, 30] }); // 10 removed (insert/remove before cursor)
    expect(result.current.cursorId).toBe(30); // still item 30, not index-shifted
  });

  it("SSE-safety: removing the cursored id falls back to the nearest surviving index", () => {
    const { result, rerender } = renderHook(({ ids }) => useListCursor({ ids }), {
      initialProps: { ids: [10, 20, 30] },
    });
    act(() => result.current.setCursor(20));
    rerender({ ids: [10, 30] }); // the cursored item itself removed
    expect(result.current.cursorId).toBe(30); // nearest surviving index: was index 1 → still index 1
  });

  it("toggleSelect maintains a Set orthogonal to the cursor", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30] }));
    act(() => result.current.toggleSelect(20));
    expect(result.current.selected.has(20)).toBe(true);
    expect(result.current.cursorId).toBe(10); // cursor did NOT move
    act(() => result.current.toggleSelect(20));
    expect(result.current.selected.has(20)).toBe(false);
  });

  it("stepperProps reflects position and end-disabled state", () => {
    const { result } = renderHook(() => useListCursor({ ids: [10, 20, 30], stepperLabel: "Sample", stepperTestIdBase: "sample" }));
    const s = result.current.stepperProps();
    expect(s.count).toBe("1 / 3");
    expect(s.prevDisabled).toBe(true);
    expect(s.nextDisabled).toBe(false);
  });
});
