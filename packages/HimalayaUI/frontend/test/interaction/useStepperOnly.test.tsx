import { describe, it, expect, vi } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { useStepperOnly } from "../../src/print/interaction/useStepperOnly";

const IDS = [10, 20, 30];

describe("useStepperOnly", () => {
  it("onPrev calls onGo with the previous id", () => {
    const onGo = vi.fn();
    const { result } = renderHook(() =>
      useStepperOnly({ ids: IDS, currentId: 20, onGo, label: "Sample", testIdBase: "sample" }),
    );
    act(() => { result.current.onPrev(); });
    expect(onGo).toHaveBeenCalledWith(10);
  });

  it("onNext calls onGo with the next id", () => {
    const onGo = vi.fn();
    const { result } = renderHook(() =>
      useStepperOnly({ ids: IDS, currentId: 20, onGo, label: "Sample", testIdBase: "sample" }),
    );
    act(() => { result.current.onNext(); });
    expect(onGo).toHaveBeenCalledWith(30);
  });

  it("prevDisabled is true at the first id", () => {
    const { result } = renderHook(() =>
      useStepperOnly({ ids: IDS, currentId: 10, onGo: vi.fn(), label: "Sample", testIdBase: "sample" }),
    );
    expect(result.current.prevDisabled).toBe(true);
    expect(result.current.nextDisabled).toBe(false);
  });

  it("nextDisabled is true at the last id", () => {
    const { result } = renderHook(() =>
      useStepperOnly({ ids: IDS, currentId: 30, onGo: vi.fn(), label: "Sample", testIdBase: "sample" }),
    );
    expect(result.current.prevDisabled).toBe(false);
    expect(result.current.nextDisabled).toBe(true);
  });

  it("count reflects position correctly for a middle id", () => {
    const { result } = renderHook(() =>
      useStepperOnly({ ids: IDS, currentId: 20, onGo: vi.fn(), label: "Sample", testIdBase: "sample" }),
    );
    expect(result.current.count).toBe("2 / 3");
  });

  it("count is '0 / 0' when ids is empty", () => {
    const { result } = renderHook(() =>
      useStepperOnly({ ids: [], currentId: undefined, onGo: vi.fn(), label: "Sample", testIdBase: "sample" }),
    );
    expect(result.current.count).toBe("0 / 0");
  });

  it("onPrev is a no-op at the first id", () => {
    const onGo = vi.fn();
    const { result } = renderHook(() =>
      useStepperOnly({ ids: IDS, currentId: 10, onGo, label: "Sample", testIdBase: "sample" }),
    );
    act(() => { result.current.onPrev(); });
    expect(onGo).not.toHaveBeenCalled();
  });

  it("onNext is a no-op at the last id", () => {
    const onGo = vi.fn();
    const { result } = renderHook(() =>
      useStepperOnly({ ids: IDS, currentId: 30, onGo, label: "Sample", testIdBase: "sample" }),
    );
    act(() => { result.current.onNext(); });
    expect(onGo).not.toHaveBeenCalled();
  });

  it("returns the props label, axis, and testIdBase verbatim", () => {
    const { result } = renderHook(() =>
      useStepperOnly({ ids: IDS, currentId: 20, onGo: vi.fn(), label: "Frame", testIdBase: "frame", axis: "horizontal" }),
    );
    expect(result.current.label).toBe("Frame");
    expect(result.current.axis).toBe("horizontal");
    expect(result.current.testIdBase).toBe("frame");
  });
});
