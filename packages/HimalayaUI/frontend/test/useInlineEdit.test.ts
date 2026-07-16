import { renderHook, act } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { useInlineEdit } from "../src/hooks/useInlineEdit";

describe("useInlineEdit", () => {
  it("begin seeds the draft and opens the key", () => {
    const { result } = renderHook(() => useInlineEdit<string>(() => {}));
    act(() => result.current.begin("row-1", "alpha"));
    expect(result.current.editingKey).toBe("row-1");
    expect(result.current.draft).toBe("alpha");
  });

  it("commit trims and fires onCommit only when the value changed", () => {
    const onCommit = vi.fn();
    const { result } = renderHook(() => useInlineEdit<string>(onCommit));
    act(() => result.current.begin("k", "alpha"));
    act(() => result.current.setDraft("  beta  "));
    act(() => result.current.commit());
    expect(onCommit).toHaveBeenCalledWith("k", "beta");
    expect(result.current.editingKey).toBeNull();
  });

  it("commit does NOT fire onCommit for an unchanged value (no no-op write)", () => {
    const onCommit = vi.fn();
    const { result } = renderHook(() => useInlineEdit<string>(onCommit));
    act(() => result.current.begin("k", "alpha"));
    act(() => result.current.setDraft("  alpha  ")); // same after trim
    act(() => result.current.commit());
    expect(onCommit).not.toHaveBeenCalled();
    expect(result.current.editingKey).toBeNull();
  });

  it("commit does NOT fire onCommit when a whitespace-padded seed is left unchanged", () => {
    const onCommit = vi.fn();
    const { result } = renderHook(() => useInlineEdit<string>(onCommit));
    act(() => result.current.begin("k", "  alpha  ")); // padded seed
    // draft left exactly as seeded (user opened the field and hit Enter)
    act(() => result.current.commit());
    expect(onCommit).not.toHaveBeenCalled();
  });

  it("double commit (Enter then blur) fires onCommit exactly once", () => {
    const onCommit = vi.fn();
    const { result } = renderHook(() => useInlineEdit<string>(onCommit));
    act(() => result.current.begin("k", "alpha"));
    act(() => result.current.setDraft("beta"));
    act(() => {
      result.current.commit();
      result.current.commit(); // the blur right after Enter
    });
    expect(onCommit).toHaveBeenCalledTimes(1);
  });

  it("cancel closes without committing", () => {
    const onCommit = vi.fn();
    const { result } = renderHook(() => useInlineEdit<string>(onCommit));
    act(() => result.current.begin("k", "alpha"));
    act(() => result.current.setDraft("beta"));
    act(() => result.current.cancel());
    expect(onCommit).not.toHaveBeenCalled();
    expect(result.current.editingKey).toBeNull();
  });

  it("boolean-style callers can use a sentinel key", () => {
    const onCommit = vi.fn();
    const { result } = renderHook(() => useInlineEdit<true>(onCommit));
    act(() => result.current.begin(true, "x"));
    expect(result.current.editingKey).toBe(true);
    act(() => result.current.setDraft("y"));
    act(() => result.current.commit());
    expect(onCommit).toHaveBeenCalledWith(true, "y");
  });
});
