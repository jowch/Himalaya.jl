import { afterEach, describe, it, expect } from "vitest";
import { act, renderHook } from "@testing-library/react";
import { useSampleSort } from "../../src/workspace/sampleSortChoice";

afterEach(() => window.localStorage.clear());

describe("useSampleSort", () => {
  it("starts on measured order and remembers a choice for each experiment on its own", () => {
    const { result, rerender } = renderHook(({ id }) => useSampleSort(id), { initialProps: { id: 1 } });
    expect(result.current[0]).toBe("measured");
    act(() => result.current[1]("name"));
    expect(result.current[0]).toBe("name");
    rerender({ id: 2 });
    expect(result.current[0]).toBe("measured");
    rerender({ id: 1 });
    expect(result.current[0]).toBe("name");
  });

  it("forgets the stored choice when measured order is picked again", () => {
    const { result } = renderHook(() => useSampleSort(3));
    act(() => result.current[1]("name"));
    expect(window.localStorage.getItem("himalaya-ui:sample-sort:3")).toBe("name");
    act(() => result.current[1]("measured"));
    expect(window.localStorage.getItem("himalaya-ui:sample-sort:3")).toBeNull();
  });
});
