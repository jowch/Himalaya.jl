import { describe, it, expect } from "vitest";
import { renderHook } from "@testing-library/react";
import { useCompareDraftDirty } from "../src/hooks/useCompareDraftDirty";
import { useAppState } from "../src/state";

describe("useCompareDraftDirty — Compare UX C-6", () => {
  it("returns false when there is no draft", () => {
    useAppState.setState({ activeDraft: null });
    const { result } = renderHook(() => useCompareDraftDirty());
    expect(result.current).toBe(false);
  });

  it("returns true when there is any draft", () => {
    useAppState.setState({
      activeDraft: {
        id: undefined,
        baseHash: undefined,
        title: "x",
        description: "",
        members: [],
        forkedFromId: undefined,
        forkedAtHash: undefined,
        viewGroupingMode: undefined,
        viewShowPeakTicks: undefined,
        viewShowPeakLabels: undefined,
      },
    });
    const { result } = renderHook(() => useCompareDraftDirty());
    expect(result.current).toBe(true);
  });
});
