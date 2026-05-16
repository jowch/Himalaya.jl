import { describe, it, expect, afterEach } from "vitest";
import { renderHook } from "@testing-library/react";
import { useCompareMode } from "../src/hooks/useCompareMode";
import { useAppState } from "../src/state";
import type { Comparison } from "../src/api";
import type { ActiveDraft } from "../src/lib/comparison/draft";

const baseComparison: Comparison = {
  id: 1,
  title: "x",
  description: null,
  content_hash: "h",
  created_by: 5,
  created_at: null,
  updated_at: null,
  forked_from_id: null,
  forked_at_hash: null,
  forked_from_title: null,
  view_grouping_mode: null,
  view_show_peak_ticks: null,
  view_show_peak_labels: null,
  last_event_at: null,
  members: [],
};

const blankDraft: ActiveDraft = {
  id: undefined,
  baseHash: undefined,
  title: "n",
  description: "",
  members: [],
  forkedFromId: undefined,
  forkedAtHash: undefined,
  viewGroupingMode: undefined,
  viewShowPeakTicks: undefined,
  viewShowPeakLabels: undefined,
};

const editingDraft: ActiveDraft = {
  id: 1,
  baseHash: "h",
  title: "x",
  description: "",
  members: [],
  forkedFromId: undefined,
  forkedAtHash: undefined,
  viewGroupingMode: undefined,
  viewShowPeakTicks: undefined,
  viewShowPeakLabels: undefined,
};

describe("useCompareMode — Compare UX C-5", () => {
  afterEach(() => {
    useAppState.setState({ activeDraft: null });
  });

  it("returns 'viewing' when no draft + no comparison", () => {
    useAppState.setState({ activeDraft: null });
    const { result } = renderHook(() =>
      useCompareMode({ comparison: undefined, currentUserId: 5 }),
    );
    expect(result.current.kind).toBe("viewing");
  });

  it("returns 'viewing' when comparison loaded and no draft", () => {
    useAppState.setState({ activeDraft: null });
    const { result } = renderHook(() =>
      useCompareMode({ comparison: baseComparison, currentUserId: 5 }),
    );
    expect(result.current.kind).toBe("viewing");
  });

  it("returns 'editing-mine' when draft id matches and user is author", () => {
    useAppState.setState({ activeDraft: editingDraft });
    const { result } = renderHook(() =>
      useCompareMode({ comparison: baseComparison, currentUserId: 5 }),
    );
    expect(result.current.kind).toBe("editing-mine");
  });

  it("returns 'editing-as-fork-of' when draft id matches but user is NOT author", () => {
    useAppState.setState({ activeDraft: editingDraft });
    const { result } = renderHook(() =>
      useCompareMode({ comparison: baseComparison, currentUserId: 99 }),
    );
    expect(result.current.kind).toBe("editing-as-fork-of");
  });

  it("returns 'creating-blank' when draft id is undefined and no fork lineage", () => {
    useAppState.setState({ activeDraft: blankDraft });
    const { result } = renderHook(() =>
      useCompareMode({ comparison: undefined, currentUserId: 5 }),
    );
    expect(result.current.kind).toBe("creating-blank");
  });

  it("returns 'creating-from-fork' when draft id is undefined and forkedFromId is set", () => {
    useAppState.setState({
      activeDraft: {
        id: undefined,
        baseHash: undefined,
        title: "Copy of x",
        description: "",
        members: [],
        forkedFromId: 42,
        forkedAtHash: "h_parent",
        viewGroupingMode: undefined,
        viewShowPeakTicks: undefined,
        viewShowPeakLabels: undefined,
      },
    });
    const { result } = renderHook(() =>
      useCompareMode({ comparison: undefined, currentUserId: 5 }),
    );
    expect(result.current.kind).toBe("creating-from-fork");
    if (result.current.kind === "creating-from-fork") {
      expect(result.current.parentId).toBe(42);
    }
  });

  it("returns 'editing-as-fork-of' when draft has id but comparison is still loading (undefined)", () => {
    useAppState.setState({ activeDraft: editingDraft });
    const { result } = renderHook(() =>
      useCompareMode({ comparison: undefined, currentUserId: 5 }),
    );
    // comparison undefined → author is null → isAuthor is false → falls through to editing-as-fork-of
    expect(result.current.kind).toBe("editing-as-fork-of");
  });

  it("returns 'viewing-stale' when no draft and staleAgainstHash is set", () => {
    useAppState.setState({ activeDraft: null });
    const { result } = renderHook(() =>
      useCompareMode({
        comparison: baseComparison,
        currentUserId: 5,
        staleAgainstHash: "h_prev",
      }),
    );
    expect(result.current.kind).toBe("viewing-stale");
  });
});
