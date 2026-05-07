/**
 * peakCycle tests (Plan §Phase 8, Task 8.1).
 *
 * Cycle semantics for the edit-mode peak click:
 *
 *   shown (default) → labeled → hidden → shown
 *
 * `alt+click` jumps directly to `hidden` regardless of starting state.
 *
 * This file tests the pure transformation function (`cyclePeakDisplay`)
 * in isolation, plus the Zustand integration through `updateMember`. The
 * pure function is the load-bearing one — UI handlers are thin shims.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { cyclePeakDisplay } from "../src/lib/comparison/peakCycle";
import { useAppState } from "../src/state";
import { emptyDraft } from "../src/lib/comparison/draft";
import type { DraftMember } from "../src/lib/comparison/draft";

function makeMember(over: Partial<DraftMember> = {}): DraftMember {
  return {
    id: undefined,
    exposure_id: 42,
    display_order: 0,
    band_height: 1,
    y_offset: 0,
    normalization: "qwindow",
    color_override: undefined,
    label_override: undefined,
    q_window_min: undefined,
    q_window_max: undefined,
    peak_display: undefined,
    snapshot: undefined,
    ...over,
  };
}

// ── pure cyclePeakDisplay ────────────────────────────────────────────────────

describe("cyclePeakDisplay (pure function)", () => {
  it("shown → labeled (peak not in hidden, not in labeled) on regular click", () => {
    const out = cyclePeakDisplay(undefined, 11, false);
    expect(out.labeled).toEqual([11]);
    expect(out.hidden).toEqual([]);
  });

  it("treats null peak_display the same as undefined (shown default)", () => {
    const out = cyclePeakDisplay(undefined, 11, false);
    expect(out.labeled).toEqual([11]);
  });

  it("labeled → hidden (peak in labeled) on regular click", () => {
    const out = cyclePeakDisplay({ hidden: [], labeled: [11] }, 11, false);
    expect(out.labeled).toEqual([]);
    expect(out.hidden).toEqual([11]);
  });

  it("hidden → shown (peak in hidden) on regular click", () => {
    const out = cyclePeakDisplay({ hidden: [11], labeled: [] }, 11, false);
    expect(out.hidden).toEqual([]);
    expect(out.labeled).toEqual([]);
  });

  it("full cycle: shown → labeled → hidden → shown", () => {
    let s: { hidden: number[]; labeled: number[] } | undefined = undefined;
    s = cyclePeakDisplay(s, 11, false);
    expect(s).toEqual({ hidden: [], labeled: [11] });
    s = cyclePeakDisplay(s, 11, false);
    expect(s).toEqual({ hidden: [11], labeled: [] });
    s = cyclePeakDisplay(s, 11, false);
    expect(s).toEqual({ hidden: [], labeled: [] });
  });

  it("alt+click on shown peak → hidden directly", () => {
    const out = cyclePeakDisplay(undefined, 11, true);
    expect(out.hidden).toEqual([11]);
    expect(out.labeled).toEqual([]);
  });

  it("alt+click on labeled peak → hidden (removes from labeled)", () => {
    const out = cyclePeakDisplay({ hidden: [], labeled: [11] }, 11, true);
    expect(out.hidden).toEqual([11]);
    expect(out.labeled).toEqual([]);
  });

  it("alt+click on hidden peak → stays hidden (idempotent)", () => {
    const out = cyclePeakDisplay({ hidden: [11], labeled: [] }, 11, true);
    expect(out.hidden).toEqual([11]);
    expect(out.labeled).toEqual([]);
  });

  it("does not mutate the input — returns a new object", () => {
    const input = { hidden: [11], labeled: [12] };
    const out = cyclePeakDisplay(input, 13, false);
    expect(input).toEqual({ hidden: [11], labeled: [12] });
    expect(out).not.toBe(input);
    expect(out.hidden).not.toBe(input.hidden);
    expect(out.labeled).not.toBe(input.labeled);
  });

  it("preserves other peaks' display state when cycling one peak", () => {
    const out = cyclePeakDisplay(
      { hidden: [13], labeled: [14] },
      11,
      false,
    );
    expect(out.hidden).toContain(13);
    expect(out.labeled).toContain(14);
    expect(out.labeled).toContain(11);
  });
});

// ── Zustand integration ──────────────────────────────────────────────────────

describe("Zustand cyclePeakDisplayForMember action", () => {
  beforeEach(() => {
    // Reset Zustand draft state to a clean slate.
    const member = makeMember({ id: 1, exposure_id: 100 });
    const draft = { ...emptyDraft(), members: [member] };
    useAppState.setState({ activeDraft: draft });
  });

  it("persists cycle into draft.members[i].peak_display", () => {
    useAppState.getState().cyclePeakDisplayForMember(0, 11, false);
    const peakDisplay = useAppState.getState().activeDraft?.members[0]?.peak_display;
    expect(peakDisplay).toEqual({ hidden: [], labeled: [11] });
  });

  it("alt+click sets hidden directly", () => {
    useAppState.getState().cyclePeakDisplayForMember(0, 11, true);
    const peakDisplay = useAppState.getState().activeDraft?.members[0]?.peak_display;
    expect(peakDisplay).toEqual({ hidden: [11], labeled: [] });
  });

  it("multi-step cycle persists state correctly across calls", () => {
    const get = useAppState.getState;
    get().cyclePeakDisplayForMember(0, 11, false); // → labeled
    get().cyclePeakDisplayForMember(0, 11, false); // → hidden
    get().cyclePeakDisplayForMember(0, 11, false); // → shown
    const peakDisplay = get().activeDraft?.members[0]?.peak_display;
    expect(peakDisplay).toEqual({ hidden: [], labeled: [] });
  });

  it("no-op when memberIdx is out of range", () => {
    useAppState.getState().cyclePeakDisplayForMember(99, 11, false);
    const peakDisplay = useAppState.getState().activeDraft?.members[0]?.peak_display;
    expect(peakDisplay).toBeUndefined();
  });

  it("no-op when no active draft", () => {
    useAppState.setState({ activeDraft: null });
    expect(() =>
      useAppState.getState().cyclePeakDisplayForMember(0, 11, false),
    ).not.toThrow();
  });
});
