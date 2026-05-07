/**
 * Member reorder tests (Plan §Phase 7, Task 7.3).
 *
 * The drag-handle on the metadata row dispatches a reorder via
 * `reorderMembers(newOrder)` from Zustand. The action contiguously
 * re-numbers `display_order` (0, 1, 2, ...).
 *
 * The plot + gutter share `computeYBands` math, so testing reorder of the
 * draft's `members` array proves alignment for both halves: any consumer
 * iterating `members.sort(by display_order)` and applying `computeYBands`
 * will produce matching y-positions.
 *
 * The actual HTML5 drag-and-drop wiring lives in the gutter component
 * (Task 7.4 — `MemberMetaGutter`); this test exercises the
 * container-level drop handler via simulated drag events.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { useAppState } from "../src/state";
import { emptyDraft } from "../src/lib/comparison/draft";
import { MemberMetaGutter } from "../src/components/MemberMetaGutter";
import { computeYBands } from "../src/components/MultiTracePlot";
import type { ComparisonMember } from "../src/api";

function seedDraft(memberCount: number) {
  useAppState.setState({
    activeDraft: {
      ...emptyDraft(),
      members: Array.from({ length: memberCount }, (_, i) => ({
        id: i + 1,
        exposure_id: 100 + i,
        display_order: i,
        band_height: 1,
        y_offset: 0,
        normalization: "qwindow" as const,
        color_override: undefined,
        label_override: undefined,
        q_window_min: undefined,
        q_window_max: undefined,
        peak_display: undefined,
        snapshot: undefined,
      })),
    },
  });
}

function makeMembers(memberCount: number): ComparisonMember[] {
  return Array.from({ length: memberCount }, (_, i) => ({
    id: i + 1,
    comparison_id: 100,
    exposure_id: 100 + i,
    display_order: i,
    band_height: 1,
    y_offset: 0,
    normalization: "qwindow",
    color_override: null,
    label_override: null,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot: null,
    is_stale: false,
    created_by: null,
    created_at: null,
  }));
}

describe("Member reorder via metadata-row drag handle", () => {
  beforeEach(() => seedDraft(3));

  it("reorderMembers re-sequences display_order contiguously", () => {
    // Move row 2 → position 0. New order: [old idx 2, 0, 1].
    useAppState.getState().reorderMembers([2, 0, 1]);
    const m = useAppState.getState().activeDraft!.members;
    expect(m[0]!.id).toBe(3);  // member that was at idx 2
    expect(m[1]!.id).toBe(1);
    expect(m[2]!.id).toBe(2);
    expect(m.map((x) => x.display_order)).toEqual([0, 1, 2]);
  });

  it("MemberMetaGutter dispatches reorder when drag-and-drop completes", () => {
    render(
      <MemberMetaGutter
        members={makeMembers(3)}
        panelHeight={300}
        mode="edit"
      />,
    );
    // dragstart on grip-1 (member at idx 0); dragover/drop on grip-3 (idx 2).
    const grips = screen.getAllByTestId("member-reorder-grip");
    expect(grips).toHaveLength(3);

    const dt = {
      types: ["text/plain"],
      data: {} as Record<string, string>,
      setData(this: typeof dt, k: string, v: string) { this.data[k] = v; },
      getData(this: typeof dt, k: string) { return this.data[k] ?? ""; },
      effectAllowed: "move",
      dropEffect: "move",
    };
    fireEvent.dragStart(grips[0]!, { dataTransfer: dt });
    // Drop onto the row that hosts grip-3 (the third row).
    const targetRow = grips[2]!.closest("[data-testid='member-meta-row']")!;
    fireEvent.dragOver(targetRow, { dataTransfer: dt });
    fireEvent.drop(targetRow, { dataTransfer: dt });

    // Member 1 moved into position 2; members 2 and 3 shifted up.
    const m = useAppState.getState().activeDraft!.members;
    expect(m.map((x) => x.id)).toEqual([2, 3, 1]);
  });

  it("computeYBands produces aligned y-positions across plot + gutter after reorder", () => {
    // After reorder, both consumers (plot & gutter) call computeYBands with
    // the same `members` array (sorted by display_order). The bands match
    // by construction — this test pins it.
    const before = computeYBands([1, 1, 1], 300);
    expect(before.map(([t]) => t)).toEqual([0, 100, 200]);

    useAppState.getState().reorderMembers([2, 0, 1]);
    const m = useAppState.getState().activeDraft!.members;
    const after = computeYBands(m.map((mm) => mm.band_height), 300);
    // band_heights are all 1 → yBands shouldn't shift even after reorder.
    expect(after).toEqual(before);
  });

  it("does NOT reorder when dragging onto the same row", () => {
    render(
      <MemberMetaGutter
        members={makeMembers(3)}
        panelHeight={300}
        mode="edit"
      />,
    );
    const grips = screen.getAllByTestId("member-reorder-grip");
    const dt = {
      types: ["text/plain"],
      data: {} as Record<string, string>,
      setData(this: typeof dt, k: string, v: string) { this.data[k] = v; },
      getData(this: typeof dt, k: string) { return this.data[k] ?? ""; },
      effectAllowed: "move",
      dropEffect: "move",
    };
    fireEvent.dragStart(grips[1]!, { dataTransfer: dt });
    const targetRow = grips[1]!.closest("[data-testid='member-meta-row']")!;
    fireEvent.dragOver(targetRow, { dataTransfer: dt });
    fireEvent.drop(targetRow, { dataTransfer: dt });

    const m = useAppState.getState().activeDraft!.members;
    expect(m.map((x) => x.id)).toEqual([1, 2, 3]);  // unchanged
  });

  it("review mode hides drag handles", () => {
    render(
      <MemberMetaGutter
        members={makeMembers(3)}
        panelHeight={300}
        mode="review"
      />,
    );
    expect(screen.queryAllByTestId("member-reorder-grip")).toHaveLength(0);
  });
});
