/**
 * Hover-driven phase coloring tests (Plan §Phase 9, Task 9.5).
 *
 * Coordinates between two components:
 *   - MemberMetaRow (hover source): mouseenter/mouseleave/click/focus/blur/key
 *     handlers dispatch through `setHighlightedCompareMemberId`.
 *   - MemberTraceLayer (render target): peak marks recolor when the
 *     member's id matches the highlight, but that's already covered by
 *     the Phase 6 `highlightedIndexId` test. Here we test the hover-source
 *     coordination + lifecycle (click-to-pin, Esc to clear, etc.).
 *
 * Members with `confirmed_index === null` have no hover affordance —
 * nothing to highlight; handlers are inert.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemberMetaRow } from "../src/components/MemberMetaRow";
import { useAppState } from "../src/state";
import type { SeriesMember } from "../src/api";

function makeMember(over: Partial<SeriesMember> = {}): SeriesMember {
  return {
    id: 1,
    series_id: 100,
    exposure_id: 42,
    display_order: 0,
    band_height: 1,
    y_offset: 0,
    normalization: "qwindow",
    color_override: null,
    label_override: null,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot: {
      effective_peaks: [
        { id: 11, q: 0.30, intensity: 50, sharpness: 1, source: "auto" },
      ],
      confirmed_index: {
        id: 7, phase: "Pn3m", lattice_d: 12, r_squared: 0.99, ngc: -1.5,
        peak_ids: [11],
      },
      analysis_inputs_hash: "abc",
    },
    is_stale: false,
    created_by: null,
    created_at: null,
    ...over,
  };
}

beforeEach(() => {
  useAppState.setState({ highlightedCompareMemberId: undefined });
});

describe("hover-driven phase coloring — review mode", () => {
  it("mouseenter sets the highlight to this member's id", () => {
    render(<MemberMetaRow member={makeMember({ id: 42 })} top={0} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />);
    const row = screen.getByTestId("member-meta-row");
    fireEvent.mouseEnter(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(42);
  });

  it("mouseleave clears the highlight when not pinned", () => {
    render(<MemberMetaRow member={makeMember({ id: 42 })} top={0} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />);
    const row = screen.getByTestId("member-meta-row");
    fireEvent.mouseEnter(row);
    fireEvent.mouseLeave(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBeUndefined();
  });

  it("clicking pins the highlight; mouseleave does NOT clear after pin", () => {
    render(<MemberMetaRow member={makeMember({ id: 42 })} top={0} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />);
    const row = screen.getByTestId("member-meta-row");
    fireEvent.click(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(42);
    // After pinning, leaving the row should NOT clear.
    fireEvent.mouseLeave(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(42);
  });

  it("clicking a pinned row unpins it (toggle)", () => {
    render(<MemberMetaRow member={makeMember({ id: 42 })} top={0} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />);
    const row = screen.getByTestId("member-meta-row");
    fireEvent.click(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(42);
    fireEvent.click(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBeUndefined();
  });

  it("focus + Enter pins; Esc clears", () => {
    render(<MemberMetaRow member={makeMember({ id: 42 })} top={0} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />);
    const row = screen.getByTestId("member-meta-row");
    expect(row).toHaveAttribute("tabindex", "0");
    fireEvent.focus(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(42);
    fireEvent.keyDown(row, { key: "Enter" });
    // Pin makes it sticky.
    fireEvent.blur(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(42);
    fireEvent.keyDown(row, { key: "Escape" });
    expect(useAppState.getState().highlightedCompareMemberId).toBeUndefined();
  });

  it("focus on a row sets the highlight (transient until pin)", () => {
    render(<MemberMetaRow member={makeMember({ id: 42 })} top={0} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />);
    const row = screen.getByTestId("member-meta-row");
    fireEvent.focus(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(42);
    // Without pin, blurring clears.
    fireEvent.blur(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBeUndefined();
  });
});

describe("hover-driven phase coloring — no confirmed_index → inert", () => {
  it("members with confirmed_index === null are NOT focusable and don't highlight on hover", () => {
    const m = makeMember({
      id: 99,
      snapshot: {
        effective_peaks: [],
        confirmed_index: null,
        analysis_inputs_hash: "abc",
      },
    });
    render(<MemberMetaRow member={m} top={0} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />);
    const row = screen.getByTestId("member-meta-row");
    expect(row).not.toHaveAttribute("tabindex", "0");
    fireEvent.mouseEnter(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBeUndefined();
  });

  it("members with snapshot === null also stay inert", () => {
    const m = makeMember({ id: 99, snapshot: null });
    render(<MemberMetaRow member={m} top={0} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />);
    const row = screen.getByTestId("member-meta-row");
    fireEvent.mouseEnter(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBeUndefined();
    // Click also doesn't pin.
    fireEvent.click(row);
    expect(useAppState.getState().highlightedCompareMemberId).toBeUndefined();
  });
});

describe("hover-driven phase coloring — switching members updates highlight", () => {
  it("hovering row B replaces highlight from row A", () => {
    const A = makeMember({ id: 1 });
    const B = makeMember({ id: 2 });
    render(
      <>
        <MemberMetaRow member={A} top={0} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />
        <MemberMetaRow member={B} top={50} height={50} mode="review" displayLabel="row-label" expanded={false} onToggleExpand={() => {}} />
      </>,
    );
    const rows = screen.getAllByTestId("member-meta-row");
    const rowA = rows[0]!;
    const rowB = rows[1]!;
    fireEvent.mouseEnter(rowA);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(1);
    fireEvent.mouseEnter(rowB);
    expect(useAppState.getState().highlightedCompareMemberId).toBe(2);
  });
});
