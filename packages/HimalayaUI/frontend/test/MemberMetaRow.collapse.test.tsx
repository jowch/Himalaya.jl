/**
 * MemberMetaRow collapse/expand tests — Compare UX Phase E, Task E-2.
 *
 * The row renders collapsed by default (label + meta + disclosure caret).
 * Clicking `member-meta-row-body` toggles expansion via the controlled
 * `onToggleExpand` prop — the row does NOT own the expand state.
 *
 * The single-expanded invariant lives on `MemberMetaGutter`, which holds
 * one `expandedMemberId` and threads `expanded` + `onToggleExpand` into
 * each row. Expanding row B collapses row A.
 */
import { useState } from "react";
import { describe, it, expect } from "vitest";
import { render, screen, fireEvent, cleanup } from "@testing-library/react";
import { MemberMetaRow } from "../src/components/MemberMetaRow";
import { MemberMetaGutter } from "../src/components/MemberMetaGutter";
import { useAppState } from "../src/state";
import { emptyDraft } from "../src/lib/comparison/draft";
import type { ComparisonMember } from "../src/api";

function makeMember(over: Partial<ComparisonMember> = {}): ComparisonMember {
  return {
    id: 1,
    comparison_id: 100,
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
        { id: 12, q: 0.50, intensity: 80, sharpness: 1, source: "auto" },
      ],
      confirmed_index: {
        id: 7, phase: "Pn3m", lattice_d: 12.345, r_squared: 0.992, ngc: -1.51,
        peak_ids: [11, 12],
      },
      analysis_inputs_hash: "abc123",
    },
    is_stale: false,
    created_by: null,
    created_at: null,
    ...over,
  };
}

/**
 * Controlled-prop harness mirroring `MemberMetaGutter`'s single-expanded
 * invariant for direct `MemberMetaRow` mounts that don't need the full
 * gutter (drag wiring, y-bands, etc.).
 */
function ControlledRow(props: { member?: ComparisonMember; mode?: "review" | "edit" }) {
  const member = props.member ?? makeMember();
  const [expandedId, setExpandedId] = useState<number | null>(null);
  return (
    <MemberMetaRow
      member={member}
      top={0}
      height={50}
      mode={props.mode ?? "review"}
      memberIndex={0}
      displayLabel="row-label"
      expanded={expandedId === member.id}
      onToggleExpand={() =>
        setExpandedId((cur) => (cur === member.id ? null : member.id))
      }
    />
  );
}

describe("MemberMetaRow collapse/expand — Compare UX E-2", () => {
  it("renders collapsed view by default", () => {
    render(<ControlledRow />);
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute(
      "data-expanded",
      "false",
    );
  });

  it("expands when clicked on the row body", () => {
    render(<ControlledRow />);
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute(
      "data-expanded",
      "true",
    );
  });

  it("collapses when clicked again", () => {
    render(<ControlledRow />);
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute(
      "data-expanded",
      "false",
    );
  });

  it("collapsed row hides the per-member control widgets", () => {
    render(<ControlledRow mode="edit" member={makeMember()} />);
    // Collapsed: edit-mode controls are not rendered.
    expect(screen.queryByTestId("member-meta-label-input")).toBeNull();
    expect(screen.queryByTestId("member-meta-detail")).toBeNull();
    // Expanding reveals them.
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    expect(screen.getByTestId("member-meta-label-input")).toBeInTheDocument();
  });

  it("carries data-member-id and data-interactable on the root", () => {
    render(<ControlledRow />);
    const row = screen.getByTestId("member-meta-row");
    expect(row).toHaveAttribute("data-member-id", "1");
    expect(row).toHaveAttribute("data-interactable", "expand");
  });

  it("only one row expanded at a time across siblings", () => {
    cleanup();
    useAppState.setState({
      activeDraft: {
        ...emptyDraft(),
        members: [
          { id: 1, exposure_id: 100, display_order: 0, band_height: 1, y_offset: 0,
            normalization: "qwindow", color_override: undefined, label_override: undefined,
            q_window_min: undefined, q_window_max: undefined, peak_display: undefined,
            snapshot: undefined },
          { id: 2, exposure_id: 200, display_order: 1, band_height: 1, y_offset: 0,
            normalization: "qwindow", color_override: undefined, label_override: undefined,
            q_window_min: undefined, q_window_max: undefined, peak_display: undefined,
            snapshot: undefined },
        ],
      },
    });
    const members: ComparisonMember[] = [
      makeMember({ id: 1, exposure_id: 100 }),
      makeMember({ id: 2, exposure_id: 200 }),
    ];
    const labelMap = new Map<number, string>(
      members.map((m) => [m.id, `row-${m.id}`]),
    );
    render(
      <MemberMetaGutter
        members={members}
        panelHeight={300}
        mode="review"
        displayLabelByMemberId={labelMap}
      />,
    );
    const rows = screen.getAllByTestId("member-meta-row");
    const rowA = rows.find((r) => r.getAttribute("data-member-id") === "1")!;
    const rowB = rows.find((r) => r.getAttribute("data-member-id") === "2")!;
    const bodyA = rowA.querySelector("[data-testid='member-meta-row-body']")!;
    const bodyB = rowB.querySelector("[data-testid='member-meta-row-body']")!;

    // Expand A.
    fireEvent.click(bodyA);
    expect(rowA).toHaveAttribute("data-expanded", "true");
    expect(rowB).toHaveAttribute("data-expanded", "false");

    // Expand B — A must collapse (single-expanded invariant).
    fireEvent.click(bodyB);
    expect(rowA).toHaveAttribute("data-expanded", "false");
    expect(rowB).toHaveAttribute("data-expanded", "true");
  });
});
