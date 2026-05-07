/**
 * MemberMetaRow tests (Plan §Phase 7, Task 7.1).
 *
 * Review-mode rendering: label / phase chip / d / R² / NGC for cubics,
 * stale flag, hover/click expansion.
 *
 * Edit-mode controls: drag handle, label override input, normalization
 * dropdown, q-window numeric inputs (focus-gated per `QNumInput` pattern),
 * "Reset color" button.
 *
 * Zustand wiring: edit-mode controls dispatch through `updateMember`.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import type { ComparisonMember } from "../src/api";
import { MemberMetaRow } from "../src/components/MemberMetaRow";
import { useAppState } from "../src/state";
import { emptyDraft } from "../src/lib/comparison/draft";

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

describe("MemberMetaRow — review mode", () => {
  it("renders the testid + member id", () => {
    render(<MemberMetaRow member={makeMember()} top={0} height={50} mode="review" />);
    const row = screen.getByTestId("member-meta-row");
    expect(row).toHaveAttribute("data-member-id", "1");
  });

  it("shows phase chip / d / R² / NGC for cubic phase", () => {
    render(<MemberMetaRow member={makeMember()} top={0} height={50} mode="review" />);
    expect(screen.getByText("Pn3m")).toBeInTheDocument();
    // lattice d = 12.345 → toFixed(2) → "12.35"
    expect(screen.getByText(/12\.35/)).toBeInTheDocument();
    expect(screen.getByText(/0\.992/)).toBeInTheDocument();   // R²
    // κ (NGC) appears for cubic phases only
    expect(screen.getByTestId("member-meta-ngc")).toBeInTheDocument();
  });

  it("hides NGC for non-cubic phases", () => {
    const m = makeMember({
      snapshot: {
        effective_peaks: [],
        confirmed_index: {
          id: 8, phase: "Hexagonal", lattice_d: 30, r_squared: 0.99, ngc: 0,
          peak_ids: [],
        },
        analysis_inputs_hash: "h",
      },
    });
    render(<MemberMetaRow member={m} top={0} height={50} mode="review" />);
    expect(screen.queryByTestId("member-meta-ngc")).toBeNull();
  });

  it("shows 'no index' when confirmed_index is null", () => {
    const m = makeMember({
      snapshot: {
        effective_peaks: [],
        confirmed_index: null,
        analysis_inputs_hash: "h",
      },
    });
    render(<MemberMetaRow member={m} top={0} height={50} mode="review" />);
    expect(screen.getByText(/no index/i)).toBeInTheDocument();
  });

  it("data-stale='true' when member.is_stale", () => {
    const m = makeMember({ is_stale: true });
    render(<MemberMetaRow member={m} top={0} height={50} mode="review" />);
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute("data-stale", "true");
    // Inline stale indicator
    expect(screen.getByTestId("member-meta-stale-icon")).toBeInTheDocument();
  });

  it("data-stale='false' when member.is_stale === false", () => {
    render(<MemberMetaRow member={makeMember()} top={0} height={50} mode="review" />);
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute("data-stale", "false");
  });

  it("uses label_override when present, else falls back to a default label", () => {
    const m = makeMember({ label_override: "Sample A1 24h" });
    render(<MemberMetaRow member={m} top={0} height={50} mode="review" />);
    expect(screen.getByTestId("member-meta-label").textContent).toContain("Sample A1 24h");
  });

  it("clicking the row expands a detail card overlay with peak count", () => {
    render(<MemberMetaRow member={makeMember()} top={0} height={50} mode="review" />);
    const row = screen.getByTestId("member-meta-row");
    fireEvent.click(row);
    const detail = screen.getByTestId("member-meta-detail");
    expect(detail).toBeInTheDocument();
    expect(detail.textContent).toContain("2");  // 2 effective_peaks
    expect(detail.textContent?.toLowerCase()).toContain("peak");
  });

  it("positions the row at the supplied top + height", () => {
    render(<MemberMetaRow member={makeMember()} top={120} height={45} mode="review" />);
    const row = screen.getByTestId("member-meta-row");
    expect(row.style.top).toBe("120px");
    expect(row.style.height).toBe("45px");
  });

  it("shows no edit-mode controls in review mode", () => {
    render(<MemberMetaRow member={makeMember()} top={0} height={50} mode="review" />);
    expect(screen.queryByTestId("member-reorder-grip")).toBeNull();
    expect(screen.queryByTestId("member-meta-label-input")).toBeNull();
    expect(screen.queryByTestId("member-meta-normalization")).toBeNull();
  });
});

describe("MemberMetaRow — edit mode", () => {
  beforeEach(() => {
    useAppState.setState({
      activeDraft: {
        ...emptyDraft(),
        members: [
          {
            id: 1,
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
          },
        ],
      },
    });
  });

  it("shows the drag handle", () => {
    render(<MemberMetaRow member={makeMember()} top={0} height={50} mode="edit" memberIndex={0} />);
    const grip = screen.getByTestId("member-reorder-grip");
    expect(grip).toBeInTheDocument();
    expect(grip).toHaveAttribute("data-member-id", "1");
  });

  it("label override input updates draft via updateMember", async () => {
    render(<MemberMetaRow member={makeMember()} top={0} height={50} mode="edit" memberIndex={0} />);
    const input = screen.getByTestId("member-meta-label-input") as HTMLInputElement;
    await userEvent.type(input, "My label");
    fireEvent.blur(input);
    const m = useAppState.getState().activeDraft!.members[0]!;
    expect(m.label_override).toBe("My label");
  });

  it("normalization dropdown updates draft", () => {
    render(<MemberMetaRow member={makeMember()} top={0} height={50} mode="edit" memberIndex={0} />);
    const sel = screen.getByTestId("member-meta-normalization") as HTMLSelectElement;
    fireEvent.change(sel, { target: { value: "max" } });
    const m = useAppState.getState().activeDraft!.members[0]!;
    expect(m.normalization).toBe("max");
  });

  it("q-window inputs commit on blur (focus-gated, mid-edit external state ignored)", () => {
    const m = makeMember({ q_window_min: 0.05, q_window_max: 0.5 });
    const { rerender } = render(
      <MemberMetaRow member={m} top={0} height={50} mode="edit" memberIndex={0} />,
    );
    const minIn = screen.getByTestId("member-meta-qwindow-min") as HTMLInputElement;
    expect(minIn.value).toBe("0.050");

    // Focus the input and start editing — external value change should NOT
    // overwrite the draft while focused (QNumInput pattern).
    act(() => minIn.focus());
    fireEvent.change(minIn, { target: { value: "0.123" } });
    rerender(
      <MemberMetaRow
        member={makeMember({ q_window_min: 0.999, q_window_max: 0.5 })}
        top={0}
        height={50}
        mode="edit"
        memberIndex={0}
      />,
    );
    expect(minIn.value).toBe("0.123");

    // Blur commits to draft via updateMember
    fireEvent.blur(minIn);
    const draftMember = useAppState.getState().activeDraft!.members[0]!;
    expect(draftMember.q_window_min).toBeCloseTo(0.123);
  });

  it("'Reset color' button appears only when color_override is set", () => {
    const m1 = makeMember({ color_override: null });
    const { rerender } = render(
      <MemberMetaRow member={m1} top={0} height={50} mode="edit" memberIndex={0} />,
    );
    expect(screen.queryByTestId("member-meta-reset-color")).toBeNull();

    const m2 = makeMember({ color_override: "#ff00aa" });
    rerender(<MemberMetaRow member={m2} top={0} height={50} mode="edit" memberIndex={0} />);
    expect(screen.getByTestId("member-meta-reset-color")).toBeInTheDocument();
  });

  it("clicking 'Reset color' clears the color_override in draft", () => {
    useAppState.setState({
      activeDraft: {
        ...useAppState.getState().activeDraft!,
        members: [
          {
            ...useAppState.getState().activeDraft!.members[0]!,
            color_override: "#ff00aa",
          },
        ],
      },
    });
    const m = makeMember({ color_override: "#ff00aa" });
    render(<MemberMetaRow member={m} top={0} height={50} mode="edit" memberIndex={0} />);
    fireEvent.click(screen.getByTestId("member-meta-reset-color"));
    expect(useAppState.getState().activeDraft!.members[0]!.color_override).toBeUndefined();
  });
});
