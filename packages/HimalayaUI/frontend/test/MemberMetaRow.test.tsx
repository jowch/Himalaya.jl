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
import { useState } from "react";
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, fireEvent, act, cleanup } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import type { ComparisonMember } from "../src/api";
import { MemberMetaRow } from "../src/components/MemberMetaRow";
import type { MemberMetaRowProps } from "../src/components/MemberMetaRow";
import { useAppState } from "../src/state";
import { emptyDraft } from "../src/lib/comparison/draft";

/**
 * Compare UX E-2 — collapse/expand is a controlled prop on `MemberMetaRow`
 * (the gutter owns the single-expanded invariant). This harness threads a
 * local `expanded`/`onToggleExpand` pair so existing tests can mount the
 * row directly. `startExpanded` opens the row immediately for tests that
 * exercise the per-member control widgets (label / norm / q-window /
 * color picker), which now live behind the `{expanded && …}` block.
 */
type RowMountProps = Omit<MemberMetaRowProps, "expanded" | "onToggleExpand">;
function RowHarness({
  startExpanded = false,
  ...props
}: RowMountProps & { startExpanded?: boolean }): JSX.Element {
  const [expanded, setExpanded] = useState(startExpanded);
  return (
    <MemberMetaRow
      {...props}
      expanded={expanded}
      onToggleExpand={() => setExpanded((e) => !e)}
    />
  );
}

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
    render(<RowHarness member={makeMember()} top={0} height={50} mode="review" displayLabel="row-label" />);
    const row = screen.getByTestId("member-meta-row");
    expect(row).toHaveAttribute("data-member-id", "1");
  });

  it("shows phase chip / d / R² / NGC for cubic phase", () => {
    render(<RowHarness member={makeMember()} top={0} height={50} mode="review" displayLabel="row-label" />);
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
    render(<RowHarness member={m} top={0} height={50} mode="review" displayLabel="row-label" />);
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
    render(<RowHarness member={m} top={0} height={50} mode="review" displayLabel="row-label" />);
    expect(screen.getByText(/no index/i)).toBeInTheDocument();
  });

  it("data-stale present when member.is_stale", () => {
    const m = makeMember({ is_stale: true });
    render(<RowHarness member={m} top={0} height={50} mode="review" displayLabel="row-label" />);
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute("data-stale");
    // Inline stale indicator
    expect(screen.getByTestId("member-meta-stale-icon")).toBeInTheDocument();
  });

  it("data-stale absent when member.is_stale === false", () => {
    render(<RowHarness member={makeMember()} top={0} height={50} mode="review" displayLabel="row-label" />);
    expect(screen.getByTestId("member-meta-row")).not.toHaveAttribute("data-stale");
  });

  it("uses label_override when present, else falls back to a default label", () => {
    // The label_override fallback chain now lives in resolveDisplayLabels
    // (lib/comparison/labels.ts) — MemberMetaRow simply renders displayLabel.
    const m = makeMember({ label_override: "Sample A1 24h" });
    render(<RowHarness member={m} top={0} height={50} mode="review" displayLabel="Sample A1 24h" />);
    expect(screen.getByTestId("member-meta-label").textContent).toContain("Sample A1 24h");
  });

  it("renders the displayLabel prop verbatim (no internal fallback chain)", () => {
    render(
      <RowHarness
        member={makeMember()}
        top={0}
        height={50}
        mode="review"
        displayLabel="JC068P · run-007.dat"
      />,
    );
    expect(screen.getByTestId("member-meta-label").textContent).toBe(
      "JC068P · run-007.dat",
    );
  });

  it("clicking the row body expands a detail card overlay with peak count", () => {
    render(<RowHarness member={makeMember()} top={0} height={50} mode="review" displayLabel="row-label" />);
    // Compare UX E-2 — the disclosure affordance is the row body.
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    const detail = screen.getByTestId("member-meta-detail");
    expect(detail).toBeInTheDocument();
    expect(detail.textContent).toContain("2");  // 2 effective_peaks
    expect(detail.textContent?.toLowerCase()).toContain("peak");
  });

  it("positions the row at the supplied top + height", () => {
    render(<RowHarness member={makeMember()} top={120} height={45} mode="review" displayLabel="row-label" />);
    const row = screen.getByTestId("member-meta-row");
    expect(row.style.top).toBe("120px");
    expect(row.style.height).toBe("45px");
  });

  it("shows no edit-mode controls in review mode", () => {
    render(<RowHarness member={makeMember()} top={0} height={50} mode="review" displayLabel="row-label" />);
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
    render(<RowHarness member={makeMember()} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />);
    const grip = screen.getByTestId("member-reorder-grip");
    expect(grip).toBeInTheDocument();
    expect(grip).toHaveAttribute("data-member-id", "1");
  });

  it("label override input updates draft via updateMember", async () => {
    // Compare UX E-2 — label-override input lives in the expanded block.
    render(<RowHarness startExpanded member={makeMember()} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />);
    const input = screen.getByTestId("member-meta-label-input") as HTMLInputElement;
    await userEvent.type(input, "My label");
    fireEvent.blur(input);
    const m = useAppState.getState().activeDraft!.members[0]!;
    expect(m.label_override).toBe("My label");
  });

  it("normalization dropdown updates draft", () => {
    // Compare UX E-2 — normalization dropdown lives in the expanded block.
    render(<RowHarness startExpanded member={makeMember()} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />);
    const sel = screen.getByTestId("member-meta-normalization") as HTMLSelectElement;
    fireEvent.change(sel, { target: { value: "max" } });
    const m = useAppState.getState().activeDraft!.members[0]!;
    expect(m.normalization).toBe("max");
  });

  it("q-window inputs commit on blur (focus-gated, mid-edit external state ignored)", () => {
    // Compare UX E-2 — q-window inputs live behind the expanded block, so
    // the harness mounts the row pre-expanded. `RowHarness` keeps the
    // expand state across rerenders, so the input stays mounted.
    const m = makeMember({ q_window_min: 0.05, q_window_max: 0.5 });
    const { rerender } = render(
      <RowHarness startExpanded member={m} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />,
    );
    const minIn = screen.getByTestId("member-meta-qwindow-min") as HTMLInputElement;
    expect(minIn.value).toBe("0.050");

    // Focus the input and start editing — external value change should NOT
    // overwrite the draft while focused (QNumInput pattern).
    act(() => minIn.focus());
    fireEvent.change(minIn, { target: { value: "0.123" } });
    rerender(
      <RowHarness
        startExpanded
        member={makeMember({ q_window_min: 0.999, q_window_max: 0.5 })}
        top={0}
        height={50}
        mode="edit"
        memberIndex={0}
        displayLabel="row-label"
      />,
    );
    expect(minIn.value).toBe("0.123");

    // Blur commits to draft via updateMember
    fireEvent.blur(minIn);
    const draftMember = useAppState.getState().activeDraft!.members[0]!;
    expect(draftMember.q_window_min).toBeCloseTo(0.123);
  });

  it("'Reset color' button appears only when color_override is set", () => {
    // Compare UX E-2 — the "Reset color" button lives in the expanded
    // controls block, so the row is mounted pre-expanded.
    const m1 = makeMember({ color_override: null });
    const { rerender } = render(
      <RowHarness startExpanded member={m1} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />,
    );
    expect(screen.queryByTestId("member-meta-reset-color")).toBeNull();

    const m2 = makeMember({ color_override: "#ff00aa" });
    rerender(<RowHarness startExpanded member={m2} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />);
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
    // Compare UX E-2 — mount pre-expanded so the "Reset color" button shows.
    render(<RowHarness startExpanded member={m} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />);
    fireEvent.click(screen.getByTestId("member-meta-reset-color"));
    expect(useAppState.getState().activeDraft!.members[0]!.color_override).toBeUndefined();
  });

  // ─── Color-override picker (Phase 9, Task 9.4) ────────────────────────

  it("color picker swatch grid renders with palette colors", () => {
    render(<RowHarness member={makeMember()} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />);
    // Compare UX E-2 — expand via the row body (the disclosure affordance).
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    const grid = screen.getByTestId("member-color-picker-grid");
    expect(grid).toBeInTheDocument();
    const swatches = screen.getAllByTestId(/^member-color-picker-swatch-/);
    expect(swatches.length).toBeGreaterThanOrEqual(10);
    expect(swatches.length).toBeLessThanOrEqual(12);
  });

  it("clicking a swatch sets color_override to that hex/oklch color", () => {
    render(<RowHarness member={makeMember()} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />);
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    const swatches = screen.getAllByTestId(/^member-color-picker-swatch-/);
    const first = swatches[0]!;
    const colorValue = first.getAttribute("data-color")!;
    expect(colorValue).toBeTruthy();
    fireEvent.click(first);
    expect(useAppState.getState().activeDraft!.members[0]!.color_override).toBe(colorValue);
  });

  it("color picker swatch grid is hidden in review mode", () => {
    render(<RowHarness member={makeMember()} top={0} height={50} mode="review" displayLabel="row-label" />);
    // Even when expanded, review mode does not render the swatch grid.
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    expect(screen.queryByTestId("member-color-picker-grid")).toBeNull();
  });

  it("active swatch (matching color_override) is marked", () => {
    // First swatch's data-color is what we'll claim is the override.
    render(<RowHarness member={makeMember()} top={0} height={50} mode="edit" memberIndex={0} displayLabel="row-label" />);
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    const swatches = screen.getAllByTestId(/^member-color-picker-swatch-/);
    const targetColor = swatches[2]!.getAttribute("data-color")!;
    cleanup();
    // Re-render with that color set as override.
    render(
      <RowHarness
        member={makeMember({ color_override: targetColor })}
        top={0}
        height={50}
        mode="edit"
        memberIndex={0}
        displayLabel="row-label"
      />,
    );
    fireEvent.click(screen.getByTestId("member-meta-row-body"));
    const reSwatches = screen.getAllByTestId(/^member-color-picker-swatch-/);
    const target = reSwatches[2]!;
    expect(target).toHaveAttribute("data-active", "true");
    // Other swatches should be inactive.
    expect(reSwatches[0]!).toHaveAttribute("data-active", "false");
  });
});
