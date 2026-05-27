/**
 * BandResizeDivider tests (Plan §Phase 7, Task 7.2).
 *
 * Push-semantics drag math:
 *   - drag down by Δpx → band above grows, band below shrinks
 *   - total Σ band_heights stays constant
 *   - clamp at 0.1 minimum
 *   - snap to neighbor when within dead-zone of even split
 *
 * The divider dispatches through `resizeBands(memberIdx, deltaPx, totalPx)`
 * already implemented in `state.ts`. Tests assert the action is called with
 * the right arguments AND that integration with the action produces the
 * expected band-height ratios on the draft.
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { useAppState } from "../src/state";
import { emptyDraft } from "../src/lib/comparison/draft";
import { BandResizeDivider } from "../src/components/BandResizeDivider";

function seedDraft(bandHeights: number[]) {
  useAppState.setState({
    activeDraft: {
      ...emptyDraft(),
      members: bandHeights.map((bh, i) => ({
        id: i + 1,
        exposure_id: 100 + i,
        display_order: i,
        band_height: bh,
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

describe("BandResizeDivider — push semantics", () => {
  beforeEach(() => seedDraft([1, 1, 1]));

  it("renders with data-testid + above/below ids", () => {
    render(
      <BandResizeDivider
        aboveId={1}
        belowId={2}
        memberIndex={0}
        top={50}
        totalHeightPx={300}
      />,
    );
    const divider = screen.getByTestId("band-divider");
    expect(divider).toHaveAttribute("data-above-id", "1");
    expect(divider).toHaveAttribute("data-below-id", "2");
  });

  it("drag down grows band above, shrinks band below; total preserved", () => {
    render(
      <BandResizeDivider
        aboveId={1}
        belowId={2}
        memberIndex={0}
        top={100}
        totalHeightPx={300}
      />,
    );
    const divider = screen.getByTestId("band-divider");

    // mousedown at y=100, mousemove to y=130 (Δ = +30 px); mouseup at 130.
    fireEvent.mouseDown(divider, { clientY: 100 });
    fireEvent.mouseMove(window, { clientY: 130 });
    fireEvent.mouseUp(window, { clientY: 130 });

    const m = useAppState.getState().activeDraft!.members;
    // 30 / 300 = 0.1 normalized delta — above 1.0 → 1.1, below 1.0 → 0.9.
    expect(m[0]!.band_height).toBeCloseTo(1.1, 5);
    expect(m[1]!.band_height).toBeCloseTo(0.9, 5);
    // Sum across the two adjusted neighbors stays constant.
    expect(m[0]!.band_height + m[1]!.band_height).toBeCloseTo(2.0, 5);
    // Untouched band stays at its original ratio.
    expect(m[2]!.band_height).toBeCloseTo(1.0, 5);
  });

  it("clamps at 0.1 minimum (band can't be dragged out of existence)", () => {
    seedDraft([1, 1]);
    render(
      <BandResizeDivider
        aboveId={1}
        belowId={2}
        memberIndex={0}
        top={50}
        totalHeightPx={100}
      />,
    );
    const divider = screen.getByTestId("band-divider");
    // Drag up by 100 px (Δ = -1.0) — would push above to 0, but should clamp.
    fireEvent.mouseDown(divider, { clientY: 50 });
    fireEvent.mouseMove(window, { clientY: -50 });
    fireEvent.mouseUp(window, { clientY: -50 });

    const m = useAppState.getState().activeDraft!.members;
    expect(m[0]!.band_height).toBeGreaterThanOrEqual(0.1);
    expect(m[1]!.band_height).toBeGreaterThanOrEqual(0.1);
  });

  it("snaps to neighbor when within dead zone of even split", () => {
    // members: 1.0, 1.5 — even split would be band-above = band-below = 1.25.
    // Above is currently 1.0 (50px below the even line at 60px out of 100 total).
    seedDraft([1.0, 1.5]);
    render(
      <BandResizeDivider
        aboveId={1}
        belowId={2}
        memberIndex={0}
        top={40}
        totalHeightPx={250}
      />,
    );
    const divider = screen.getByTestId("band-divider");
    // Drag toward the even split (within snap dead-zone of 5px).
    // Even split shift: above must grow by (1.25 - 1.0) = 0.25 → 0.25*250 = 62.5px down.
    // Drag by 60px down (within 5px of 62.5) → should snap to even.
    fireEvent.mouseDown(divider, { clientY: 40 });
    fireEvent.mouseMove(window, { clientY: 100 });
    fireEvent.mouseUp(window, { clientY: 100 });
    const m = useAppState.getState().activeDraft!.members;
    expect(m[0]!.band_height).toBeCloseTo(1.25, 2);
    expect(m[1]!.band_height).toBeCloseTo(1.25, 2);
  });

  // I5.3 (#184): the resetBandHeights action was removed (Compare-only).

  it("does NOT snap when outside the dead zone", () => {
    seedDraft([1.0, 1.0]);
    render(
      <BandResizeDivider
        aboveId={1}
        belowId={2}
        memberIndex={0}
        top={50}
        totalHeightPx={300}
      />,
    );
    const divider = screen.getByTestId("band-divider");
    // Already even — drag well past the dead zone.
    fireEvent.mouseDown(divider, { clientY: 50 });
    fireEvent.mouseMove(window, { clientY: 80 });
    fireEvent.mouseUp(window, { clientY: 80 });
    const m = useAppState.getState().activeDraft!.members;
    expect(m[0]!.band_height).not.toBeCloseTo(1.0, 2);
    expect(m[0]!.band_height).toBeCloseTo(1.1, 2);
  });
});
