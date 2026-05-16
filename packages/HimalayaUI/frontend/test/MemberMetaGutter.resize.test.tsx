/**
 * Inter-row resize gap tests — Compare UX E-4.
 *
 * JSDOM has no layout engine: `getComputedStyle().height` returns "" for
 * elements without an inline style and `getBoundingClientRect` returns zero
 * rects. So the gap component drives its own height via INLINE `style.height`
 * ("4px" rest / "12px" hover) and exposes a component-controlled `data-state`
 * ("rest" | "hover"). Tests assert on those — never on Tailwind class strings
 * (CLAUDE.md / test/AGENTS.md).
 *
 * Band-tint coupling is verified through the real `ActiveBandContext`: the
 * gutter publishes the active member-id and a real `MemberBandOverlay`
 * subscriber (the same component `MultiTracePlot` mounts) sets
 * `data-active-band` when its `data-member-id` matches.
 *
 * JSDOM note: JSDOM's `PointerEvent` constructor does NOT carry
 * `clientX`/`clientY` from its init dict, so `fireEvent.pointerMove` delivers
 * `undefined` coordinates and the drag math sees `NaN`. `MouseEvent` *does*
 * carry them, and React's synthetic pointer handlers listen for the DOM
 * event *type*, so dispatching a `MouseEvent` with a pointer-event type
 * triggers the handler with real coordinates. Same `pointer()` helper as
 * `MemberMetaRow.drag.test.tsx`.
 */
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, fireEvent } from "@testing-library/react";
import { useAppState } from "../src/state";
import { emptyDraft } from "../src/lib/comparison/draft";
import type { ComparisonMember } from "../src/api";
import { MemberMetaGutter } from "../src/components/MemberMetaGutter";
import {
  ActiveBandProvider,
  useActiveBand,
} from "../src/components/ActiveBandContext";

/**
 * Dispatch a pointer-typed event carrying real clientX/Y (see JSDOM note).
 * Routed through `fireEvent` so the resulting React state update flushes
 * inside `act()`.
 */
function pointer(
  el: Element,
  type: "pointerdown" | "pointermove" | "pointerup" | "pointercancel",
  clientX: number,
  clientY: number,
): void {
  fireEvent(
    el,
    new MouseEvent(type, { bubbles: true, cancelable: true, clientX, clientY }),
  );
}

/** Two-member fixture; the gap renders between member id 7 (above) and 9. */
function makeMembers(): ComparisonMember[] {
  const base = (over: Partial<ComparisonMember>): ComparisonMember => ({
    id: 0,
    comparison_id: 100,
    exposure_id: 0,
    display_order: 0,
    band_height: 1,
    y_offset: 0,
    normalization: "qwindow",
    color_override: null,
    label_override: null,
    q_window_min: null,
    q_window_max: null,
    peak_display: null,
    snapshot: undefined,
    is_stale: false,
    created_by: null,
    created_at: null,
    ...over,
  });
  return [
    base({ id: 7, exposure_id: 70, display_order: 0 }),
    base({ id: 9, exposure_id: 90, display_order: 1 }),
  ];
}

const labels = new Map<number, string>([
  [7, "member-7"],
  [9, "member-9"],
]);

/**
 * Real overlay subscriber — the same shape `MultiTracePlot` mounts for each
 * band. Reads the active member-id from the shared context and reflects it
 * onto `data-active-band` when it matches.
 */
function MemberBandOverlay(props: { memberId: number }): JSX.Element {
  const { activeBandMemberId } = useActiveBand();
  return (
    <div
      data-testid="member-trace"
      data-member-id={String(props.memberId)}
      {...(activeBandMemberId === props.memberId
        ? { "data-active-band": String(props.memberId) }
        : {})}
    />
  );
}

function seedDraft(members: ComparisonMember[]): void {
  useAppState.setState({
    activeDraft: {
      ...emptyDraft(),
      members: members.map((m) => ({
        id: m.id,
        exposure_id: m.exposure_id,
        display_order: m.display_order,
        band_height: m.band_height,
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

beforeEach(() => {
  seedDraft(makeMembers());
});

describe("Inter-row resize gap — Compare UX E-4", () => {
  it("renders the gap between members in rest state", () => {
    const { getByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
      />,
    );
    const gap = getByTestId("member-resize-gap");
    expect(gap).toHaveAttribute("data-state", "rest");
    expect(gap.style.height).toBe("4px");
  });

  it("flips to hover state on pointer-enter and resets on pointer-leave", () => {
    const { getByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
      />,
    );
    const gap = getByTestId("member-resize-gap");
    fireEvent.pointerEnter(gap);
    expect(gap).toHaveAttribute("data-state", "hover");
    expect(gap.style.height).toBe("12px");
    fireEvent.pointerLeave(gap);
    expect(gap).toHaveAttribute("data-state", "rest");
    expect(gap.style.height).toBe("4px");
  });

  it("dispatches band-height update on drag (Δy → onResize)", () => {
    const onResize = vi.fn();
    const { getByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
        onResize={onResize}
      />,
    );
    const gap = getByTestId("member-resize-gap");
    pointer(gap, "pointerdown", 0, 100);
    pointer(gap, "pointermove", 0, 110);
    pointer(gap, "pointerup", 0, 110);
    expect(onResize).toHaveBeenCalled();
    const last = onResize.mock.calls.at(-1);
    expect(last?.[0]).toMatchObject({ dy: 10 });
  });

  it("drag dispatches the existing resizeBands band-height logic", () => {
    const { getByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
      />,
    );
    const gap = getByTestId("member-resize-gap");
    // Drag down 30px on a 300px panel → +0.1 normalized: above 1.0 → 1.1.
    pointer(gap, "pointerdown", 0, 100);
    pointer(gap, "pointermove", 0, 130);
    pointer(gap, "pointerup", 0, 130);
    const m = useAppState.getState().activeDraft!.members;
    expect(m[0]!.band_height).toBeCloseTo(1.1, 5);
    expect(m[1]!.band_height).toBeCloseTo(0.9, 5);
  });

  it("tints the band on hover via data-active-band on the plot overlay", () => {
    const { getByTestId } = render(
      <ActiveBandProvider>
        <MemberMetaGutter
          members={makeMembers()}
          panelHeight={300}
          mode="edit"
          displayLabelByMemberId={labels}
        />
        <MemberBandOverlay memberId={7} />
      </ActiveBandProvider>,
    );
    const gap = getByTestId("member-resize-gap");
    const overlay = getByTestId("member-trace");
    expect(overlay).not.toHaveAttribute("data-active-band");
    fireEvent.pointerEnter(gap);
    expect(overlay).toHaveAttribute("data-active-band", "7");
    fireEvent.pointerLeave(gap);
    expect(overlay).not.toHaveAttribute("data-active-band");
  });

  it("aborts the drag on pointercancel — clears tint, no resizeBands, back to rest", () => {
    const onResize = vi.fn();
    const { getByTestId } = render(
      <ActiveBandProvider>
        <MemberMetaGutter
          members={makeMembers()}
          panelHeight={300}
          mode="edit"
          displayLabelByMemberId={labels}
          onResize={onResize}
        />
        <MemberBandOverlay memberId={7} />
      </ActiveBandProvider>,
    );
    const gap = getByTestId("member-resize-gap");
    const overlay = getByTestId("member-trace");
    const beforeAbove = useAppState.getState().activeDraft!.members[0]!
      .band_height;
    const beforeBelow = useAppState.getState().activeDraft!.members[1]!
      .band_height;

    // Begin a drag, move, then have the browser cancel the gesture.
    pointer(gap, "pointerdown", 0, 100);
    expect(gap).toHaveAttribute("data-state", "hover");
    expect(overlay).toHaveAttribute("data-active-band", "7");
    pointer(gap, "pointermove", 0, 130);
    onResize.mockClear();
    pointer(gap, "pointercancel", 0, 130);

    // Tint cleared, gap reset to rest.
    expect(overlay).not.toHaveAttribute("data-active-band");
    expect(gap).toHaveAttribute("data-state", "rest");
    // A cancelled drag commits no resize and fires no final onResize.
    const m = useAppState.getState().activeDraft!.members;
    expect(m[0]!.band_height).toBeCloseTo(beforeAbove, 5);
    expect(m[1]!.band_height).toBeCloseTo(beforeBelow, 5);
    expect(onResize).not.toHaveBeenCalled();
  });

  it("clears the published tint when the gutter unmounts mid-drag", () => {
    function Harness(props: { showGutter: boolean }): JSX.Element {
      return (
        <ActiveBandProvider>
          {props.showGutter ? (
            <MemberMetaGutter
              members={makeMembers()}
              panelHeight={300}
              mode="edit"
              displayLabelByMemberId={labels}
            />
          ) : null}
          <MemberBandOverlay memberId={7} />
        </ActiveBandProvider>
      );
    }
    const { getByTestId, rerender } = render(<Harness showGutter />);
    const gap = getByTestId("member-resize-gap");
    const overlay = getByTestId("member-trace");

    // Begin a drag — the gap publishes its active band id.
    pointer(gap, "pointerdown", 0, 100);
    expect(overlay).toHaveAttribute("data-active-band", "7");

    // Unmount the gutter mid-drag (member removed / edit mode exits).
    rerender(<Harness showGutter={false} />);

    // The unmount-cleanup effect cleared the leaked tint.
    expect(overlay).not.toHaveAttribute("data-active-band");
  });
});
