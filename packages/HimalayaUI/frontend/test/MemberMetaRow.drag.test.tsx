/**
 * MemberMetaRow grab-anywhere drag-vs-click tests — Compare UX Phase E, Task E-3.
 *
 * E-2 made collapse/expand a controlled prop toggled by clicking the row
 * body. E-3 makes the whole row "grab-anywhere" draggable: a pointer
 * gesture is a CLICK (→ toggle expand) when the pointer moved ≤4px before
 * release, or a DRAG (→ reorder) when it moved >4px.
 *
 * The 4px threshold + Manhattan-distance semantics come from
 * `makeDragThresholdState` in `lib/comparison/dragThreshold.ts` — that
 * helper is the source of truth (it uses `|dx| + |dy| > thresholdPx`).
 *
 * JSDOM note: JSDOM's `PointerEvent` constructor exists but does NOT
 * carry `clientX`/`clientY` from its init dict, so `fireEvent.pointerMove`
 * delivers `undefined` coordinates — the threshold math would see `NaN`.
 * `MouseEvent` *does* carry `clientX`/`clientY` in JSDOM, and React's
 * synthetic `onPointerDown`/`Move`/`Up` listen for the DOM event *type*
 * (`"pointerdown"` …), so dispatching a `MouseEvent` with a pointer-event
 * type triggers the React handler with real coordinates. This mirrors the
 * `new MouseEvent(...)` pattern already used in `peakTooltip.test.tsx`.
 */
import { useState } from "react";
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemberMetaRow } from "../src/components/MemberMetaRow";
import type { ComparisonMember } from "../src/api";

/**
 * Dispatch a pointer-typed event carrying real clientX/Y (see JSDOM note).
 * Routed through `fireEvent` so any resulting React state update flushes
 * inside `act()` (a bare `dispatchEvent` leaves the re-render pending and
 * trips the act(...) warning).
 */
function pointer(
  el: Element,
  type: "pointerdown" | "pointermove" | "pointerup",
  clientX: number,
  clientY: number,
): void {
  fireEvent(
    el,
    new MouseEvent(type, { bubbles: true, cancelable: true, clientX, clientY }),
  );
}

function makeMember(over: Partial<ComparisonMember> = {}): ComparisonMember {
  return {
    id: 7,
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
        id: 9, phase: "Pn3m", lattice_d: 12.345, r_squared: 0.992, ngc: -1.51,
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
 * Controlled-prop harness — the row never owns expand state (the gutter
 * does). Tests pass spies for `onToggleExpand` / `onDragStart` directly.
 */
function DragRow(props: {
  onToggleExpand: () => void;
  onDragStart: (memberId: number) => void;
  member?: ComparisonMember;
}): JSX.Element {
  const [expanded, setExpanded] = useState(false);
  return (
    <MemberMetaRow
      member={props.member ?? makeMember()}
      top={0}
      height={50}
      mode="edit"
      memberIndex={0}
      displayLabel="row-label"
      expanded={expanded}
      onToggleExpand={() => {
        setExpanded((e) => !e);
        props.onToggleExpand();
      }}
      onDragStart={props.onDragStart}
    />
  );
}

describe("MemberMetaRow drag-vs-click threshold — Compare UX E-3", () => {
  it("treats a pointer-up within 4px as toggle-expand", () => {
    const onToggleExpand = vi.fn();
    const onDragStart = vi.fn();
    render(<DragRow onToggleExpand={onToggleExpand} onDragStart={onDragStart} />);
    const body = screen.getByTestId("member-meta-row-body");
    pointer(body, "pointerdown", 10, 10);
    pointer(body, "pointermove", 12, 10);
    pointer(body, "pointerup", 12, 10);
    expect(onToggleExpand).toHaveBeenCalledTimes(1);
    expect(onDragStart).not.toHaveBeenCalled();
  });

  it("treats a pointer-move beyond 4px as drag (NOT toggle-expand)", () => {
    const onToggleExpand = vi.fn();
    const onDragStart = vi.fn();
    render(<DragRow onToggleExpand={onToggleExpand} onDragStart={onDragStart} />);
    const body = screen.getByTestId("member-meta-row-body");
    pointer(body, "pointerdown", 10, 10);
    pointer(body, "pointermove", 20, 10); // Δx=10>4
    pointer(body, "pointerup", 20, 10);
    expect(onDragStart).toHaveBeenCalledWith(7);
    expect(onToggleExpand).not.toHaveBeenCalled();
  });

  it("uses Manhattan distance — diagonal movement past threshold counts as drag", () => {
    // dragThreshold.ts uses Manhattan distance (|dx| + |dy| > thresholdPx),
    // so Δx=3 + Δy=3 = 6 > 4 → drag.
    const onToggleExpand = vi.fn();
    const onDragStart = vi.fn();
    render(<DragRow onToggleExpand={onToggleExpand} onDragStart={onDragStart} />);
    const body = screen.getByTestId("member-meta-row-body");
    pointer(body, "pointerdown", 10, 10);
    pointer(body, "pointermove", 13, 13); // Δx+Δy=6>4
    pointer(body, "pointerup", 13, 13);
    expect(onDragStart).toHaveBeenCalledTimes(1);
    expect(onDragStart).toHaveBeenCalledWith(7);
    expect(onToggleExpand).not.toHaveBeenCalled();
  });

  it("a drag detected only at pointer-up (no intervening move) still suppresses expand", () => {
    // The helper re-checks distance on pointer-up: a gesture with no move
    // event but a >4px down→up delta resolves to "drag", not "click".
    const onToggleExpand = vi.fn();
    const onDragStart = vi.fn();
    render(<DragRow onToggleExpand={onToggleExpand} onDragStart={onDragStart} />);
    const body = screen.getByTestId("member-meta-row-body");
    pointer(body, "pointerdown", 10, 10);
    pointer(body, "pointerup", 30, 10); // Δx=20>4
    expect(onToggleExpand).not.toHaveBeenCalled();
  });

  it("a pointer-up exactly at the down point (Δ=0) is a click — toggles expand", () => {
    const onToggleExpand = vi.fn();
    const onDragStart = vi.fn();
    render(<DragRow onToggleExpand={onToggleExpand} onDragStart={onDragStart} />);
    const body = screen.getByTestId("member-meta-row-body");
    pointer(body, "pointerdown", 50, 50);
    pointer(body, "pointerup", 50, 50);
    expect(onToggleExpand).toHaveBeenCalledTimes(1);
    expect(onDragStart).not.toHaveBeenCalled();
  });

  it("resets between gestures — a tap after a drag still resolves as a click", () => {
    // The threshold machine must clear its `dragging`/down-point state on
    // pointer-up; otherwise a >4px drag would leak into the next gesture
    // and a subsequent zero-displacement tap would mis-resolve as a drag.
    const onToggleExpand = vi.fn();
    const onDragStart = vi.fn();
    render(<DragRow onToggleExpand={onToggleExpand} onDragStart={onDragStart} />);
    const body = screen.getByTestId("member-meta-row-body");
    // Gesture 1 — a >4px drag.
    pointer(body, "pointerdown", 10, 10);
    pointer(body, "pointermove", 30, 10); // Δx=20>4
    pointer(body, "pointerup", 30, 10);
    expect(onDragStart).toHaveBeenCalledTimes(1);
    expect(onToggleExpand).not.toHaveBeenCalled();
    // Gesture 2 — a zero-displacement tap on the same row → CLICK.
    pointer(body, "pointerdown", 60, 60);
    pointer(body, "pointerup", 60, 60);
    expect(onToggleExpand).toHaveBeenCalledTimes(1);
    expect(onDragStart).toHaveBeenCalledTimes(1); // unchanged — no new drag
  });
});
