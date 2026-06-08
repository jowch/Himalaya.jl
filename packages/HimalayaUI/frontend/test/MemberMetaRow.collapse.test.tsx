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
 *
 * Compare UX E-3 — the row body's expand affordance moved off a plain
 * `onClick` onto a pointerdown→up gesture gated by the 4px drag
 * threshold. A zero-displacement gesture (down + up at the same coords)
 * resolves to "click" → toggle expand. `tapBody` reproduces that here;
 * a bare `fireEvent.click` no longer dispatches the pointer events the
 * handler now listens for.
 */
import { useState } from "react";
import { describe, it, expect } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { MemberMetaRow } from "../src/components/MemberMetaRow";
import type { SeriesMember } from "../src/api";

/**
 * Compare UX E-3 — a zero-displacement pointer gesture on the row body.
 * JSDOM's `PointerEvent` drops `clientX`/`clientY`, but `MouseEvent`
 * carries them and React's `onPointerDown/Up` listen by event *type*, so
 * a `MouseEvent` typed `"pointerdown"`/`"pointerup"` triggers the handler
 * with real coords. Down == up coords → Manhattan distance 0 → "click".
 * Dispatched via `fireEvent` so the resulting React state update flushes
 * inside `act()` (a bare `dispatchEvent` would leave the re-render pending).
 */
function tapBody(el: Element): void {
  fireEvent(el, new MouseEvent("pointerdown", { bubbles: true, clientX: 5, clientY: 5 }));
  fireEvent(el, new MouseEvent("pointerup", { bubbles: true, clientX: 5, clientY: 5 }));
}

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
function ControlledRow(props: { member?: SeriesMember; mode?: "review" | "edit" }) {
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
    tapBody(screen.getByTestId("member-meta-row-body"));
    expect(screen.getByTestId("member-meta-row")).toHaveAttribute(
      "data-expanded",
      "true",
    );
  });

  it("collapses when clicked again", () => {
    render(<ControlledRow />);
    tapBody(screen.getByTestId("member-meta-row-body"));
    tapBody(screen.getByTestId("member-meta-row-body"));
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
    tapBody(screen.getByTestId("member-meta-row-body"));
    expect(screen.getByTestId("member-meta-label-input")).toBeInTheDocument();
  });

  it("carries data-member-id and data-interactable on the root", () => {
    render(<ControlledRow />);
    const row = screen.getByTestId("member-meta-row");
    expect(row).toHaveAttribute("data-member-id", "1");
    expect(row).toHaveAttribute("data-interactable", "expand");
  });
});
