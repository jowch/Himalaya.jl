/**
 * Drop indicator + plot mirror on reorder drag — Compare UX E-5.
 *
 * Reorder in this codebase COMMITS through the platform's HTML5
 * drag-and-drop path: `MemberMetaRow`'s body keeps its `draggable`
 * attribute (Task E-3, commit 563f048), the gutter's native
 * `onDragOver`/`onDrop` (`handleDrop`) compute the new order and dispatch
 * `reorderMembers`. The POINTER events from E-3 exist only for the
 * click-vs-drag 4px threshold; E-5's pointer-move tracking drives the 2px
 * `<div data-testid="drop-indicator">` and publishes the drop-target
 * member id to `ActiveBandContext` so the plot's per-band overlay mirrors
 * it with `data-drop-target="true"`.
 *
 * Real-browser gesture ordering — the regression E-5 originally shipped:
 *   pointerdown → pointermove(>4px) → native dragstart → pointercancel → drop
 * The browser fires `pointercancel` the instant a `draggable` element
 * begins a native drag. The fix makes `pointercancel` VISUAL-ONLY teardown
 * so it never nulls `dragSourceIdxRef`; the native `drop` stays the SINGLE
 * reorder-commit path and can still read the source index. These tests
 * drive that exact ordering and assert `reorderMembers` fires exactly once.
 *
 * JSDOM note: JSDOM's `PointerEvent` constructor does NOT carry
 * `clientX`/`clientY` from its init dict, so `fireEvent.pointerMove`
 * delivers `undefined` coordinates and the threshold math sees `NaN`.
 * `MouseEvent` *does* carry them, and React's synthetic pointer handlers
 * listen for the DOM event *type*, so dispatching a `MouseEvent` with a
 * pointer-event type triggers the handler with real coordinates. Same
 * `pointer()` helper as `MemberMetaRow.drag.test.tsx`.
 *
 * JSDOM has no layout engine, so the gutter cannot read row geometry from
 * `getBoundingClientRect` (it returns zero rects). The drop-position math
 * works off the `computeYBands` envelopes the component already derives
 * from `panelHeight` — those are plain numbers, layout-engine-free. The
 * test maps a pointer clientY directly onto that band space.
 *
 * RTL `fireEvent` DOES support synthetic `dragOver`/`drop` in JSDOM, so the
 * HTML5 commit path is fully testable.
 */
import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, fireEvent } from "@testing-library/react";
import { useAppState } from "../src/state";
import { emptyDraft } from "../src/lib/comparison/draft";
import type { SeriesMember } from "../src/api";
import { MemberMetaGutter } from "../src/components/MemberMetaGutter";
import {
  ActiveBandProvider,
  useActiveBand,
} from "../src/components/ActiveBandContext";

/** Dispatch a pointer-typed event carrying real clientX/Y (see JSDOM note). */
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

/**
 * Minimal `dataTransfer` stub — JSDOM's synthetic `dragOver`/`drop` events
 * carry no `dataTransfer`, so `handleDragOver`'s `e.dataTransfer.dropEffect`
 * write would throw. Same stub shape as `MemberReorder.test.tsx`.
 */
function makeDataTransfer(): Record<string, unknown> {
  return {
    types: ["text/plain"],
    data: {} as Record<string, string>,
    setData(this: { data: Record<string, string> }, k: string, v: string) {
      this.data[k] = v;
    },
    getData(this: { data: Record<string, string> }, k: string) {
      return this.data[k] ?? "";
    },
    effectAllowed: "move",
    dropEffect: "move",
  };
}

/** Three-member fixture — three equal 100px bands on a 300px panel. */
function makeMembers(): SeriesMember[] {
  const base = (over: Partial<SeriesMember>): SeriesMember => ({
    id: 0,
    series_id: 100,
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
    base({ id: 11, exposure_id: 110, display_order: 2 }),
  ];
}

const labels = new Map<number, string>([
  [7, "member-7"],
  [9, "member-9"],
  [11, "member-11"],
]);

/**
 * Real overlay subscriber — the same shape `MultiTracePlot` mounts for each
 * band. Reflects the published drop-target id onto `data-drop-target`.
 */
function MemberBandOverlay(props: { memberId: number }): JSX.Element {
  const { dropTargetMemberId } = useActiveBand();
  return (
    <div
      data-testid={`overlay-${props.memberId}`}
      data-member-id={String(props.memberId)}
      {...(dropTargetMemberId === props.memberId
        ? { "data-drop-target": "true" }
        : {})}
    />
  );
}

function seedDraft(members: SeriesMember[]): void {
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

/**
 * Begin a reorder drag on a row body: pointerdown then a >4px pointermove
 * crosses the threshold so `MemberMetaRow` calls the gutter's `onDragStart`,
 * which primes `dragSourceIdxRef`. Mirrors what a real browser does just
 * before it escalates the gesture to a native `dragstart`.
 */
function startReorderDrag(rowBody: Element, startY: number): void {
  pointer(rowBody, "pointerdown", 10, startY);
  pointer(rowBody, "pointermove", 30, startY); // Δx=20 > 4 → drag-start
}

describe("Drop indicator + plot mirror — Compare UX E-5", () => {
  it("renders a 2px drop indicator at the hover position during a reorder drag", () => {
    const { getByTestId, queryByTestId, getAllByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
      />,
    );
    // No drag in progress → no indicator.
    expect(queryByTestId("drop-indicator")).toBeNull();

    // Drag row 0 (member 7), then move the pointer down into band 2 (member
    // 11 occupies y∈[200,300]).
    const rowBody = getAllByTestId("member-meta-row-body")[0]!;
    startReorderDrag(rowBody, 50);
    pointer(rowBody, "pointermove", 30, 250);

    const indicator = getByTestId("drop-indicator");
    // 2px tall, inline-styled so JSDOM can read it without a layout pass.
    expect(indicator.style.height).toBe("2px");
    // Positioned at the inter-row boundary the pointer is over (bottom of
    // band 2 = 300px), not left at 0.
    expect(parseFloat(indicator.style.top)).toBeGreaterThan(0);
  });

  it("draws the indicator at the TOP edge when dropping ABOVE the source row", () => {
    // Drag row index 1 (member 9), hover into band 0 (member 7, y∈[0,100]).
    // Dropping above the source → indicator at the TARGET band's TOP edge
    // (band 0 top = 0; indicator `top` is `dropIndicatorY - 1` = -1px). The
    // edge decision reads the reactive `dragSourceIdx` STATE, not the ref.
    const { getAllByTestId, getByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
      />,
    );
    const rowBody = getAllByTestId("member-meta-row-body")[1]!;
    startReorderDrag(rowBody, 150);
    pointer(rowBody, "pointermove", 30, 50); // into band 0

    const indicator = getByTestId("drop-indicator");
    // Band 0 top edge = 0 → inline `top` = 0 - 1 = -1px.
    expect(parseFloat(indicator.style.top)).toBe(-1);
  });

  it("draws the indicator at the BOTTOM edge when dropping BELOW the source row", () => {
    // Drag row index 1 (member 9), hover into band 2 (member 11, y∈[200,300]).
    // Dropping below the source → indicator at the TARGET band's BOTTOM edge
    // (band 2 bottom = 300; indicator `top` = 300 - 1 = 299px).
    const { getAllByTestId, getByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
      />,
    );
    const rowBody = getAllByTestId("member-meta-row-body")[1]!;
    startReorderDrag(rowBody, 150);
    pointer(rowBody, "pointermove", 30, 250); // into band 2

    const indicator = getByTestId("drop-indicator");
    // Band 2 bottom edge = 300 → inline `top` = 300 - 1 = 299px.
    expect(parseFloat(indicator.style.top)).toBe(299);
  });

  it("plot band overlay mirrors the indicator with data-drop-target=true", () => {
    const { getAllByTestId, getByTestId } = render(
      <ActiveBandProvider>
        <MemberMetaGutter
          members={makeMembers()}
          panelHeight={300}
          mode="edit"
          displayLabelByMemberId={labels}
        />
        <MemberBandOverlay memberId={7} />
        <MemberBandOverlay memberId={9} />
        <MemberBandOverlay memberId={11} />
      </ActiveBandProvider>,
    );
    // Nothing mirrored before a drag.
    expect(getByTestId("overlay-9")).not.toHaveAttribute("data-drop-target");

    // Drag row 0 and hover over band index 1 (member 9, y∈[100,200]).
    const rowBody = getAllByTestId("member-meta-row-body")[0]!;
    startReorderDrag(rowBody, 50);
    pointer(rowBody, "pointermove", 30, 150);

    expect(getByTestId("overlay-9")).toHaveAttribute("data-drop-target", "true");
    // Only the hovered band mirrors.
    expect(getByTestId("overlay-7")).not.toHaveAttribute("data-drop-target");
    expect(getByTestId("overlay-11")).not.toHaveAttribute("data-drop-target");
  });

  it("commits the reorder via the native HTML5 drop — exactly once", () => {
    const reorderMembers = vi.fn();
    useAppState.setState({ reorderMembers });

    const { getAllByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
      />,
    );
    const rowBody = getAllByTestId("member-meta-row-body")[0]!;
    // Prime the drag source (member 7, index 0).
    startReorderDrag(rowBody, 50);
    // Native HTML5 path: dragOver then drop on the band-2 wrapper (index 2).
    const bandWrappers = getAllByTestId("member-meta-row").map(
      (row) => row.parentElement!,
    );
    const dt = makeDataTransfer();
    fireEvent.dragOver(bandWrappers[2]!, { dataTransfer: dt });
    fireEvent.drop(bandWrappers[2]!, { dataTransfer: dt });

    expect(reorderMembers).toHaveBeenCalledTimes(1);
    // Lift index 0 out, insert at index 2 → [1, 2, 0].
    expect(reorderMembers).toHaveBeenCalledWith([1, 2, 0]);
  });

  it("REGRESSION GUARD: pointercancel before the native drop must NOT abort the commit", () => {
    // Real-browser gesture ordering: pointerdown → pointermove(drag-start) →
    // native dragstart → pointercancel → dragOver → drop. The browser fires
    // `pointercancel` the instant a `draggable` element begins a native
    // drag. E-5 originally nulled `dragSourceIdxRef` in the pointercancel
    // handler, so the native `drop` saw `null` and the reorder silently
    // never committed. The fix makes `pointercancel` visual-only.
    const reorderMembers = vi.fn();
    useAppState.setState({ reorderMembers });

    const { getAllByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
      />,
    );
    const rowBody = getAllByTestId("member-meta-row-body")[0]!;
    const bandWrappers = getAllByTestId("member-meta-row").map(
      (row) => row.parentElement!,
    );

    startReorderDrag(rowBody, 50);
    // Browser fires pointercancel at native dragstart — BEFORE the drop.
    pointer(rowBody, "pointercancel", 30, 50);
    // Native drop still lands.
    const dt = makeDataTransfer();
    fireEvent.dragOver(bandWrappers[2]!, { dataTransfer: dt });
    fireEvent.drop(bandWrappers[2]!, { dataTransfer: dt });

    // The reorder STILL commits, exactly once.
    expect(reorderMembers).toHaveBeenCalledTimes(1);
    expect(reorderMembers).toHaveBeenCalledWith([1, 2, 0]);
  });

  it("clears the indicator and the mirror on pointercancel (visual-only teardown)", () => {
    const { getAllByTestId, getByTestId, queryByTestId } = render(
      <ActiveBandProvider>
        <MemberMetaGutter
          members={makeMembers()}
          panelHeight={300}
          mode="edit"
          displayLabelByMemberId={labels}
        />
        <MemberBandOverlay memberId={9} />
      </ActiveBandProvider>,
    );
    const rowBody = getAllByTestId("member-meta-row-body")[0]!;
    startReorderDrag(rowBody, 50);
    pointer(rowBody, "pointermove", 30, 150);
    expect(queryByTestId("drop-indicator")).not.toBeNull();
    expect(getByTestId("overlay-9")).toHaveAttribute("data-drop-target", "true");

    // pointercancel tears down the VISUALS (indicator + mirror) but does not
    // touch `dragSourceIdxRef` — the native drop still owns the commit.
    pointer(rowBody, "pointercancel", 30, 150);

    expect(queryByTestId("drop-indicator")).toBeNull();
    expect(getByTestId("overlay-9")).not.toHaveAttribute("data-drop-target");
  });

  it("clears the indicator and the mirror on native dragend", () => {
    const { getAllByTestId, getByTestId, queryByTestId } = render(
      <ActiveBandProvider>
        <MemberMetaGutter
          members={makeMembers()}
          panelHeight={300}
          mode="edit"
          displayLabelByMemberId={labels}
        />
        <MemberBandOverlay memberId={9} />
      </ActiveBandProvider>,
    );
    const rowBody = getAllByTestId("member-meta-row-body")[0]!;
    startReorderDrag(rowBody, 50);
    pointer(rowBody, "pointermove", 30, 150);
    expect(queryByTestId("drop-indicator")).not.toBeNull();
    expect(getByTestId("overlay-9")).toHaveAttribute("data-drop-target", "true");

    // A drag aborted off-target fires `dragend` with no `drop` — it clears
    // both `dragSourceIdxRef` and the visuals so nothing leaks.
    fireEvent.dragEnd(rowBody);

    expect(queryByTestId("drop-indicator")).toBeNull();
    expect(getByTestId("overlay-9")).not.toHaveAttribute("data-drop-target");
  });

  it("does not commit a reorder when the native drag ends with dragend and no drop", () => {
    const reorderMembers = vi.fn();
    useAppState.setState({ reorderMembers });

    const { getAllByTestId } = render(
      <MemberMetaGutter
        members={makeMembers()}
        panelHeight={300}
        mode="edit"
        displayLabelByMemberId={labels}
      />,
    );
    const rowBody = getAllByTestId("member-meta-row-body")[0]!;
    startReorderDrag(rowBody, 50);
    pointer(rowBody, "pointermove", 30, 150);
    // Drag aborted: dragend fires, no drop.
    fireEvent.dragEnd(rowBody);

    expect(reorderMembers).not.toHaveBeenCalled();
  });

  it("does not leak a stale drop-target into the context on mid-drag unmount", () => {
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
          <MemberBandOverlay memberId={9} />
        </ActiveBandProvider>
      );
    }
    const { getAllByTestId, getByTestId, rerender } = render(
      <Harness showGutter />,
    );
    const rowBody = getAllByTestId("member-meta-row-body")[0]!;
    startReorderDrag(rowBody, 50);
    pointer(rowBody, "pointermove", 30, 150);
    expect(getByTestId("overlay-9")).toHaveAttribute("data-drop-target", "true");

    // Unmount the gutter mid-drag (member removed / edit mode exits).
    rerender(<Harness showGutter={false} />);

    expect(getByTestId("overlay-9")).not.toHaveAttribute("data-drop-target");
  });
});
