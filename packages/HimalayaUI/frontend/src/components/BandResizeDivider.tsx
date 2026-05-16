/**
 * BandResizeGap — visible inter-row resize handle between two adjacent
 * metadata rows in edit mode (Compare UX E-4; supersedes the Phase 7.2
 * 1px hairline `BandResizeDivider`).
 *
 * Visible affordance:
 *   - 4px tall at rest, 12px tall on hover — driven by INLINE `style.height`
 *     ("4px" / "12px") so JSDOM tests can assert without a layout pass and
 *     so the height never depends on a Tailwind class string.
 *   - Component-controlled `data-state` = "rest" | "hover".
 *   - `cursor: row-resize` while hovering.
 *   - On hover OR active drag, publishes the member id ABOVE the gap to
 *     `ActiveBandContext` so the plot's per-band overlay tints accent.
 *
 * Push semantics in normalized space (unchanged from Phase 7.2):
 *   above grows by Δ, below shrinks by Δ → Σ band_heights stays constant.
 *   Clamp at 0.1 (handled in `resizeBands`).
 *   Snap to neighbor: when the drag puts the boundary within a 5px dead
 *   zone of the even-split delta, snap to it.
 *
 * Drag math:
 *   - pointerDown captures startY.
 *   - pointerMove tracks the cumulative dy = currentY - startY.
 *   - pointerUp invokes `resizeBands(memberIdx, snappedDy, totalHeightPx)`
 *     AND forwards `{ dy }` to the optional `onResize` callback (E-4 test
 *     hook + parent observability).
 *   - pointerCancel (touch interruption, OS gesture) ABORTS the drag: it
 *     resets `dragging`, clears `startYRef`, releases capture, and clears
 *     the tint — but dispatches NO `resizeBands` mutation and fires NO
 *     final `onResize` (a cancelled drag is not a committed resize).
 *   - Unmounting mid-drag clears the published tint iff this instance is
 *     the active publisher (`ownsActiveBandRef`).
 *
 * Test selector contract:
 *   `data-testid="member-resize-gap"`, `data-state`, `data-above-id`,
 *   `data-below-id`.
 *
 * The legacy `BandResizeDivider` (Phase 7.2 hairline,
 * `data-testid="band-divider"`, mouse-driven) is retained at the bottom of
 * this file. It is NOT a wrapper over `BandResizeGap` — it is a full
 * retained copy of the original Phase 7.2 implementation, kept only so
 * `test/BandResizeDivider.test.tsx` stays green. Slated for deletion once
 * that test is migrated to the gap surface. New call sites use
 * `BandResizeGap`; the gutter no longer mounts the hairline.
 */
import { useCallback, useEffect, useRef, useState } from "react";
import { useAppState } from "../state";
import { useActiveBand } from "./ActiveBandContext";

const SNAP_DEAD_ZONE_PX = 5;

/** Resize-gap heights (px). Inline so JSDOM tests can read them directly. */
const GAP_HEIGHT_REST_PX = "4px";
const GAP_HEIGHT_HOVER_PX = "12px";

export interface BandResizeGapProps {
  /** Member id above the gap (its band shrinks when the gap drags up). */
  aboveId: number;
  /** Member id below the gap. */
  belowId: number;
  /** Index of the *above* member in `activeDraft.members`. */
  memberIndex: number;
  /** Pixel y-position of the band boundary in the gutter. */
  top: number;
  /** Total panel height in pixels — denominator for normalized deltas. */
  totalHeightPx: number;
  /**
   * Optional observer of the cumulative pointer delta. Receives `{ dy }`
   * on every pointerMove and a final call on pointerUp carrying the total
   * dy. Genuinely optional (parents that only need the Zustand side
   * effect omit it) — see `exactOptionalPropertyTypes`.
   */
  onResize?: (delta: { dy: number }) => void;
}

export function BandResizeGap(props: BandResizeGapProps): JSX.Element {
  const { aboveId, belowId, memberIndex, top, totalHeightPx, onResize } = props;
  const startYRef = useRef<number | null>(null);
  const [dragging, setDragging] = useState(false);
  const [hovering, setHovering] = useState(false);
  const resizeBands = useAppState((s) => s.resizeBands);
  const { setActiveBandMemberId } = useActiveBand();

  // Tracks whether THIS instance currently owns the published active band
  // id. The unmount-cleanup effect below reads it to clear the tint only
  // when this gap is the live publisher — never stomping another gap's id.
  const ownsActiveBandRef = useRef(false);

  // Active = hovered OR mid-drag. Both states tint the band above and grow
  // the gap to its 12px hit target.
  const active = hovering || dragging;

  // Publish (or clear) the active band id and keep `ownsActiveBandRef` in
  // sync so unmount cleanup knows whether this instance is the publisher.
  const publishActiveBand = useCallback(
    (id: number | null) => {
      ownsActiveBandRef.current = id !== null;
      setActiveBandMemberId(id);
    },
    [setActiveBandMemberId],
  );

  const handlePointerEnter = useCallback(() => {
    setHovering(true);
    publishActiveBand(aboveId);
  }, [aboveId, publishActiveBand]);

  const handlePointerLeave = useCallback(() => {
    setHovering(false);
    // Don't clear the tint mid-drag — the pointer can leave the 4/12px
    // strip while still dragging; pointerUp owns the clear in that case.
    if (!dragging) publishActiveBand(null);
  }, [dragging, publishActiveBand]);

  const handlePointerDown = useCallback(
    (e: React.PointerEvent) => {
      e.preventDefault();
      e.stopPropagation();
      startYRef.current = e.clientY;
      setDragging(true);
      publishActiveBand(aboveId);
      // Capture the pointer so move/up keep firing on this element even
      // when the cursor leaves the thin strip mid-drag.
      try {
        e.currentTarget.setPointerCapture(e.pointerId);
      } catch {
        // JSDOM may not implement setPointerCapture — non-fatal; the
        // element still receives move/up while the button is held.
      }
    },
    [aboveId, publishActiveBand],
  );

  const handlePointerMove = useCallback(
    (e: React.PointerEvent) => {
      if (startYRef.current === null) return;
      const dy = e.clientY - startYRef.current;
      onResize?.({ dy });
    },
    [onResize],
  );

  const handlePointerUp = useCallback(
    (e: React.PointerEvent) => {
      const startY = startYRef.current;
      startYRef.current = null;
      setDragging(false);
      try {
        e.currentTarget.releasePointerCapture(e.pointerId);
      } catch {
        // Non-fatal — see handlePointerDown.
      }
      if (startY === null) return;

      const rawDy = e.clientY - startY;

      // Snap target: the dy that equalizes the two neighbor bands.
      // even_above = even_below = (above + below) / 2.
      // Required normalized delta = (below - above) / 2.
      // In pixel space: snapTargetPx = ((below - above) / 2) * totalHeightPx.
      const draft = useAppState.getState().activeDraft;
      let dy = rawDy;
      if (draft && draft.members[memberIndex] && draft.members[memberIndex + 1]) {
        const above = draft.members[memberIndex]!.band_height;
        const below = draft.members[memberIndex + 1]!.band_height;
        const snapTargetPx = ((below - above) / 2) * totalHeightPx;
        if (Math.abs(rawDy - snapTargetPx) <= SNAP_DEAD_ZONE_PX) {
          dy = snapTargetPx;
        }
      }

      resizeBands(memberIndex, dy, totalHeightPx);
      // Final call delivers the total (possibly snapped) dy.
      onResize?.({ dy });
      // Drag finished — clear the tint unless the pointer is still hovering.
      if (!hovering) publishActiveBand(null);
    },
    [
      memberIndex,
      totalHeightPx,
      resizeBands,
      onResize,
      hovering,
      publishActiveBand,
    ],
  );

  // pointercancel — the browser aborted the gesture (touch interruption, OS
  // gesture). pointerUp will never arrive, so reset the drag state here.
  // A cancelled drag is NOT a committed resize: dispatch no `resizeBands`
  // mutation and fire no final `onResize`.
  const handlePointerCancel = useCallback(
    (e: React.PointerEvent) => {
      startYRef.current = null;
      setDragging(false);
      try {
        e.currentTarget.releasePointerCapture(e.pointerId);
      } catch {
        // Non-fatal — see handlePointerDown.
      }
      // Clear the tint unless the pointer is still hovering the strip.
      if (!hovering) publishActiveBand(null);
    },
    [hovering, publishActiveBand],
  );

  // Unmount cleanup — if this gap is torn down mid-drag (member removed,
  // edit mode exits) pointerUp/Cancel never fire, so the published tint
  // would leak. Clear it, but ONLY when this instance is the live
  // publisher, so we never stomp a sibling gap's active id.
  useEffect(() => {
    return () => {
      if (ownsActiveBandRef.current) setActiveBandMemberId(null);
    };
  }, [setActiveBandMemberId]);

  return (
    <div
      data-testid="member-resize-gap"
      data-state={active ? "hover" : "rest"}
      data-above-id={String(aboveId)}
      data-below-id={String(belowId)}
      onPointerEnter={handlePointerEnter}
      onPointerLeave={handlePointerLeave}
      onPointerDown={handlePointerDown}
      onPointerMove={handlePointerMove}
      onPointerUp={handlePointerUp}
      onPointerCancel={handlePointerCancel}
      style={{
        position: "absolute",
        left: 0,
        right: 0,
        // Center the strip on the band boundary; the visible height grows
        // downward+upward symmetrically when hovered.
        top: `${top - (active ? 6 : 2)}px`,
        height: active ? GAP_HEIGHT_HOVER_PX : GAP_HEIGHT_REST_PX,
        cursor: "row-resize",
      }}
      className={[
        "z-10 rounded-full transition-[height,background-color]",
        active ? "bg-accent/40" : "bg-border/60",
      ].join(" ")}
      aria-label="Resize band"
      role="separator"
    />
  );
}

// ── Legacy hairline divider (Phase 7.2) ─────────────────────────────────────
//
// Retained only so `test/BandResizeDivider.test.tsx` stays green through the
// E-4 transition. Not mounted by `MemberMetaGutter` anymore — `BandResizeGap`
// above is the live component. Safe to delete once that legacy test is
// migrated to the gap surface.

export interface BandResizeDividerProps {
  /** Member id above the divider (its band shrinks when divider drags up). */
  aboveId: number;
  /** Member id below the divider. */
  belowId: number;
  /** Index of the *above* member in `activeDraft.members`. */
  memberIndex: number;
  /** Pixel y-position of the divider in the gutter (band boundary). */
  top: number;
  /** Total panel height in pixels — denominator for normalized deltas. */
  totalHeightPx: number;
}

export function BandResizeDivider(props: BandResizeDividerProps): JSX.Element {
  const { aboveId, belowId, memberIndex, top, totalHeightPx } = props;
  const startYRef = useRef<number | null>(null);
  const [dragging, setDragging] = useState(false);
  const resizeBands = useAppState((s) => s.resizeBands);

  const handleMouseDown = useCallback((e: React.MouseEvent) => {
    e.preventDefault();
    e.stopPropagation();
    startYRef.current = e.clientY;
    setDragging(true);
  }, []);

  useEffect(() => {
    if (!dragging) return;

    const onMove = (e: MouseEvent) => {
      void e;
    };

    const onUp = (e: MouseEvent) => {
      const startY = startYRef.current;
      startYRef.current = null;
      setDragging(false);
      if (startY === null) return;

      const rawDeltaPx = e.clientY - startY;
      const draft = useAppState.getState().activeDraft;
      let deltaPx = rawDeltaPx;
      if (draft && draft.members[memberIndex] && draft.members[memberIndex + 1]) {
        const above = draft.members[memberIndex]!.band_height;
        const below = draft.members[memberIndex + 1]!.band_height;
        const snapTargetPx = ((below - above) / 2) * totalHeightPx;
        if (Math.abs(rawDeltaPx - snapTargetPx) <= SNAP_DEAD_ZONE_PX) {
          deltaPx = snapTargetPx;
        }
      }

      resizeBands(memberIndex, deltaPx, totalHeightPx);
    };

    window.addEventListener("mousemove", onMove);
    window.addEventListener("mouseup", onUp);
    return () => {
      window.removeEventListener("mousemove", onMove);
      window.removeEventListener("mouseup", onUp);
    };
  }, [dragging, memberIndex, totalHeightPx, resizeBands]);

  return (
    <div
      data-testid="band-divider"
      data-above-id={String(aboveId)}
      data-below-id={String(belowId)}
      onMouseDown={handleMouseDown}
      style={{
        position: "absolute",
        left: 0,
        right: 0,
        top: `${top - 3}px`,
        height: "6px",
        cursor: "row-resize",
      }}
      className={["z-10", dragging ? "bg-accent/40" : "hover:bg-accent/20"].join(
        " ",
      )}
      aria-label="Resize band"
      role="separator"
    />
  );
}
