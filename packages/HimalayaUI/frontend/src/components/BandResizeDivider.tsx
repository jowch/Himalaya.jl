/**
 * BandResizeDivider — drag handle between two adjacent metadata rows in
 * edit mode (Plan §Phase 7, Task 7.2).
 *
 * Push semantics in normalized space:
 *   above grows by Δ, below shrinks by Δ → Σ band_heights stays constant.
 * Clamp at 0.1 (handled in `resizeBands`).
 * Snap to neighbor: when the drag puts the divider within a 5px dead zone
 * of "match the band on the other side," snap to the even-split delta.
 *
 * Drag math:
 *   - mousedown captures startY + the divider's `top` (the parent supplies
 *     it from `computeYBands`).
 *   - mousemove tracks deltaPx = currentY - startY.
 *   - mouseup invokes `resizeBands(memberIdx, snappedDeltaPx, totalHeightPx)`.
 *
 * The dead-zone snap fires when the chosen `deltaPx` would land within 5px
 * of the even-split target. The even-split target is the `deltaPx` that
 * would equalize the two neighbor bands' ratios.
 *
 * Test selector contract:
 *   `data-testid="band-divider"`, `data-above-id`, `data-below-id` — used
 *   by Playwright + vitest tests per spec selector table.
 */
import { useCallback, useEffect, useRef, useState } from "react";
import { useAppState } from "../state";

const SNAP_DEAD_ZONE_PX = 5;

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

  // Read the current pair of band heights at drag-start so we can compute
  // the even-split snap target. We read from Zustand directly inside the
  // mousedown handler so the snap math always reflects live draft state.
  const handleMouseDown = useCallback(
    (e: React.MouseEvent) => {
      e.preventDefault();
      e.stopPropagation();
      startYRef.current = e.clientY;
      setDragging(true);
    },
    [],
  );

  useEffect(() => {
    if (!dragging) return;

    const onMove = (e: MouseEvent) => {
      // Visual feedback only — actual band update happens on mouseup so
      // intermediate moves don't dispatch a rerender storm into Zustand.
      void e;
    };

    const onUp = (e: MouseEvent) => {
      const startY = startYRef.current;
      startYRef.current = null;
      setDragging(false);
      if (startY === null) return;

      const rawDeltaPx = e.clientY - startY;

      // Compute snap target: deltaPx that equalizes the two neighbor bands.
      // even_above = even_below = (above + below) / 2.
      // Required normalized delta = (below - above) / 2.
      // In pixel space: snapTargetPx = ((below - above) / 2) * totalHeightPx.
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
      className={[
        "z-10",
        dragging ? "bg-accent/40" : "hover:bg-accent/20",
      ].join(" ")}
      aria-label="Resize band"
      role="separator"
    />
  );
}
