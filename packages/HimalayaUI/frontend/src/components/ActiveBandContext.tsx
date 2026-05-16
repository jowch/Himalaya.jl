/**
 * ActiveBandContext — minimal gutter↔plot coupling for the inter-row resize
 * gap (Compare UX E-4) and the reorder drop target (Compare UX E-5).
 *
 * `MemberMetaGutter` and `MultiTracePlot` are siblings under `Compare.tsx`
 * with no shared state. Two distinct gutter gestures need a plot-side
 * visual mirror:
 *
 *   - E-4 (resize): when a resize gap is hovered or dragged, the band
 *     directly ABOVE the gap tints accent. Carried by `activeBandMemberId`.
 *   - E-5 (reorder): while a row is being pointer-dragged to a new
 *     position, the band the row would drop INTO is highlighted as the
 *     drop target. Carried by `dropTargetMemberId`.
 *
 * The two fields are kept SEPARATE on purpose — they are distinct gestures
 * with independent lifecycles (a resize tint and a reorder drop target can
 * never be live at the same time, but conflating them into one field would
 * couple their cleanup paths). Rather than lift drag state into the page
 * (it has no other use for it) this context carries the two nullable
 * member ids.
 *
 * Contract:
 *   - The gutter is the sole publisher: it calls the setters.
 *   - The plot overlay is the sole subscriber: each per-band overlay sets
 *     `data-active-band` / `data-drop-target` when its `data-member-id`
 *     matches the corresponding active id.
 *   - Outside a provider all reads no-op (`null` / `() => {}`), so the
 *     gutter still renders standalone (review mode, isolated tests).
 *
 * Kept colocated + tiny on purpose — see components/AGENTS.md "minimal
 * context" guidance.
 */
import { createContext, useContext, useMemo, useState } from "react";

interface ActiveBandValue {
  /** Member id of the band currently hovered/dragged for tint, or null. */
  activeBandMemberId: number | null;
  /** Publisher — called by `MemberMetaGutter` only. */
  setActiveBandMemberId: (id: number | null) => void;
  /**
   * Compare UX E-5 — member id of the band a reorder drag would drop the
   * dragged row INTO, or null when no reorder drag is active. Distinct
   * from `activeBandMemberId`: that is the resize-tint gesture, this is the
   * reorder drop-target gesture.
   */
  dropTargetMemberId: number | null;
  /** Publisher — called by `MemberMetaGutter` only. */
  setDropTargetMemberId: (id: number | null) => void;
}

const ActiveBandContext = createContext<ActiveBandValue>({
  activeBandMemberId: null,
  setActiveBandMemberId: () => {},
  dropTargetMemberId: null,
  setDropTargetMemberId: () => {},
});

export function ActiveBandProvider(props: {
  children: React.ReactNode;
}): JSX.Element {
  const [activeBandMemberId, setActiveBandMemberId] = useState<number | null>(
    null,
  );
  const [dropTargetMemberId, setDropTargetMemberId] = useState<number | null>(
    null,
  );
  const value = useMemo<ActiveBandValue>(
    () => ({
      activeBandMemberId,
      setActiveBandMemberId,
      dropTargetMemberId,
      setDropTargetMemberId,
    }),
    [activeBandMemberId, dropTargetMemberId],
  );
  return (
    <ActiveBandContext.Provider value={value}>
      {props.children}
    </ActiveBandContext.Provider>
  );
}

export function useActiveBand(): ActiveBandValue {
  return useContext(ActiveBandContext);
}
