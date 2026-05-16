/**
 * ActiveBandContext — minimal gutter↔plot coupling for the inter-row resize
 * gap (Compare UX E-4).
 *
 * `MemberMetaGutter` and `MultiTracePlot` are siblings under `Compare.tsx`
 * with no shared state. When a resize gap is hovered or dragged, the band
 * directly ABOVE the gap must tint accent on the plot side. Rather than
 * lift drag state into the page (it has no other use for it) this context
 * carries a single nullable member id: the band currently "active" for the
 * tint coupling.
 *
 * Contract:
 *   - The gutter is the sole publisher: it calls `setActiveBandMemberId`.
 *   - The plot overlay is the sole subscriber: each per-band overlay sets
 *     `data-active-band` when its `data-member-id` matches the active id.
 *   - Outside a provider both reads no-op (`null` / `() => {}`), so the
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
}

const ActiveBandContext = createContext<ActiveBandValue>({
  activeBandMemberId: null,
  setActiveBandMemberId: () => {},
});

export function ActiveBandProvider(props: {
  children: React.ReactNode;
}): JSX.Element {
  const [activeBandMemberId, setActiveBandMemberId] = useState<number | null>(
    null,
  );
  const value = useMemo<ActiveBandValue>(
    () => ({ activeBandMemberId, setActiveBandMemberId }),
    [activeBandMemberId],
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
