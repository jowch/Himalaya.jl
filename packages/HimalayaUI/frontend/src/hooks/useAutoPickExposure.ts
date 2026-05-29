import { useEffect } from "react";
import { useAppState } from "../state";
import { useExposures } from "../queries";
import type { Exposure } from "../api";

/**
 * Exposures that can be indexed — anything not explicitly rejected. The query
 * returns the unfiltered list (one cache row shared with Inspect); the focus
 * page is the only consumer that wants acceptable-only. Single source of truth
 * for "which exposures count", reused by the auto-pick effect below and by
 * PlotCard's zero-exposure empty state.
 */
export function acceptableExposures(
  exposures: readonly Exposure[] | undefined,
): Exposure[] {
  return (exposures ?? []).filter((e) => e.status !== "rejected");
}

/**
 * Whether the active sample has no usable exposure to index, and if so why.
 *
 * `noUsable` stays false until the exposures have actually loaded (undefined
 * data) so a caller doesn't flash an empty state during the fetch. When there
 * are zero acceptable exposures, `allRejected` distinguishes "the sample has
 * exposures but all are rejected" (true → offer to restore in the loupe) from
 * "the sample has no exposures at all" (false).
 */
export function noUsableExposureState(
  exposures: readonly Exposure[] | undefined,
): { noUsable: boolean; allRejected: boolean } {
  if (exposures === undefined) return { noUsable: false, allRejected: false };
  const noUsable = acceptableExposures(exposures).length === 0;
  return { noUsable, allRejected: noUsable && exposures.length > 0 };
}

/**
 * Auto-select the right exposure for the active sample (focus workspace):
 *
 *   1. If the current `activeExposureId` is already a valid acceptable
 *      exposure on this sample, KEEP it. This is what lets the R5 (#228, F-11)
 *      representative-exposure switcher hold a deliberate in-session switch —
 *      a switch sets a valid exposure, and the auto-pick must not yank it back
 *      to the representative on the next render.
 *   2. Otherwise (cold arrival with undefined, a stale id left over from
 *      another sample, or the current exposure was rejected) adopt the flagged
 *      representative — the `selected` exposure set by
 *      `PATCH /api/exposures/:id/select`. This preserves the load-bearing
 *      Inspect→focus handoff: arriving with `activeExposureId === undefined`
 *      lands on the marked-for-indexing exposure.
 *   3. Otherwise fall back to the first acceptable exposure.
 *
 * Extracted from PlotCard so the regression tests subscribe to the same effect
 * the production tree runs (`test/useAutoPickExposure.test.tsx`), eliminating
 * drift between a hand-maintained mirror and the real auto-pick logic.
 */
export function useAutoPickExposure(activeSampleId: number | undefined): void {
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const setActiveExposure = useAppState((s) => s.setActiveExposure);
  const exposuresQ = useExposures(activeSampleId);
  useEffect(() => {
    // Acceptable-only (see acceptableExposures): the shared cache row is the
    // unfiltered Inspect list, but the focus auto-pick wants indexable ones.
    const exposures = acceptableExposures(exposuresQ.data);
    if (exposures.length === 0) return;
    // 1. A valid current pick wins — never clobber a deliberate switch (F-11).
    if (exposures.some((e) => e.id === activeExposureId)) return;
    // 2. No valid current pick: adopt the flagged representative if present
    //    (the Inspect→focus marked-for-indexing handoff), else
    // 3. the first acceptable exposure.
    const flagged = exposures.find((e) => e.selected);
    setActiveExposure((flagged ?? exposures[0]!).id);
  }, [exposuresQ.data, activeExposureId, setActiveExposure]);
}
