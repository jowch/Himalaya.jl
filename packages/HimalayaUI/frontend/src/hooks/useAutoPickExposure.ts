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
 * Resolve WHICH exposure the focus workspace should show, given the currently
 * stored pick and the sample's (possibly still-loading) exposures. Pure, so the
 * page can resolve it DURING RENDER — not just in the auto-pick effect — which
 * is what lets a navigation onto a cached sample paint its trace in the same
 * frame instead of flashing a skeleton through the store→effect seam
 * (FO-NAV-SKELETON). The selection rule:
 *
 *   1. If the stored `currentId` is already a valid acceptable exposure on this
 *      sample, KEEP it. This is what lets the R5 (#228, F-11) representative
 *      switcher hold a deliberate in-session switch — a switch sets a valid
 *      exposure, and the resolver must not yank it back to the representative.
 *   2. Otherwise (cold arrival with undefined, a stale id left over from another
 *      sample, or the current exposure was rejected) adopt the flagged
 *      representative — the `selected` exposure set by
 *      `PATCH /api/exposures/:id/select`. This preserves the load-bearing
 *      Inspect→focus handoff: arriving with `currentId === undefined` lands on
 *      the marked-for-indexing exposure.
 *   3. Otherwise fall back to the first acceptable exposure.
 *   4. No acceptable exposures (none, or not yet loaded) → undefined. The page
 *      reads that as "still resolving" (skeleton) or "no exposures" (empty),
 *      per whether the exposures have actually loaded.
 */
export function resolveActiveExposure(
  currentId: number | undefined,
  exposures: readonly Exposure[] | undefined,
): number | undefined {
  const acceptable = acceptableExposures(exposures);
  if (acceptable.length === 0) return undefined;
  if (acceptable.some((e) => e.id === currentId)) return currentId;
  const flagged = acceptable.find((e) => e.selected);
  return (flagged ?? acceptable[0]!).id;
}

/**
 * Auto-select the right exposure for the active sample (focus workspace) and
 * PERSIST it to the store, so consumers that read `activeExposureId` directly
 * (the exposure stepper, deliberate switches, the loupe) stay coherent. The
 * selection rule lives in `resolveActiveExposure` (above) — the page resolves
 * the same value at render time, so this effect only reconciles the store.
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
    const next = resolveActiveExposure(activeExposureId, exposuresQ.data);
    // Only write when the rule picks a DIFFERENT exposure (an undefined result
    // means "nothing to adopt yet" — leave the store alone so we don't thrash).
    if (next !== undefined && next !== activeExposureId) setActiveExposure(next);
  }, [exposuresQ.data, activeExposureId, setActiveExposure]);
}
