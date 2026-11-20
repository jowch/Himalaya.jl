import { useEffect } from "react";
import { useAppState } from "../state";
import { useExposures } from "../queries";
import type { Exposure } from "../api";
import { keptExposures, representativeExposure } from "../lib/exposure";

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
  const noUsable = keptExposures(exposures).length === 0;
  return { noUsable, allRejected: noUsable && exposures.length > 0 };
}

export function resolveActiveExposure(
  currentId: number | undefined,
  exposures: readonly Exposure[] | undefined,
): number | undefined {
  if (keptExposures(exposures).some((e) => e.id === currentId)) return currentId;
  return representativeExposure(exposures)?.id;
}

/**
 * Auto-select the right exposure for the active sample (focus workspace) and
 * PERSIST it to the store, so consumers that read `activeExposureId` directly
 * (the exposure stepper, deliberate switches, the loupe) stay coherent. The
 * selection rule lives in `resolveActiveExposure` — the page resolves
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
